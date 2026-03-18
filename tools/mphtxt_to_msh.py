#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path
from typing import Dict, List, Tuple


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip() if "#" in line else line.strip()


def next_nonempty(lines: List[str], idx: int) -> Tuple[str, int]:
    while idx < len(lines):
        s = strip_comment(lines[idx])
        if s:
            return s, idx + 1
        idx += 1
    raise ValueError("Unexpected EOF")


def parse_mphtxt(path: Path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    vertex_count = None
    for i, line in enumerate(lines):
        if "number of mesh vertices" in line:
            val = strip_comment(line)
            vertex_count = int(val) if val else int(next_nonempty(lines, i + 1)[0])
            break
    if vertex_count is None:
        raise ValueError("Cannot find number of mesh vertices")

    coord_start = None
    for i, line in enumerate(lines):
        if "# Mesh vertex coordinates" in line:
            coord_start = i + 1
            break
    if coord_start is None:
        raise ValueError("Cannot find vertex coordinates section")

    points: List[Tuple[float, float, float]] = []
    idx = coord_start
    while len(points) < vertex_count:
        s = strip_comment(lines[idx])
        if s:
            x, y, z = s.split()
            points.append((float(x), float(y), float(z)))
        idx += 1

    sections: Dict[str, Dict[str, List]] = {}
    while idx < len(lines):
        if "# Type #" not in lines[idx]:
            idx += 1
            continue

        idx += 1
        type_line, idx = next_nonempty(lines, idx)
        tokens = type_line.split()
        type_name = tokens[1]

        nverts = None
        nelems = None
        while idx < len(lines):
            raw = lines[idx]
            s = strip_comment(raw)
            idx += 1
            if not s:
                continue
            if "number of vertices per element" in raw:
                nverts = int(s)
            elif "number of elements" in raw:
                nelems = int(s)
                break
        if nverts is None or nelems is None:
            raise ValueError(f"Malformed section {type_name}")

        while idx < len(lines) and "# Elements" not in lines[idx]:
            idx += 1
        idx += 1

        elems: List[Tuple[int, ...]] = []
        for _ in range(nelems):
            s, idx = next_nonempty(lines, idx)
            conn = tuple(int(v) for v in s.split())
            if len(conn) != nverts:
                raise ValueError(f"Unexpected node count in {type_name}")
            elems.append(conn)

        while idx < len(lines) and "# Geometric entity indices" not in lines[idx]:
            idx += 1
        idx += 1

        geo: List[int] = []
        for _ in range(nelems):
            s, idx = next_nonempty(lines, idx)
            geo.append(int(s.split()[0]))

        sections[type_name] = {"elements": elems, "geo": geo}

    return points, sections


def write_msh(out_path: Path, points, sections, scale: float):
    elements: List[Tuple[int, int, int, Tuple[int, ...]]] = []
    nonpositive_tag_count = 0
    shifted_boundary_tag_count = 0
    skipped_interior_boundary_faces = 0

    def sq_dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        return dx * dx + dy * dy + dz * dz

    def vec_sub(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def vec_cross(u: Tuple[float, float, float], v: Tuple[float, float, float]) -> Tuple[float, float, float]:
        return (
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        )

    def vec_dot(u: Tuple[float, float, float], v: Tuple[float, float, float]) -> float:
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    def midpoint(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
        return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5)

    def signed_tet_volume6(c0: int, c1: int, c2: int, c3: int) -> float:
        p0, p1, p2, p3 = points[c0], points[c1], points[c2], points[c3]
        return vec_dot(vec_cross(vec_sub(p1, p0), vec_sub(p2, p0)), vec_sub(p3, p0))

    def signed_prism_volume6(c0: int, c1: int, c2: int, c3: int, c4: int, c5: int) -> float:
        # Decompose prism into 3 tetrahedra with Gmsh prism corner order.
        return (
            signed_tet_volume6(c0, c1, c2, c3)
            + signed_tet_volume6(c1, c2, c4, c3)
            + signed_tet_volume6(c2, c4, c5, c3)
        )

    def reorder_edges_by_geometry(corners: Tuple[int, ...], mids: Tuple[int, ...],
                                  edge_pairs: Tuple[Tuple[int, int], ...]) -> Tuple[int, ...]:
        corner_pts = [points[idx] for idx in corners]
        mid_pts = [points[idx] for idx in mids]
        target_midpoints = [midpoint(corner_pts[a], corner_pts[b]) for a, b in edge_pairs]

        best_perm = None
        best_err = math.inf
        for perm in itertools.permutations(range(len(mids))):
            err = 0.0
            for gmsh_slot, mid_slot in enumerate(perm):
                err += sq_dist(mid_pts[mid_slot], target_midpoints[gmsh_slot])
            if err < best_err:
                best_err = err
                best_perm = perm

        assert best_perm is not None
        return tuple(mids[i] for i in best_perm)

    def reorder_tri2_from_comsol(raw_conn: Tuple[int, ...], corners: Tuple[int, int, int]) -> Tuple[int, ...]:
        raw_corners = raw_conn[:3]
        raw_edge_to_mid = {
            (0, 1): raw_conn[3],
            (0, 2): raw_conn[4],
            (1, 2): raw_conn[5],
        }
        corner_to_raw = [raw_corners.index(c) for c in corners]

        gmsh_edges = ((0, 1), (1, 2), (2, 0))
        mids = []
        for a, b in gmsh_edges:
            ra = corner_to_raw[a]
            rb = corner_to_raw[b]
            mids.append(raw_edge_to_mid[tuple(sorted((ra, rb)))])
        return tuple(corners) + tuple(mids)

    def reorder_tet2_from_comsol(raw_conn: Tuple[int, ...],
                                 corners: Tuple[int, int, int, int]) -> Tuple[int, ...]:
        raw_corners = raw_conn[:4]
        raw_edge_to_mid = {
            (0, 1): raw_conn[4],
            (0, 2): raw_conn[5],
            (1, 2): raw_conn[6],
            (0, 3): raw_conn[7],
            (1, 3): raw_conn[8],
            (2, 3): raw_conn[9],
        }
        corner_to_raw = [raw_corners.index(c) for c in corners]

        gmsh_edges = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
        mids = []
        for a, b in gmsh_edges:
            ra = corner_to_raw[a]
            rb = corner_to_raw[b]
            mids.append(raw_edge_to_mid[tuple(sorted((ra, rb)))])
        return tuple(corners) + tuple(mids)

    def reorder_quadratic_nodes(name: str, conn: Tuple[int, ...]) -> Tuple[int, ...]:
        if name == "tri2":
            corners = conn[:3]
            return reorder_tri2_from_comsol(conn, corners)
        if name == "tet2":
            c0, c1, c2, c3 = conn[:4]
            if signed_tet_volume6(c0, c1, c2, c3) < 0.0:
                corners = (c1, c0, c2, c3)
            else:
                corners = (c0, c1, c2, c3)
            return reorder_tet2_from_comsol(conn, corners)
        return conn

    element_type_map = {
        "tri": (2, 3),
        "quad": (3, 4),
        "tri2": (9, 6),
        "tet": (4, 4),
        "hex": (5, 8),
        "tet2": (11, 10),
        "prism": (6, 6),
    }

    element_dim_map = {
        "tri": 2,
        "quad": 2,
        "tri2": 2,
        "tet": 3,
        "hex": 3,
        "tet2": 3,
        "prism": 3,
    }

    max_dim = 0
    for name in element_type_map:
        if name in sections and sections[name].get("elements"):
            max_dim = max(max_dim, element_dim_map[name])

    # For 3D meshes, enforce globally consistent tet orientation on shared faces first.
    oriented_tet_blocks: Dict[str, List[Tuple[int, ...]]] = {}
    oriented_prism_block: List[Tuple[int, ...]] = []
    oriented_hex_block: List[Tuple[int, ...]] = []
    if max_dim == 3:
        tet_entries: List[Tuple[str, int, Tuple[int, ...], Tuple[int, int, int, int]]] = []
        for tet_name in ("tet", "tet2"):
            if tet_name not in sections:
                continue
            for local_idx, conn in enumerate(sections[tet_name]["elements"]):
                c0, c1, c2, c3 = conn[:4]
                if signed_tet_volume6(c0, c1, c2, c3) < 0.0:
                    base_corners = (c1, c0, c2, c3)
                else:
                    base_corners = (c0, c1, c2, c3)
                tet_entries.append((tet_name, local_idx, conn, base_corners))

        def face_corners(c: Tuple[int, int, int, int], face_id: int) -> Tuple[int, int, int]:
            # MFEM Geometry::Constants<Geometry::TETRAHEDRON>::FaceVert
            # {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}}
            c0, c1, c2, c3 = c
            if face_id == 0:
                return (c1, c2, c3)
            if face_id == 1:
                return (c0, c3, c2)
            if face_id == 2:
                return (c0, c1, c3)
            return (c0, c2, c1)

        def rotate3(t: Tuple[int, int, int], k: int) -> Tuple[int, int, int]:
            if k == 0:
                return t
            if k == 1:
                return (t[1], t[2], t[0])
            return (t[2], t[0], t[1])

        def same_oriented(a: Tuple[int, int, int], b: Tuple[int, int, int]) -> bool:
            return any(rotate3(a, k) == b for k in (0, 1, 2))

        def opposite_oriented(a: Tuple[int, int, int], b: Tuple[int, int, int]) -> bool:
            ar = (a[0], a[2], a[1])
            return any(rotate3(ar, k) == b for k in (0, 1, 2))

        # Build shared-face adjacency with required parity relation.
        face_map: Dict[Tuple[int, int, int], List[Tuple[int, int]]] = {}
        for ei, (_, _, _, corners) in enumerate(tet_entries):
            for fid in range(4):
                f = face_corners(corners, fid)
                key = tuple(sorted(f))
                face_map.setdefault(key, []).append((ei, fid))

        adjacency: List[List[Tuple[int, int]]] = [[] for _ in tet_entries]
        for key, refs in face_map.items():
            if len(refs) != 2:
                continue
            (e1, f1), (e2, f2) = refs
            c1 = tet_entries[e1][3]
            c2 = tet_entries[e2][3]
            base1 = face_corners(c1, f1)
            base2 = face_corners(c2, f2)

            ok = []
            for s1 in (0, 1):
                cc1 = (c1[1], c1[0], c1[2], c1[3]) if s1 else c1
                ff1 = face_corners(cc1, f1)
                for s2 in (0, 1):
                    cc2 = (c2[1], c2[0], c2[2], c2[3]) if s2 else c2
                    ff2 = face_corners(cc2, f2)
                    if opposite_oriented(ff1, ff2):
                        ok.append((s1, s2))

            if not ok:
                continue
            same_state_required = all(a == b for a, b in ok)
            req = 0 if same_state_required else 1
            adjacency[e1].append((e2, req))
            adjacency[e2].append((e1, req))

        # Solve parity constraints per connected component.
        state = [-1] * len(tet_entries)
        for start in range(len(tet_entries)):
            if state[start] != -1:
                continue
            state[start] = 0
            queue = [start]
            qh = 0
            while qh < len(queue):
                u = queue[qh]
                qh += 1
                for v, req in adjacency[u]:
                    sv = state[u] ^ req
                    if state[v] == -1:
                        state[v] = sv
                        queue.append(v)

        # Materialize oriented tet blocks with Gmsh node order: nodes -> edges.
        oriented_tet_blocks = {"tet": [], "tet2": []}
        for ei, (name, _, conn, base_corners) in enumerate(tet_entries):
            c0, c1, c2, c3 = base_corners
            if state[ei] == 1:
                corners = (c1, c0, c2, c3)
            else:
                corners = (c0, c1, c2, c3)

            if name == "tet":
                oriented_tet_blocks["tet"].append(corners)
            else:
                oriented_tet_blocks["tet2"].append(reorder_tet2_from_comsol(conn, corners))

        if "prism" in sections:
            for conn in sections["prism"]["elements"]:
                c0, c1, c2, c3, c4, c5 = conn[:6]
                if signed_prism_volume6(c0, c1, c2, c3, c4, c5) < 0.0:
                    # Flip orientation while preserving vertical edge pairing.
                    oriented_prism_block.append((c1, c0, c2, c4, c3, c5))
                else:
                    oriented_prism_block.append((c0, c1, c2, c3, c4, c5))

        if "hex" in sections:
            for conn in sections["hex"]["elements"]:
                oriented_hex_block.append(tuple(conn[:8]))

    # Build boundary-face lookup from oriented volume elements.
    # Value stores one interior probe point and one ordered face loop
    # per incident volume element.
    face_to_probe: Dict[frozenset, List[Tuple[float, float, float]]] = {}
    face_to_order: Dict[frozenset, List[Tuple[int, ...]]] = {}

    def add_face_probe(face: Tuple[int, ...], center: Tuple[float, float, float]) -> None:
        key = frozenset(face)
        face_to_probe.setdefault(key, []).append(center)
        face_to_order.setdefault(key, []).append(tuple(face))

    if max_dim == 3:
        for tet_name in ("tet", "tet2"):
            block = oriented_tet_blocks.get(tet_name, [])
            for conn in block:
                c0, c1, c2, c3 = conn[:4]
                center = (
                    (points[c0][0] + points[c1][0] + points[c2][0] + points[c3][0]) * 0.25,
                    (points[c0][1] + points[c1][1] + points[c2][1] + points[c3][1]) * 0.25,
                    (points[c0][2] + points[c1][2] + points[c2][2] + points[c3][2]) * 0.25,
                )
                faces = [
                    (c1, c2, c3),
                    (c0, c3, c2),
                    (c0, c1, c3),
                    (c0, c2, c1),
                ]
                for face in faces:
                    add_face_probe(face, center)

        for conn in oriented_prism_block:
            c0, c1, c2, c3, c4, c5 = conn[:6]
            center = (
                (points[c0][0] + points[c1][0] + points[c2][0] + points[c3][0] + points[c4][0] + points[c5][0])
                / 6.0,
                (points[c0][1] + points[c1][1] + points[c2][1] + points[c3][1] + points[c4][1] + points[c5][1])
                / 6.0,
                (points[c0][2] + points[c1][2] + points[c2][2] + points[c3][2] + points[c4][2] + points[c5][2])
                / 6.0,
            )
            for face in ((c0, c1, c2), (c3, c4, c5)):
                add_face_probe(face, center)

            for face in ((c0, c1, c4, c3), (c1, c2, c5, c4), (c2, c0, c3, c5)):
                add_face_probe(face, center)

        for conn in oriented_hex_block:
            c0, c1, c2, c3, c4, c5, c6, c7 = conn[:8]
            center = (
                (points[c0][0] + points[c1][0] + points[c2][0] + points[c3][0] + points[c4][0] + points[c5][0] + points[c6][0] + points[c7][0])
                / 8.0,
                (points[c0][1] + points[c1][1] + points[c2][1] + points[c3][1] + points[c4][1] + points[c5][1] + points[c6][1] + points[c7][1])
                / 8.0,
                (points[c0][2] + points[c1][2] + points[c2][2] + points[c3][2] + points[c4][2] + points[c5][2] + points[c6][2] + points[c7][2])
                / 8.0,
            )
            for face in (
                (c0, c1, c2, c3),
                (c4, c5, c6, c7),
                (c0, c1, c5, c4),
                (c1, c2, c6, c5),
                (c2, c3, c7, c6),
                (c3, c0, c4, c7),
            ):
                add_face_probe(face, center)

    tet_cursors = {"tet": 0, "tet2": 0}
    prism_cursor = 0
    hex_cursor = 0
    for name, (gmsh_type, expected_nodes) in element_type_map.items():
        if name not in sections:
            continue
        is_boundary_block = element_dim_map[name] < max_dim
        for local_idx, (conn, geo) in enumerate(zip(sections[name]["elements"], sections[name]["geo"])):
            tag = geo + 1 if is_boundary_block else geo
            if is_boundary_block:
                shifted_boundary_tag_count += 1
            if tag <= 0:
                tag = 1
                nonpositive_tag_count += 1
            if len(conn) != expected_nodes:
                raise ValueError(
                    f"Unexpected node count for {name}: got {len(conn)}, expected {expected_nodes}"
                )

            if is_boundary_block and max_dim == 3 and name in ("tri", "tri2"):
                key = frozenset(conn[:3])
                probe_candidates = face_to_probe.get(key, [])
                order_candidates = face_to_order.get(key, [])
                if len(probe_candidates) != 1 or len(order_candidates) != 1:
                    skipped_interior_boundary_faces += 1
                    continue
                corners = list(order_candidates[0])
                pa, pb, pc = points[corners[0]], points[corners[1]], points[corners[2]]
                pd = probe_candidates[0]
                normal = vec_cross(vec_sub(pb, pa), vec_sub(pc, pa))
                if vec_dot(normal, vec_sub(pd, pa)) > 0.0:
                    corners[1], corners[2] = corners[2], corners[1]
                if name == "tri":
                    conn = tuple(corners)
                else:
                    conn = tuple(corners) + conn[3:]

            if is_boundary_block and max_dim == 3 and name == "quad":
                key = frozenset(conn[:4])
                probe_candidates = face_to_probe.get(key, [])
                order_candidates = face_to_order.get(key, [])
                if len(probe_candidates) != 1 or len(order_candidates) != 1:
                    skipped_interior_boundary_faces += 1
                    continue
                corners = list(order_candidates[0])
                pa, pb, pc = points[corners[0]], points[corners[1]], points[corners[2]]
                pd = probe_candidates[0]
                normal = vec_cross(vec_sub(pb, pa), vec_sub(pc, pa))
                if vec_dot(normal, vec_sub(pd, pa)) > 0.0:
                    corners = [corners[0], corners[3], corners[2], corners[1]]
                conn = tuple(corners)

            if max_dim == 3 and name in ("tet", "tet2"):
                conn = oriented_tet_blocks[name][tet_cursors[name]]
                tet_cursors[name] += 1
            elif max_dim == 3 and name == "prism":
                conn = oriented_prism_block[prism_cursor]
                prism_cursor += 1
            elif max_dim == 3 and name == "hex":
                conn = oriented_hex_block[hex_cursor]
                hex_cursor += 1
            else:
                conn = reorder_quadratic_nodes(name, conn)
            elements.append((gmsh_type, tag, tag, tuple(v + 1 for v in conn)))

    if not elements:
        raise ValueError("No tri/quad/tri2/tet/hex/tet2/prism elements found")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        f.write("$Nodes\n")
        f.write(f"{len(points)}\n")
        for i, (x, y, z) in enumerate(points, start=1):
            f.write(f"{i} {x * scale:.16g} {y * scale:.16g} {z * scale:.16g}\n")
        f.write("$EndNodes\n")

        f.write("$Elements\n")
        f.write(f"{len(elements)}\n")
        for eid, (etype, ptag, etag, conn) in enumerate(elements, start=1):
            f.write(f"{eid} {etype} 2 {ptag} {etag} {' '.join(map(str, conn))}\n")
        f.write("$EndElements\n")

    return nonpositive_tag_count, shifted_boundary_tag_count, skipped_interior_boundary_faces


def main():
    parser = argparse.ArgumentParser(
        description="Convert COMSOL mphtxt to Gmsh v2 msh (single file or batch)."
    )
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument(
        "--mesh-dir",
        type=Path,
        default=Path("mesh_data"),
        help="Default directory for mesh files (default: mesh_data)",
    )
    parser.add_argument("--scale", type=float, default=1.0)
    args = parser.parse_args()

    mesh_dir = args.mesh_dir

    def resolve_in_mesh_dir(path: Path) -> Path:
        if path.is_absolute() or path.parent != Path("."):
            return path
        return mesh_dir / path

    def convert_one(input_path: Path, output_path: Path):
        points, sections = parse_mphtxt(input_path)
        nonpositive_tag_count, shifted_boundary_tag_count, skipped_interior_boundary_faces = write_msh(
            output_path, points, sections, args.scale
        )

        tri_n = len(sections.get("tri", {}).get("elements", [])) + len(
            sections.get("tri2", {}).get("elements", [])
        )
        quad_n = len(sections.get("quad", {}).get("elements", []))
        tet_n = len(sections.get("tet", {}).get("elements", [])) + len(
            sections.get("tet2", {}).get("elements", [])
        )
        hex_n = len(sections.get("hex", {}).get("elements", []))
        prism_n = len(sections.get("prism", {}).get("elements", []))
        print(
            f"Converted {input_path.name}: nodes={len(points)}, tri={tri_n}, quad={quad_n}, "
            f"tet={tet_n}, hex={hex_n}, prism={prism_n}"
        )
        if shifted_boundary_tag_count > 0:
            print(
                "  Shifted boundary geometric indices by +1 "
                f"for {shifted_boundary_tag_count} boundary elements"
            )
        if skipped_interior_boundary_faces > 0:
            print(
                "  Skipped interior tri/tri2 faces misclassified as boundary: "
                f"{skipped_interior_boundary_faces}"
            )
        if nonpositive_tag_count > 0:
            print(f"  Adjusted non-positive element tags: {nonpositive_tag_count}")
        print(f"  Output: {output_path}")

    if args.input is not None:
        input_path = resolve_in_mesh_dir(args.input)
        output_path = resolve_in_mesh_dir(args.output) if args.output is not None else input_path.with_suffix(
            ".msh"
        )
        convert_one(input_path, output_path)
        return

    if not mesh_dir.exists() or not mesh_dir.is_dir():
        raise ValueError(f"Mesh directory not found: {mesh_dir}")

    inputs = sorted(mesh_dir.glob("*.mphtxt"))
    if not inputs:
        print(f"No .mphtxt files found under: {mesh_dir}")
        return

    print(f"Found {len(inputs)} .mphtxt file(s) in {mesh_dir}")
    for input_path in inputs:
        output_path = input_path.with_suffix(".msh")
        convert_one(input_path, output_path)


if __name__ == "__main__":
    main()
