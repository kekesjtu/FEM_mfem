#!/usr/bin/env python3
from __future__ import annotations

import argparse
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

    def signed_tet_volume6(c0: int, c1: int, c2: int, c3: int) -> float:
        p0, p1, p2, p3 = points[c0], points[c1], points[c2], points[c3]
        return vec_dot(vec_cross(vec_sub(p1, p0), vec_sub(p2, p0)), vec_sub(p3, p0))

    def signed_prism_volume6(c0: int, c1: int, c2: int, c3: int, c4: int, c5: int) -> float:
        return (
            signed_tet_volume6(c0, c1, c2, c3)
            + signed_tet_volume6(c1, c2, c4, c3)
            + signed_tet_volume6(c2, c4, c5, c3)
        )

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

        gmsh_edges = ((0, 1), (1, 2), (2, 0), (0, 3), (2, 3), (1, 3))
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

    # For 3D meshes: orient volume elements in-place and build boundary face lookup.
    # Positive signed volume guarantees globally consistent orientation for
    # conforming meshes, so no BFS/adjacency graph is needed.
    face_to_probe: Dict[tuple, Tuple[float, float, float]] = {}
    face_to_order: Dict[tuple, Tuple[int, ...]] = {}
    face_shared: set = set()

    if max_dim == 3:
        # Orient tets in-place to have positive signed volume.
        for tet_name in ("tet", "tet2"):
            if tet_name not in sections:
                continue
            elems = sections[tet_name]["elements"]
            for i, conn in enumerate(elems):
                c0, c1, c2, c3 = conn[:4]
                if signed_tet_volume6(c0, c1, c2, c3) < 0.0:
                    corners = (c1, c0, c2, c3)
                else:
                    corners = (c0, c1, c2, c3)
                if tet_name == "tet":
                    elems[i] = corners
                else:
                    elems[i] = reorder_tet2_from_comsol(conn, corners)

        # Orient prisms in-place.
        if "prism" in sections:
            elems = sections["prism"]["elements"]
            for i, conn in enumerate(elems):
                c0, c1, c2, c3, c4, c5 = conn[:6]
                if signed_prism_volume6(c0, c1, c2, c3, c4, c5) < 0.0:
                    elems[i] = (c1, c0, c2, c4, c3, c5)
                else:
                    elems[i] = (c0, c1, c2, c3, c4, c5)

        # Normalize hex elements.
        if "hex" in sections:
            elems = sections["hex"]["elements"]
            for i, conn in enumerate(elems):
                elems[i] = tuple(conn[:8])

        # Collect boundary face keys and boundary vertices for fast pre-filtering.
        boundary_face_keys: set = set()
        boundary_vertices: set = set()
        for bname in ("tri", "tri2"):
            if bname in sections:
                for conn in sections[bname]["elements"]:
                    boundary_face_keys.add(tuple(sorted(conn[:3])))
                    boundary_vertices.add(conn[0])
                    boundary_vertices.add(conn[1])
                    boundary_vertices.add(conn[2])
        if "quad" in sections:
            for conn in sections["quad"]["elements"]:
                boundary_face_keys.add(tuple(sorted(conn[:4])))
                boundary_vertices.update(conn[:4])

        # Build face_to_probe/face_to_order ONLY for faces that match boundary elements.
        # This avoids storing data for millions of interior faces.
        bv = boundary_vertices
        bfk = boundary_face_keys

        def _register_face(face_verts: tuple, center: Tuple[float, float, float]) -> None:
            key = tuple(sorted(face_verts))
            if key not in bfk:
                return
            if key in face_to_probe:
                face_shared.add(key)
            else:
                face_to_probe[key] = center
                face_to_order[key] = face_verts

        for tet_name in ("tet", "tet2"):
            if tet_name not in sections:
                continue
            for conn in sections[tet_name]["elements"]:
                c0, c1, c2, c3 = conn[:4]
                if (c0 in bv) + (c1 in bv) + (c2 in bv) + (c3 in bv) < 3:
                    continue
                p0, p1, p2, p3 = points[c0], points[c1], points[c2], points[c3]
                center = (
                    (p0[0] + p1[0] + p2[0] + p3[0]) * 0.25,
                    (p0[1] + p1[1] + p2[1] + p3[1]) * 0.25,
                    (p0[2] + p1[2] + p2[2] + p3[2]) * 0.25,
                )
                # MFEM face vertex ordering for positive-volume tet
                _register_face((c1, c2, c3), center)
                _register_face((c0, c3, c2), center)
                _register_face((c0, c1, c3), center)
                _register_face((c0, c2, c1), center)

        if "prism" in sections:
            for conn in sections["prism"]["elements"]:
                c0, c1, c2, c3, c4, c5 = conn[:6]
                if sum(1 for v in (c0, c1, c2, c3, c4, c5) if v in bv) < 3:
                    continue
                p = [points[v] for v in (c0, c1, c2, c3, c4, c5)]
                center = (
                    sum(pt[0] for pt in p) / 6.0,
                    sum(pt[1] for pt in p) / 6.0,
                    sum(pt[2] for pt in p) / 6.0,
                )
                for face in ((c0, c1, c2), (c3, c4, c5)):
                    _register_face(face, center)
                for face in ((c0, c1, c4, c3), (c1, c2, c5, c4), (c2, c0, c3, c5)):
                    _register_face(face, center)

        if "hex" in sections:
            for conn in sections["hex"]["elements"]:
                c0, c1, c2, c3, c4, c5, c6, c7 = conn[:8]
                if sum(1 for v in conn[:8] if v in bv) < 4:
                    continue
                p = [points[v] for v in conn[:8]]
                center = (
                    sum(pt[0] for pt in p) / 8.0,
                    sum(pt[1] for pt in p) / 8.0,
                    sum(pt[2] for pt in p) / 8.0,
                )
                for face in (
                    (c0, c1, c2, c3), (c4, c5, c6, c7),
                    (c0, c1, c5, c4), (c1, c2, c6, c5),
                    (c2, c3, c7, c6), (c3, c0, c4, c7),
                ):
                    _register_face(face, center)

    for name, (gmsh_type, expected_nodes) in element_type_map.items():
        if name not in sections:
            continue
        is_boundary_block = element_dim_map[name] < max_dim
        for conn, geo in zip(sections[name]["elements"], sections[name]["geo"]):
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
                key = tuple(sorted(conn[:3]))
                if key in face_shared or key not in face_to_probe:
                    skipped_interior_boundary_faces += 1
                    continue
                corners = list(face_to_order[key])
                pa, pb, pc = points[corners[0]], points[corners[1]], points[corners[2]]
                pd = face_to_probe[key]
                normal = vec_cross(vec_sub(pb, pa), vec_sub(pc, pa))
                if vec_dot(normal, vec_sub(pd, pa)) > 0.0:
                    corners[1], corners[2] = corners[2], corners[1]
                if name == "tri":
                    conn = tuple(corners)
                else:
                    conn = reorder_tri2_from_comsol(conn, tuple(corners))

            elif is_boundary_block and max_dim == 3 and name == "quad":
                key = tuple(sorted(conn[:4]))
                if key in face_shared or key not in face_to_probe:
                    skipped_interior_boundary_faces += 1
                    continue
                corners = list(face_to_order[key])
                pa, pb, pc = points[corners[0]], points[corners[1]], points[corners[2]]
                pd = face_to_probe[key]
                normal = vec_cross(vec_sub(pb, pa), vec_sub(pc, pa))
                if vec_dot(normal, vec_sub(pd, pa)) > 0.0:
                    corners = [corners[0], corners[3], corners[2], corners[1]]
                conn = tuple(corners)

            elif max_dim == 3 and not is_boundary_block:
                pass  # Already oriented in-place
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
