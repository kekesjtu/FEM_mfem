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

    for name, gmsh_type in (("tri", 2), ("tet", 4)):
        if name not in sections:
            continue
        for conn, geo in zip(sections[name]["elements"], sections[name]["geo"]):
            tag = geo
            if tag <= 0:
                tag = 1
                nonpositive_tag_count += 1
            elements.append((gmsh_type, tag, tag, tuple(v + 1 for v in conn)))

    if not elements:
        raise ValueError("No tri/tet elements found")

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

    return nonpositive_tag_count


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
        nonpositive_tag_count = write_msh(output_path, points, sections, args.scale)

        tri_n = len(sections.get("tri", {}).get("elements", []))
        tet_n = len(sections.get("tet", {}).get("elements", []))
        print(f"Converted {input_path.name}: nodes={len(points)}, tri={tri_n}, tet={tet_n}")
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
