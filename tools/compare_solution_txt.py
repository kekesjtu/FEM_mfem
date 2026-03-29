#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field as dfield
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


@dataclass
class DataSet:
    name: str
    path: Path
    coords: List[Tuple[float, float, float]]
    values: List[List[float]]  # each row is a list of value columns
    ncols: int


def resolve_existing_path(path: Path) -> Path:
    if path.exists():
        return path

    fallback = Path("results") / path
    if fallback.exists():
        return fallback

    raise FileNotFoundError(
        f"找不到文件: {path}（也尝试了 {fallback}）"
    )


def _split_numbers(line: str) -> List[float]:
    parts = [p for p in re.split(r"[\s,]+", line.strip()) if p]
    return [float(p) for p in parts]


def load_dataset(path: Path, name: str) -> DataSet:
    coords: List[Tuple[float, float, float]] = []
    values: List[List[float]] = []
    ncols = 0

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("%") or line.startswith("#"):
                continue
            try:
                nums = _split_numbers(line)
            except ValueError:
                continue
            if len(nums) < 4:
                continue

            ncols = max(ncols, len(nums))
            x = nums[0]
            y = nums[1]
            z = nums[2] if len(nums) > 2 else 0.0

            coords.append((x, y, z))
            values.append(nums[3:])

    if not values:
        raise ValueError(f"文件 {path} 未解析到有效数值数据。")

    return DataSet(name=name, path=path, coords=coords, values=values, ncols=ncols)


def summary(ds: DataSet) -> str:
    n = len(ds.values)
    num_val_cols = len(ds.values[0]) if ds.values else 0
    all_first = [row[0] for row in ds.values]
    vmin = min(all_first)
    vmax = max(all_first)
    mean = sum(all_first) / n
    return (
        f"[{ds.name}] rows={n}, cols~{ds.ncols}, value_cols={num_val_cols}, "
        f"col0_min={vmin:.12g}, col0_max={vmax:.12g}, col0_mean={mean:.12g}"
    )


def _qkey(p: Tuple[float, float, float], tol: float) -> Tuple[int, int, int]:
    return (int(round(p[0] / tol)), int(round(p[1] / tol)), int(round(p[2] / tol)))


def compare_column(
    src_coords: List[Tuple[float, float, float]],
    src_vals: List[float],
    ref_coords: List[Tuple[float, float, float]],
    ref_vals: List[float],
    tol: float,
    label: str,
) -> str:
    buckets: Dict[Tuple[int, int, int], List[int]] = {}
    for i, p in enumerate(ref_coords):
        buckets.setdefault(_qkey(p, tol), []).append(i)

    abs_err: List[float] = []
    rel_err: List[float] = []
    miss = 0

    for p, v in zip(src_coords, src_vals):
        key = _qkey(p, tol)
        candidates = buckets.get(key)
        if not candidates:
            miss += 1
            continue

        best_i = candidates[0]
        best_d2 = float("inf")
        for i in candidates:
            rp = ref_coords[i]
            d2 = (p[0] - rp[0]) ** 2 + (p[1] - rp[1]) ** 2 + (p[2] - rp[2]) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_i = i

        dv = v - ref_vals[best_i]
        ae = abs(dv)
        abs_err.append(ae)
        denom = abs(ref_vals[best_i])
        if denom > 0.0:
            rel_err.append(ae / denom)

    matched = len(abs_err)
    total = len(src_vals)
    if matched == 0:
        return (
            f"  [{label}] 未匹配到坐标点，"
            f"建议调大 --tol（当前 {tol:g}）。"
        )

    mae = sum(abs_err) / matched
    rmse = math.sqrt(sum(e * e for e in abs_err) / matched)
    maxe = max(abs_err)
    mrel = (sum(rel_err) / len(rel_err)) if rel_err else float("nan")

    return (
        f"  [{label}] matched={matched}/{total}, miss={miss}, "
        f"MAE={mae:.12g}, RMSE={rmse:.12g}, MaxAE={maxe:.12g}, MeanRel={mrel:.12g}"
    )


def detect_mode(mine: DataSet, comsol: DataSet) -> str:
    """Auto-detect scalar vs vector based on value column count."""
    mine_vcols = len(mine.values[0]) if mine.values else 0
    comsol_vcols = len(comsol.values[0]) if comsol.values else 0
    # Vector field: our output has 4+ value cols (ux, uy, uz, magnitude)
    # COMSOL vector export typically has 3 cols (ux, uy, uz) or more
    if mine_vcols >= 4 and comsol_vcols >= 3:
        return "vector"
    return "scalar"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="比较一份 FEM 自输出 txt 与一份 COMSOL 导出 txt（按坐标匹配）。\n"
                    "支持稳态（标量/矢量自动检测）和瞬态模式。"
    )
    parser.add_argument("--mine", "--mine-a", required=True, type=Path, help="自己导出的 txt")
    parser.add_argument("--comsol", required=True, type=Path, help="COMSOL 导出的 txt")
    parser.add_argument(
        "--value-col",
        type=int,
        default=None,
        help="比较值列索引（默认自动检测；指定后强制标量模式，如 -1 表示最后一列）",
    )
    parser.add_argument(
        "--transient",
        action="store_true",
        help="瞬态模式：按 (fields_per_step × n_steps) 解析值列",
    )
    parser.add_argument(
        "--fields-per-step",
        type=int,
        default=3,
        help="瞬态模式下每时间步的场数（默认 3: V, T, disp）",
    )
    parser.add_argument(
        "--field-names",
        type=str,
        default="V,T,disp",
        help="瞬态模式下各场名称，逗号分隔（默认 V,T,disp）",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-8,
        help="坐标匹配容差（量化键基准，默认 1e-8）",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        mine_path = resolve_existing_path(args.mine)
        comsol_path = resolve_existing_path(args.comsol)

        mine = load_dataset(mine_path, "mine")
        comsol = load_dataset(comsol_path, "comsol")
    except Exception as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    print("==== 数据摘要 ====")
    print(summary(mine))
    print(summary(comsol))

    # --- Transient mode ---
    if args.transient:
        fps = args.fields_per_step
        field_names = [s.strip() for s in args.field_names.split(",")]
        if len(field_names) != fps:
            print(f"错误: --field-names 数量({len(field_names)})与 --fields-per-step({fps})不一致",
                  file=sys.stderr)
            return 1

        mine_vcols = len(mine.values[0]) if mine.values else 0
        comsol_vcols = len(comsol.values[0]) if comsol.values else 0
        n_steps_mine = mine_vcols // fps
        n_steps_comsol = comsol_vcols // fps
        n_steps = min(n_steps_mine, n_steps_comsol)
        if n_steps == 0:
            print("错误: 无法确定时间步数", file=sys.stderr)
            return 1

        print(f"\n==== 瞬态对比 (时间步数={n_steps}, 场数={fps}) ====")
        # Build coordinate index once
        buckets: Dict[Tuple[int, int, int], List[int]] = {}
        for i, p in enumerate(comsol.coords):
            buckets.setdefault(_qkey(p, args.tol), []).append(i)

        # Match coordinates
        matches: List[Tuple[int, int]] = []  # (mine_idx, comsol_idx)
        for mi, p in enumerate(mine.coords):
            key = _qkey(p, args.tol)
            candidates = buckets.get(key)
            if not candidates:
                continue
            best_i = candidates[0]
            best_d2 = float("inf")
            for ci in candidates:
                rp = comsol.coords[ci]
                d2 = (p[0]-rp[0])**2 + (p[1]-rp[1])**2 + (p[2]-rp[2])**2
                if d2 < best_d2:
                    best_d2 = d2
                    best_i = ci
            matches.append((mi, best_i))

        n_matched = len(matches)
        n_total = len(mine.coords)
        print(f"  坐标匹配: {n_matched}/{n_total}")

        if n_matched == 0:
            print("  未匹配到坐标点，建议调大 --tol。")
            return 1

        # Per-field aggregate stats
        for fi, fname in enumerate(field_names):
            step_maes: List[float] = []
            step_maxes: List[float] = []
            for si in range(n_steps):
                col = si * fps + fi
                abs_errs = []
                for mi, ci in matches:
                    mv = mine.values[mi][col] if col < len(mine.values[mi]) else 0.0
                    cv = comsol.values[ci][col] if col < len(comsol.values[ci]) else 0.0
                    abs_errs.append(abs(mv - cv))
                mae = sum(abs_errs) / len(abs_errs)
                maxe = max(abs_errs)
                step_maes.append(mae)
                step_maxes.append(maxe)

            overall_mae = sum(step_maes) / len(step_maes)
            overall_max = max(step_maxes)
            worst_step = step_maes.index(max(step_maes))

            print(f"\n  [{fname}] 全时间步平均MAE={overall_mae:.6g}, 最大MaxAE={overall_max:.6g}, 最差步={worst_step}")
            # Print per-step table
            print(f"    {'step':>5s}  {'MAE':>14s}  {'MaxAE':>14s}")
            for si in range(n_steps):
                print(f"    {si:5d}  {step_maes[si]:14.6g}  {step_maxes[si]:14.6g}")

        return 0

    # --- Steady-state modes below ---

    # Determine mode
    if args.value_col is not None:
        # Legacy single-column mode
        mine_vals = []
        for row in mine.values:
            idx = args.value_col if args.value_col >= 0 else len(row) + args.value_col
            mine_vals.append(row[idx])
        comsol_vals = []
        for row in comsol.values:
            idx = args.value_col if args.value_col >= 0 else len(row) + args.value_col
            comsol_vals.append(row[idx])

        print(f"\n==== 与 COMSOL 对比 (value_col={args.value_col}) ====")
        print(compare_column(mine.coords, mine_vals, comsol.coords, comsol_vals,
                             args.tol, f"col{args.value_col}"))
    else:
        mode = detect_mode(mine, comsol)
        if mode == "vector":
            print(f"\n==== 与 COMSOL 对比 (矢量场自动检测) ====")
            # Component names
            comp_names = ["ux", "uy", "uz"]
            num_comps = min(3, len(mine.values[0]), len(comsol.values[0]))
            for c in range(num_comps):
                mine_c = [row[c] for row in mine.values]
                comsol_c = [row[c] for row in comsol.values]
                print(compare_column(mine.coords, mine_c, comsol.coords, comsol_c,
                                     args.tol, comp_names[c]))

            # Magnitude comparison
            def magnitude(row: List[float], n: int) -> float:
                return math.sqrt(sum(row[i] ** 2 for i in range(min(n, len(row)))))

            mine_mag = [magnitude(row, num_comps) for row in mine.values]
            comsol_mag = [magnitude(row, num_comps) for row in comsol.values]
            print(compare_column(mine.coords, mine_mag, comsol.coords, comsol_mag,
                                 args.tol, "magnitude"))
        else:
            # Scalar: compare first value column (index 0 in values, which is col 3 in file)
            mine_vals = [row[0] for row in mine.values]
            comsol_vals = [row[0] for row in comsol.values]
            print(f"\n==== 与 COMSOL 对比 (标量场) ====")
            print(compare_column(mine.coords, mine_vals, comsol.coords, comsol_vals,
                                 args.tol, "value"))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
