#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


@dataclass
class DataSet:
    name: str
    path: Path
    coords: List[Tuple[float, float, float]]
    values: List[float]
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


def load_dataset(path: Path, name: str, value_col: int) -> DataSet:
    coords: List[Tuple[float, float, float]] = []
    values: List[float] = []
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

            idx = value_col if value_col >= 0 else len(nums) + value_col
            if idx < 0 or idx >= len(nums):
                raise ValueError(
                    f"文件 {path} 的 value_col={value_col} 越界，当前行列数={len(nums)}"
                )
            v = nums[idx]

            coords.append((x, y, z))
            values.append(v)

    if not values:
        raise ValueError(f"文件 {path} 未解析到有效数值数据。")

    return DataSet(name=name, path=path, coords=coords, values=values, ncols=ncols)


def summary(ds: DataSet) -> str:
    n = len(ds.values)
    vmin = min(ds.values)
    vmax = max(ds.values)
    mean = sum(ds.values) / n
    return (
        f"[{ds.name}] rows={n}, cols~{ds.ncols}, min={vmin:.12g}, "
        f"max={vmax:.12g}, mean={mean:.12g}"
    )


def _qkey(p: Tuple[float, float, float], tol: float) -> Tuple[int, int, int]:
    return (int(round(p[0] / tol)), int(round(p[1] / tol)), int(round(p[2] / tol)))


def compare_by_coord(src: DataSet, ref: DataSet, tol: float) -> str:
    buckets: Dict[Tuple[int, int, int], List[int]] = {}
    for i, p in enumerate(ref.coords):
        buckets.setdefault(_qkey(p, tol), []).append(i)

    abs_err: List[float] = []
    rel_err: List[float] = []
    miss = 0

    for p, v in zip(src.coords, src.values):
        key = _qkey(p, tol)
        candidates = buckets.get(key)
        if not candidates:
            miss += 1
            continue

        best_i = candidates[0]
        best_d2 = float("inf")
        for i in candidates:
            rp = ref.coords[i]
            d2 = (p[0] - rp[0]) ** 2 + (p[1] - rp[1]) ** 2 + (p[2] - rp[2]) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_i = i

        dv = v - ref.values[best_i]
        ae = abs(dv)
        abs_err.append(ae)
        denom = abs(ref.values[best_i])
        if denom > 0.0:
            rel_err.append(ae / denom)

    matched = len(abs_err)
    total = len(src.values)
    if matched == 0:
        return (
            f"[{src.name} vs {ref.name}] 未匹配到坐标点，"
            f"建议调大 --tol（当前 {tol:g}）。"
        )

    mae = sum(abs_err) / matched
    rmse = math.sqrt(sum(e * e for e in abs_err) / matched)
    maxe = max(abs_err)
    mrel = (sum(rel_err) / len(rel_err)) if rel_err else float("nan")

    return (
        f"[{src.name} vs {ref.name}] matched={matched}/{total}, miss={miss}, "
        f"MAE={mae:.12g}, RMSE={rmse:.12g}, MaxAE={maxe:.12g}, MeanRel={mrel:.12g}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="比较一份 FEM 自输出 txt 与一份 COMSOL 导出 txt（按坐标匹配）。"
    )
    parser.add_argument("--mine", "--mine-a", required=True, type=Path, help="自己导出的 txt")
    parser.add_argument("--comsol", required=True, type=Path, help="COMSOL 导出的 txt")
    parser.add_argument(
        "--value-col",
        type=int,
        default=-1,
        help="比较值列索引（默认 -1 表示最后一列，常用于温度/位移模）",
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

        a = load_dataset(mine_path, "mine", args.value_col)
        c = load_dataset(comsol_path, "comsol", args.value_col)
    except Exception as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    print("==== 数据摘要 ====")
    print(summary(a))
    print(summary(c))

    print("\n==== 与 COMSOL 对比 ====")
    print(compare_by_coord(a, c, args.tol))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
