#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import math
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Coord = Tuple[float, float, float]


@dataclass
class DataSet:
    name: str
    path: Path
    coords: List[Coord]
    values: List[List[float]]
    ncols: int
    header_times: List[float]


@dataclass
class PointMatch:
    target: Coord
    mine_idx: int
    comsol_idx: int
    mine_dist: float
    comsol_dist: float


@dataclass
class Metric:
    count: int
    mae: float
    rmse: float
    maxae: float
    mean_rel: float
    worst_step: int
    worst_time: float


def resolve_existing_path(path: Path) -> Path:
    if path.exists():
        return path

    fallback = Path("results") / path
    if fallback.exists():
        return fallback

    raise FileNotFoundError(f"找不到文件: {path}（也尝试了 {fallback}）")


def split_numbers(line: str) -> List[float]:
    parts = [p for p in re.split(r"[\s,]+", line.strip()) if p]
    return [float(p) for p in parts]


def parse_header_times(line: str) -> List[float]:
    # Handles both "V @ t=0.001" and "V (V) @ t=0.001".
    pattern = re.compile(
        r"@\s*t\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    )
    return [float(m.group(1)) for m in pattern.finditer(line)]


def load_dataset(path: Path, name: str, coord_dim: int) -> DataSet:
    coords: List[Coord] = []
    values: List[List[float]] = []
    ncols = 0
    header_times: List[float] = []

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("%") or line.startswith("#"):
                header_times.extend(parse_header_times(line))
                continue

            try:
                nums = split_numbers(line)
            except ValueError:
                continue

            if coord_dim == 2:
                if len(nums) < 3:
                    continue
                x, y = nums[0], nums[1]
                z = 0.0
                row_values = nums[2:]
            else:
                if len(nums) < 4:
                    continue
                x, y, z = nums[0], nums[1], nums[2]
                row_values = nums[3:]

            if not row_values:
                continue

            ncols = max(ncols, len(nums))
            coords.append((x, y, z))
            values.append(row_values)

    if not values:
        raise ValueError(f"文件 {path} 未解析到有效数值数据。")

    return DataSet(
        name=name, path=path, coords=coords, values=values, ncols=ncols, header_times=header_times
    )


def qkey(p: Coord, tol: float) -> Tuple[int, int, int]:
    return (int(round(p[0] / tol)), int(round(p[1] / tol)), int(round(p[2] / tol)))


def neighbor_keys(key: Tuple[int, int, int]) -> Iterable[Tuple[int, int, int]]:
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield (key[0] + dx, key[1] + dy, key[2] + dz)


def build_buckets(coords: Sequence[Coord], tol: float) -> Dict[Tuple[int, int, int], List[int]]:
    buckets: Dict[Tuple[int, int, int], List[int]] = {}
    for i, p in enumerate(coords):
        buckets.setdefault(qkey(p, tol), []).append(i)
    return buckets


def dist2(a: Coord, b: Coord) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def find_nearest(
    point: Coord,
    coords: Sequence[Coord],
    buckets: Dict[Tuple[int, int, int], List[int]],
    tol: float,
) -> Optional[Tuple[int, float]]:
    key = qkey(point, tol)
    best_idx: Optional[int] = None
    best_d2 = float("inf")
    tol2 = tol * tol

    for nkey in neighbor_keys(key):
        for idx in buckets.get(nkey, []):
            d2 = dist2(point, coords[idx])
            if d2 <= tol2 and d2 < best_d2:
                best_idx = idx
                best_d2 = d2

    if best_idx is None:
        return None
    return best_idx, math.sqrt(best_d2)


def parse_point_token(token: str) -> Coord:
    nums = split_numbers(token)
    if len(nums) == 2:
        return (nums[0], nums[1], 0.0)
    if len(nums) >= 3:
        return (nums[0], nums[1], nums[2])
    raise ValueError(f"采样点坐标至少需要 2 个数: {token!r}")


def parse_points(points: str) -> List[Coord]:
    if not points.strip():
        return []
    tokens = [t.strip() for t in points.split(";") if t.strip()]
    return [parse_point_token(token) for token in tokens]


def load_point_file(path: Path) -> List[Coord]:
    out: List[Coord] = []
    with path.open("r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("%"):
                continue
            try:
                out.append(parse_point_token(line))
            except Exception as exc:
                raise ValueError(f"{path}:{lineno}: {exc}") from exc
    return out


def match_all_coords(mine: DataSet, comsol: DataSet, tol: float) -> List[PointMatch]:
    comsol_buckets = build_buckets(comsol.coords, tol)
    matches: List[PointMatch] = []
    for mi, p in enumerate(mine.coords):
        found = find_nearest(p, comsol.coords, comsol_buckets, tol)
        if found is None:
            continue
        ci, d = found
        matches.append(PointMatch(target=p, mine_idx=mi, comsol_idx=ci, mine_dist=0.0, comsol_dist=d))
    return matches


def match_requested_points(
    points: Sequence[Coord], mine: DataSet, comsol: DataSet, tol: float
) -> Tuple[List[PointMatch], List[Coord]]:
    mine_buckets = build_buckets(mine.coords, tol)
    comsol_buckets = build_buckets(comsol.coords, tol)
    matches: List[PointMatch] = []
    missed: List[Coord] = []

    for p in points:
        mine_found = find_nearest(p, mine.coords, mine_buckets, tol)
        comsol_found = find_nearest(p, comsol.coords, comsol_buckets, tol)
        if mine_found is None or comsol_found is None:
            missed.append(p)
            continue
        mi, md = mine_found
        ci, cd = comsol_found
        matches.append(PointMatch(target=p, mine_idx=mi, comsol_idx=ci, mine_dist=md, comsol_dist=cd))

    return matches, missed


def value_column(step: int, field_idx: int, fields_per_step: int) -> int:
    return step * fields_per_step + field_idx


def infer_times(ds: DataSet, fields_per_step: int, n_steps: int) -> List[float]:
    if len(ds.header_times) >= fields_per_step * n_steps:
        return ds.header_times[: fields_per_step * n_steps : fields_per_step]
    if len(ds.header_times) >= n_steps:
        # Fallback for headers that list one time per step.
        return ds.header_times[:n_steps]
    return [float(i) for i in range(n_steps)]


def parse_field_names(raw: str, fields_per_step: int) -> List[str]:
    names = [s.strip() for s in raw.split(",") if s.strip()]
    if len(names) != fields_per_step:
        raise ValueError(
            f"--field-names 数量({len(names)})与 --fields-per-step({fields_per_step})不一致"
        )
    return names


def parse_field_selection(raw: str, field_names: Sequence[str]) -> List[int]:
    if not raw.strip():
        return list(range(len(field_names)))

    name_to_idx = {name.lower(): i for i, name in enumerate(field_names)}
    selected: List[int] = []
    for token in [s.strip() for s in raw.split(",") if s.strip()]:
        if re.fullmatch(r"\d+", token):
            idx = int(token)
        else:
            key = token.lower()
            if key not in name_to_idx:
                raise ValueError(f"未知字段 {token!r}，可选字段: {', '.join(field_names)}")
            idx = name_to_idx[key]

        if idx < 0 or idx >= len(field_names):
            raise ValueError(f"字段索引越界: {idx}")
        if idx not in selected:
            selected.append(idx)

    return selected


def compute_metric(errors: Sequence[float], rel_errors: Sequence[float], times: Sequence[float]) -> Metric:
    if not errors:
        return Metric(0, float("nan"), float("nan"), float("nan"), float("nan"), -1, float("nan"))
    maxae = max(errors)
    worst_step = errors.index(maxae)
    return Metric(
        count=len(errors),
        mae=sum(errors) / len(errors),
        rmse=math.sqrt(sum(e * e for e in errors) / len(errors)),
        maxae=maxae,
        mean_rel=(sum(rel_errors) / len(rel_errors)) if rel_errors else float("nan"),
        worst_step=worst_step,
        worst_time=times[worst_step] if 0 <= worst_step < len(times) else float("nan"),
    )


def fmt_coord(p: Coord) -> str:
    return f"({p[0]:.12g}, {p[1]:.12g}, {p[2]:.12g})"


def print_metric(prefix: str, metric: Metric) -> None:
    print(
        f"{prefix} count={metric.count}, MAE={metric.mae:.12g}, "
        f"RMSE={metric.rmse:.12g}, MaxAE={metric.maxae:.12g}, "
        f"MeanRel={metric.mean_rel:.12g}, worst_step={metric.worst_step}, "
        f"worst_time={metric.worst_time:.12g}"
    )


def sanitize_filename(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    return cleaned.strip("_") or "field"


def padded_range(values: Sequence[float]) -> Tuple[float, float]:
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return 0.0, 1.0
    vmin = min(finite)
    vmax = max(finite)
    if vmin == vmax:
        pad = abs(vmin) * 0.05 if vmin else 1.0
    else:
        pad = (vmax - vmin) * 0.08
    return vmin - pad, vmax + pad


def tick_values(vmin: float, vmax: float, count: int = 6) -> List[float]:
    if count <= 1 or vmin == vmax:
        return [vmin]
    step = (vmax - vmin) / (count - 1)
    return [vmin + i * step for i in range(count)]


def save_svg_plot(
    path: Path,
    field_name: str,
    point_index: int,
    target: Coord,
    mine_times: Sequence[float],
    comsol_times: Sequence[float],
    mine_series: Sequence[float],
    comsol_series: Sequence[float],
    metric: Metric,
) -> None:
    width = 980
    height = 620
    left = 92
    right = 38
    top = 76
    bottom = 86
    plot_w = width - left - right
    plot_h = height - top - bottom

    x_min, x_max = padded_range([*mine_times, *comsol_times])
    y_min, y_max = padded_range([*mine_series, *comsol_series])

    def sx(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    def polyline(times: Sequence[float], values: Sequence[float]) -> str:
        pts = []
        for x, y in zip(times, values):
            if math.isfinite(x) and math.isfinite(y):
                pts.append(f"{sx(x):.3f},{sy(y):.3f}")
        return " ".join(pts)

    def text(x: float, y: float, body: str, size: int = 14, anchor: str = "start") -> str:
        return (
            f'<text x="{x:.3f}" y="{y:.3f}" font-size="{size}" '
            f'font-family="Arial, sans-serif" text-anchor="{anchor}" fill="#1f2933">'
            f"{html.escape(body)}</text>"
        )

    x_ticks = tick_values(x_min, x_max)
    y_ticks = tick_values(y_min, y_max)

    lines: List[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        text(width / 2, 34, f"Point {point_index} - {field_name}", 22, "middle"),
        text(
            width / 2,
            58,
            f"target=({target[0]:.6g}, {target[1]:.6g}, {target[2]:.6g})",
            13,
            "middle",
        ),
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" '
        'fill="#fbfcfe" stroke="#d0d7de"/>',
    ]

    for xt in x_ticks:
        x = sx(xt)
        lines.append(f'<line x1="{x:.3f}" y1="{top}" x2="{x:.3f}" y2="{top + plot_h}" '
                     'stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(text(x, top + plot_h + 24, f"{xt:.6g}", 12, "middle"))

    for yt in y_ticks:
        y = sy(yt)
        lines.append(f'<line x1="{left}" y1="{y:.3f}" x2="{left + plot_w}" y2="{y:.3f}" '
                     'stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(text(left - 12, y + 4, f"{yt:.6g}", 12, "end"))

    lines.extend(
        [
            f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" '
            f'y2="{top + plot_h}" stroke="#111827" stroke-width="1.4"/>',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" '
            'stroke="#111827" stroke-width="1.4"/>',
            text(left + plot_w / 2, height - 26, "time", 14, "middle"),
            (
                f'<text x="24" y="{top + plot_h / 2:.3f}" font-size="14" '
                'font-family="Arial, sans-serif" text-anchor="middle" fill="#1f2933" '
                f'transform="rotate(-90 24 {top + plot_h / 2:.3f})">'
                f"{html.escape(field_name)}</text>"
            ),
            f'<polyline points="{polyline(mine_times, mine_series)}" fill="none" '
            'stroke="#1f77b4" stroke-width="2.4"/>',
            f'<polyline points="{polyline(comsol_times, comsol_series)}" fill="none" '
            'stroke="#d62728" stroke-width="2.4" stroke-dasharray="7 4"/>',
        ]
    )

    for x, y in zip(mine_times, mine_series):
        if math.isfinite(x) and math.isfinite(y):
            lines.append(f'<circle cx="{sx(x):.3f}" cy="{sy(y):.3f}" r="3" fill="#1f77b4"/>')
    for x, y in zip(comsol_times, comsol_series):
        if math.isfinite(x) and math.isfinite(y):
            lines.append(f'<circle cx="{sx(x):.3f}" cy="{sy(y):.3f}" r="3" fill="#d62728"/>')

    legend_x = left + plot_w - 178
    legend_y = top + 18
    lines.extend(
        [
            f'<rect x="{legend_x - 12}" y="{legend_y - 18}" width="178" height="62" '
            'fill="#ffffff" stroke="#d0d7de"/>',
            f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 34}" y2="{legend_y}" '
            'stroke="#1f77b4" stroke-width="2.4"/>',
            text(legend_x + 44, legend_y + 5, "FEM", 13),
            f'<line x1="{legend_x}" y1="{legend_y + 26}" x2="{legend_x + 34}" '
            f'y2="{legend_y + 26}" stroke="#d62728" stroke-width="2.4" stroke-dasharray="7 4"/>',
            text(legend_x + 44, legend_y + 31, "COMSOL", 13),
            text(
                left,
                height - 56,
                f"MAE={metric.mae:.6g}, RMSE={metric.rmse:.6g}, "
                f"MaxAE={metric.maxae:.6g}, MeanRel={metric.mean_rel:.6g}",
                13,
            ),
            "</svg>",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")


def save_png_plot(
    path: Path,
    field_name: str,
    point_index: int,
    target: Coord,
    mine_times: Sequence[float],
    comsol_times: Sequence[float],
    mine_series: Sequence[float],
    comsol_series: Sequence[float],
    metric: Metric,
    dpi: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError("输出 PNG 需要安装 matplotlib；当前环境建议使用默认 SVG。") from exc

    fig, ax = plt.subplots(figsize=(9.8, 6.2), dpi=dpi)
    ax.plot(mine_times, mine_series, marker="o", linewidth=2.2, label="FEM")
    ax.plot(comsol_times, comsol_series, marker="s", linewidth=2.2, linestyle="--", label="COMSOL")
    ax.set_title(
        f"Point {point_index} - {field_name}\n"
        f"target=({target[0]:.6g}, {target[1]:.6g}, {target[2]:.6g})"
    )
    ax.set_xlabel("time")
    ax.set_ylabel(field_name)
    ax.grid(True, alpha=0.28)
    ax.legend()
    fig.text(
        0.12,
        0.02,
        f"MAE={metric.mae:.6g}, RMSE={metric.rmse:.6g}, "
        f"MaxAE={metric.maxae:.6g}, MeanRel={metric.mean_rel:.6g}",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path)
    plt.close(fig)


def save_curve_plot(
    path: Path,
    plot_format: str,
    field_name: str,
    point_index: int,
    target: Coord,
    mine_times: Sequence[float],
    comsol_times: Sequence[float],
    mine_series: Sequence[float],
    comsol_series: Sequence[float],
    metric: Metric,
    dpi: int,
) -> None:
    if plot_format == "svg":
        save_svg_plot(
            path,
            field_name,
            point_index,
            target,
            mine_times,
            comsol_times,
            mine_series,
            comsol_series,
            metric,
        )
        return
    if plot_format == "png":
        save_png_plot(
            path,
            field_name,
            point_index,
            target,
            mine_times,
            comsol_times,
            mine_series,
            comsol_series,
            metric,
            dpi,
        )
        return
    raise ValueError(f"Unsupported plot format: {plot_format}")


def compare_curves(
    mine: DataSet,
    comsol: DataSet,
    matches: Sequence[PointMatch],
    fields_per_step: int,
    field_names: Sequence[str],
    selected_fields: Sequence[int],
    n_steps: int,
    mine_times: Sequence[float],
    comsol_times: Sequence[float],
    rel_floor: float,
    output_csv: Optional[Path],
    plot_dir: Optional[Path],
    plot_format: str,
    plot_dpi: int,
    print_points: bool,
) -> int:
    aggregate_errors: Dict[int, List[float]] = {fi: [] for fi in selected_fields}
    aggregate_rel_errors: Dict[int, List[float]] = {fi: [] for fi in selected_fields}
    all_errors: List[float] = []
    all_rel_errors: List[float] = []

    csv_rows: List[Dict[str, object]] = []
    plot_paths: List[Path] = []
    if plot_dir is not None:
        plot_dir.mkdir(parents=True, exist_ok=True)

    if print_points:
        print("\n==== 采样点曲线误差 ====")

    for pi, match in enumerate(matches):
        mine_coord = mine.coords[match.mine_idx]
        comsol_coord = comsol.coords[match.comsol_idx]
        if print_points:
            print(
                f"\n[point {pi}] target={fmt_coord(match.target)}, "
                f"mine={fmt_coord(mine_coord)}, comsol={fmt_coord(comsol_coord)}, "
                f"dist=({match.mine_dist:.4g}, {match.comsol_dist:.4g})"
            )

        for fi in selected_fields:
            errors: List[float] = []
            rel_errors: List[float] = []
            mine_series: List[float] = []
            comsol_series: List[float] = []

            for si in range(n_steps):
                col = value_column(si, fi, fields_per_step)
                mv = mine.values[match.mine_idx][col]
                cv = comsol.values[match.comsol_idx][col]
                ae = abs(mv - cv)
                denom = max(abs(cv), rel_floor)
                re = ae / denom

                errors.append(ae)
                rel_errors.append(re)
                mine_series.append(mv)
                comsol_series.append(cv)
                aggregate_errors[fi].append(ae)
                aggregate_rel_errors[fi].append(re)
                all_errors.append(ae)
                all_rel_errors.append(re)

                if output_csv is not None:
                    csv_rows.append(
                        {
                            "point": pi,
                            "target_x": match.target[0],
                            "target_y": match.target[1],
                            "target_z": match.target[2],
                            "mine_x": mine_coord[0],
                            "mine_y": mine_coord[1],
                            "mine_z": mine_coord[2],
                            "comsol_x": comsol_coord[0],
                            "comsol_y": comsol_coord[1],
                            "comsol_z": comsol_coord[2],
                            "field": field_names[fi],
                            "field_index": fi,
                            "step": si,
                            "time_mine": mine_times[si],
                            "time_comsol": comsol_times[si],
                            "mine": mv,
                            "comsol": cv,
                            "abs_error": ae,
                            "rel_error": re,
                        }
                    )

            metric = compute_metric(errors, rel_errors, mine_times)
            if print_points:
                print_metric(f"  [{field_names[fi]}]", metric)
            if plot_dir is not None:
                filename = f"point_{pi:03d}_{sanitize_filename(field_names[fi])}.{plot_format}"
                plot_path = plot_dir / filename
                save_curve_plot(
                    plot_path,
                    plot_format,
                    field_names[fi],
                    pi,
                    match.target,
                    mine_times,
                    comsol_times,
                    mine_series,
                    comsol_series,
                    metric,
                    plot_dpi,
                )
                plot_paths.append(plot_path)

    print("\n==== 字段聚合误差 ====")
    for fi in selected_fields:
        metric = compute_metric(aggregate_errors[fi], aggregate_rel_errors[fi], mine_times * len(matches))
        print_metric(f"  [{field_names[fi]}]", metric)

    print("\n==== 总体误差 ====")
    total_metric = compute_metric(all_errors, all_rel_errors, mine_times * len(matches) * len(selected_fields))
    print_metric("  [all]", total_metric)

    if output_csv is not None:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()) if csv_rows else [])
            if csv_rows:
                writer.writeheader()
                writer.writerows(csv_rows)
        print(f"\n已写出 CSV 明细: {output_csv}")

    if plot_paths:
        print(f"\n已写出曲线图片: {plot_dir} ({len(plot_paths)} files)")
        for path in plot_paths[:10]:
            print(f"  {path}")
        if len(plot_paths) > 10:
            print(f"  ... 其余 {len(plot_paths) - 10} 个省略")

    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "比较本项目 transient.txt 与 COMSOL 瞬态宽表在采样点上的时间曲线。\n"
            "默认格式为 x y z 后接 V,T,disp 按时间步交错排列。"
        )
    )
    parser.add_argument("--mine", required=True, type=Path, help="本项目导出的 transient.txt")
    parser.add_argument("--comsol", required=True, type=Path, help="COMSOL 导出的参考 txt")
    parser.add_argument(
        "--coord-dim",
        type=int,
        choices=(2, 3),
        default=3,
        help="坐标列数，默认 3（x y z）",
    )
    parser.add_argument(
        "--fields-per-step",
        type=int,
        default=3,
        help="每个时间步的字段数，默认 3",
    )
    parser.add_argument(
        "--field-names",
        type=str,
        default="V,T,disp",
        help="字段名，逗号分隔，默认 V,T,disp",
    )
    parser.add_argument(
        "--fields",
        type=str,
        default="",
        help="只比较部分字段，可填字段名或 0-based 索引，逗号分隔；默认全部字段",
    )
    parser.add_argument(
        "--random-points",
        "--sample-size",
        type=int,
        default=None,
        help="从可匹配点中随机选取 N 个点；未指定采样点时默认选 5 个",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="随机采样种子，默认 12345",
    )
    parser.add_argument(
        "--points",
        type=str,
        default="",
        help="显式采样点，如 'x,y,z; x,y,z'；2D 点可写 'x,y'",
    )
    parser.add_argument(
        "--point-file",
        type=Path,
        default=None,
        help="采样点文件，每行一个点，支持空格或逗号分隔",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-8,
        help="坐标匹配容差，默认 1e-8",
    )
    parser.add_argument(
        "--time-tol",
        type=float,
        default=1e-12,
        help="时间标签不一致 warning 阈值，默认 1e-12；比较仍按列序号进行",
    )
    parser.add_argument(
        "--rel-floor",
        type=float,
        default=1e-30,
        help="相对误差分母下限，默认 1e-30",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="写出逐点逐时间步误差明细 CSV",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=None,
        help="输出 FEM 与 COMSOL 时间曲线对比图片目录；每个采样点/字段一张图",
    )
    parser.add_argument(
        "--plot-format",
        choices=("svg", "png"),
        default="svg",
        help="图片格式，默认 svg；png 需要 matplotlib",
    )
    parser.add_argument(
        "--plot-dpi",
        type=int,
        default=160,
        help="PNG 输出 DPI，默认 160",
    )
    parser.add_argument(
        "--no-point-summary",
        action="store_true",
        help="不打印逐采样点指标，只打印字段聚合指标",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.tol <= 0.0:
        print("错误: --tol 必须 > 0", file=sys.stderr)
        return 1
    if args.fields_per_step <= 0:
        print("错误: --fields-per-step 必须 > 0", file=sys.stderr)
        return 1
    if args.rel_floor <= 0.0:
        print("错误: --rel-floor 必须 > 0", file=sys.stderr)
        return 1
    if args.plot_dpi <= 0:
        print("错误: --plot-dpi 必须 > 0", file=sys.stderr)
        return 1

    try:
        field_names = parse_field_names(args.field_names, args.fields_per_step)
        selected_fields = parse_field_selection(args.fields, field_names)

        mine_path = resolve_existing_path(args.mine)
        comsol_path = resolve_existing_path(args.comsol)
        mine = load_dataset(mine_path, "mine", args.coord_dim)
        comsol = load_dataset(comsol_path, "comsol", args.coord_dim)
    except Exception as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    mine_vcols = min(len(row) for row in mine.values)
    comsol_vcols = min(len(row) for row in comsol.values)
    n_steps_mine = mine_vcols // args.fields_per_step
    n_steps_comsol = comsol_vcols // args.fields_per_step
    n_steps = min(n_steps_mine, n_steps_comsol)

    print("==== 数据摘要 ====")
    print(f"[mine] path={mine.path}, rows={len(mine.coords)}, value_cols={mine_vcols}")
    print(f"[comsol] path={comsol.path}, rows={len(comsol.coords)}, value_cols={comsol_vcols}")
    print(
        f"fields_per_step={args.fields_per_step}, fields={','.join(field_names)}, "
        f"selected={','.join(field_names[i] for i in selected_fields)}"
    )

    if n_steps <= 0:
        print("错误: 无法从值列数量推断时间步数。", file=sys.stderr)
        return 1

    if mine_vcols % args.fields_per_step != 0:
        print(f"警告: mine 值列数 {mine_vcols} 不能被 fields_per_step 整除，尾列将被忽略。")
    if comsol_vcols % args.fields_per_step != 0:
        print(f"警告: comsol 值列数 {comsol_vcols} 不能被 fields_per_step 整除，尾列将被忽略。")
    if n_steps_mine != n_steps_comsol:
        print(f"警告: 时间步数不同，mine={n_steps_mine}, comsol={n_steps_comsol}，只比较前 {n_steps} 步。")

    mine_times = infer_times(mine, args.fields_per_step, n_steps)
    comsol_times = infer_times(comsol, args.fields_per_step, n_steps)
    max_time_diff = max(abs(a - b) for a, b in zip(mine_times, comsol_times))
    print(f"time_steps={n_steps}, t0=({mine_times[0]:.12g}, {comsol_times[0]:.12g}), "
          f"t_last=({mine_times[-1]:.12g}, {comsol_times[-1]:.12g})")
    if max_time_diff > args.time_tol:
        print(
            f"警告: mine 与 COMSOL 时间标签最大差异 {max_time_diff:.12g} > "
            f"--time-tol {args.time_tol:.12g}，当前仍按列序号比较。"
        )

    try:
        explicit_points = parse_points(args.points)
        if args.point_file is not None:
            explicit_points.extend(load_point_file(args.point_file))
    except Exception as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    if explicit_points:
        matches, missed = match_requested_points(explicit_points, mine, comsol, args.tol)
        if missed:
            print(f"警告: {len(missed)} 个指定点未在两份数据中同时匹配到。")
            for p in missed[:10]:
                print(f"  missed {fmt_coord(p)}")
            if len(missed) > 10:
                print(f"  ... 其余 {len(missed) - 10} 个省略")
    else:
        all_matches = match_all_coords(mine, comsol, args.tol)
        if not all_matches:
            print(f"错误: 未匹配到任何坐标点，建议调大 --tol（当前 {args.tol:g}）。", file=sys.stderr)
            return 1

        sample_size = args.random_points if args.random_points is not None else 5
        if sample_size <= 0:
            print("错误: --random-points 必须 > 0", file=sys.stderr)
            return 1
        if sample_size > len(all_matches):
            print(f"警告: 请求 {sample_size} 个随机点，但只匹配到 {len(all_matches)} 个点，将全部使用。")
            sample_size = len(all_matches)
        rng = random.Random(args.seed)
        matches = rng.sample(all_matches, sample_size)
        print(f"坐标匹配: {len(all_matches)}/{len(mine.coords)}，随机采样 {len(matches)} 个点，seed={args.seed}")

    if not matches:
        print("错误: 没有可比较的采样点。", file=sys.stderr)
        return 1

    return compare_curves(
        mine=mine,
        comsol=comsol,
        matches=matches,
        fields_per_step=args.fields_per_step,
        field_names=field_names,
        selected_fields=selected_fields,
        n_steps=n_steps,
        mine_times=mine_times,
        comsol_times=comsol_times,
        rel_floor=args.rel_floor,
        output_csv=args.output_csv,
        plot_dir=args.plot_dir,
        plot_format=args.plot_format,
        plot_dpi=args.plot_dpi,
        print_points=not args.no_point_summary,
    )


if __name__ == "__main__":
    raise SystemExit(main())
