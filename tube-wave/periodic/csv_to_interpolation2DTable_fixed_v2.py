#!/usr/bin/env python3
"""
Convert a CSV slice export to an OpenFOAM interpolation2DTable.

This version fixes the geometry-mask problem seen with circular slices.
The previous envelope-based mask could fail on unstructured point clouds:
if a few rounded y-bins were sparse or one-sided, interpolating z_min/z_max
between those bins produced artificial gaps (for example, only positive z
values on a row that should span both negative and positive z).

For circular / pipe slices, this script can infer a circle from the convex
hull of the original point cloud and use that circle as the mask. This keeps
all rows consistent with the actual circular slice shape.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import ConvexHull

DEFAULT_Y_CANDIDATES = ["Points:1", "y", "Y", "coordY"]
DEFAULT_Z_CANDIDATES = ["Points:2", "z", "Z", "coordZ"]
DEFAULT_VALUE_CANDIDATES = ["U:0", "Ux", "ux", "u_x", "value", "Value"]


def infer_column(df: pd.DataFrame, candidates: list[str], explicit_name: Optional[str], explicit_index: Optional[int], kind: str) -> pd.Series:
    if explicit_name is not None:
        if explicit_name not in df.columns:
            raise KeyError(f"Requested {kind} column name '{explicit_name}' not found. Available columns: {list(df.columns)}")
        return pd.to_numeric(df[explicit_name], errors="coerce")

    if explicit_index is not None:
        if explicit_index < 0 or explicit_index >= len(df.columns):
            raise IndexError(f"Requested {kind} column index {explicit_index} out of range for {len(df.columns)} columns")
        return pd.to_numeric(df.iloc[:, explicit_index], errors="coerce")

    for name in candidates:
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce")

    raise KeyError(
        f"Could not infer {kind} column automatically. Available columns: {list(df.columns)}. "
        f"Please pass --{kind}-name or --{kind}-col explicitly."
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert slice CSV to OpenFOAM interpolation2DTable")
    p.add_argument("--input", required=True, help="Input CSV file")
    p.add_argument("--output", required=True, help="Output interpolation2DTable file")

    p.add_argument("--y-name", default=None, help="Column name for y")
    p.add_argument("--z-name", default=None, help="Column name for z")
    p.add_argument("--value-name", default=None, help="Column name for the table value")
    p.add_argument("--y-col", type=int, default=None, help="0-based column index for y")
    p.add_argument("--z-col", type=int, default=None, help="0-based column index for z")
    p.add_argument("--value-col", type=int, default=None, help="0-based column index for value")

    p.add_argument("--dy", type=float, default=0.001, help="Regular output spacing in y")
    p.add_argument("--dz", type=float, default=0.001, help="Regular output spacing in z")
    p.add_argument("--y-min", type=float, default=None, help="Optional minimum y for output grid")
    p.add_argument("--y-max", type=float, default=None, help="Optional maximum y for output grid")
    p.add_argument("--z-min", type=float, default=None, help="Optional minimum z for output grid")
    p.add_argument("--z-max", type=float, default=None, help="Optional maximum z for output grid")

    p.add_argument("--fill-nearest", action="store_true", help="Fill NaNs left by linear interpolation using nearest neighbor")
    p.add_argument(
        "--mask-method",
        choices=["auto", "circle", "envelope", "none"],
        default="auto",
        help="Geometry mask used on the regular grid",
    )
    p.add_argument("--boundary-round-decimals", type=int, default=4, help="Rounding used by envelope mask")
    p.add_argument("--mask-padding", type=float, default=None, help="Extra half-width added to each valid z-range. Default = 0.5*dz")

    # Circular mask options
    p.add_argument("--circle-center-y", type=float, default=None, help="Explicit circle center y")
    p.add_argument("--circle-center-z", type=float, default=None, help="Explicit circle center z")
    p.add_argument("--circle-radius", type=float, default=None, help="Explicit circle radius")
    p.add_argument(
        "--circle-fit-tol-rel",
        type=float,
        default=0.02,
        help="Relative RMS residual threshold for accepting auto circle fit from convex hull",
    )

    p.add_argument("--fmt", default=".8g", help="Number format, e.g. .8g or .6f")
    return p.parse_args()


def read_input(args: argparse.Namespace) -> pd.DataFrame:
    raw = pd.read_csv(args.input)
    y = infer_column(raw, DEFAULT_Y_CANDIDATES, args.y_name, args.y_col, "y")
    z = infer_column(raw, DEFAULT_Z_CANDIDATES, args.z_name, args.z_col, "z")
    value = infer_column(raw, DEFAULT_VALUE_CANDIDATES, args.value_name, args.value_col, "value")

    out = pd.DataFrame({"y": y, "z": z, "value": value}).dropna().copy()
    out = out.sort_values(["y", "z"]).reset_index(drop=True)
    return out


def make_regular_axis(vmin: float, vmax: float, dv: float) -> np.ndarray:
    tol = 1.0e-9 * max(1.0, abs(vmin), abs(vmax), abs(dv))
    qmin = vmin / dv
    qmax = vmax / dv
    nmin = np.round(qmin)
    nmax = np.round(qmax)
    start = nmin * dv if abs(qmin - nmin) < tol / abs(dv) else np.ceil(qmin) * dv
    end = nmax * dv if abs(qmax - nmax) < tol / abs(dv) else np.floor(qmax) * dv
    start = np.round(start, 12)
    end = np.round(end, 12)
    n = int(np.floor((end - start) / dv + 0.5)) + 1
    return np.round(start + np.arange(n) * dv, 12)


def fit_circle_least_squares(points: np.ndarray) -> tuple[float, float, float, float]:
    x = points[:, 0]
    y = points[:, 1]
    A = np.column_stack([x, y, np.ones_like(x)])
    b = -(x * x + y * y)
    coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    a, bb, c = coef
    cy = -a / 2.0
    cz = -bb / 2.0
    radius = float(np.sqrt(max(cy * cy + cz * cz - c, 0.0)))
    r = np.sqrt((x - cy) ** 2 + (y - cz) ** 2)
    rms = float(np.sqrt(np.mean((r - radius) ** 2)))
    return float(cy), float(cz), float(radius), rms


def infer_circle_from_hull(df: pd.DataFrame) -> tuple[float, float, float, float]:
    pts = df[["y", "z"]].to_numpy(dtype=float)
    hull = ConvexHull(pts)
    hull_points = pts[hull.vertices]
    return fit_circle_least_squares(hull_points)


def build_geometry_envelope(df: pd.DataFrame, y_grid: np.ndarray, round_decimals: int) -> tuple[np.ndarray, np.ndarray]:
    tmp = df.copy()
    tmp["y_bin"] = tmp["y"].round(round_decimals)
    boundary = (
        tmp.groupby("y_bin", as_index=False)
        .agg(z_min=("z", "min"), z_max=("z", "max"))
        .sort_values("y_bin")
    )
    yb = boundary["y_bin"].to_numpy(dtype=float)
    zmin = boundary["z_min"].to_numpy(dtype=float)
    zmax = boundary["z_max"].to_numpy(dtype=float)
    if len(yb) < 2:
        raise RuntimeError("Not enough data to build an envelope mask.")
    zmin_grid = np.interp(y_grid, yb, zmin)
    zmax_grid = np.interp(y_grid, yb, zmax)
    return zmin_grid, zmax_grid


def resolve_mask(df: pd.DataFrame, y_grid: np.ndarray, args: argparse.Namespace) -> tuple[str, np.ndarray, np.ndarray, dict]:
    pad = 0.5 * args.dz if args.mask_padding is None else args.mask_padding

    if args.mask_method == "none":
        zmin_grid = np.full_like(y_grid, float(df["z"].min()))
        zmax_grid = np.full_like(y_grid, float(df["z"].max()))
        return "none", zmin_grid, zmax_grid, {"pad": pad}

    if args.mask_method == "envelope":
        zmin_grid, zmax_grid = build_geometry_envelope(df, y_grid, args.boundary_round_decimals)
        return "envelope", zmin_grid, zmax_grid, {"pad": pad}

    # circle or auto
    if args.circle_center_y is not None and args.circle_center_z is not None and args.circle_radius is not None:
        cy = float(args.circle_center_y)
        cz = float(args.circle_center_z)
        radius = float(args.circle_radius)
        rms = 0.0
        accepted = True
    else:
        cy, cz, radius, rms = infer_circle_from_hull(df)
        accepted = (rms / max(radius, 1e-30)) <= args.circle_fit_tol_rel

    if args.mask_method == "circle" and not accepted:
        raise RuntimeError(
            f"Circle fit from convex hull was not sufficiently accurate: RMS={rms:.6g}, R={radius:.6g}, rel={rms/max(radius, 1e-30):.6g}"
        )

    if args.mask_method == "circle" or (args.mask_method == "auto" and accepted):
        dy = y_grid - cy
        inside = radius * radius - dy * dy
        inside = np.where(inside >= 0.0, inside, np.nan)
        half_width = np.sqrt(inside)
        zmin_grid = cz - half_width
        zmax_grid = cz + half_width
        return "circle", zmin_grid, zmax_grid, {"pad": pad, "cy": cy, "cz": cz, "radius": radius, "rms": rms}

    zmin_grid, zmax_grid = build_geometry_envelope(df, y_grid, args.boundary_round_decimals)
    return "envelope", zmin_grid, zmax_grid, {"pad": pad, "circle_rms": rms, "circle_radius": radius}


def build_table(df: pd.DataFrame, args: argparse.Namespace) -> tuple[list[tuple[float, list[tuple[float, float]]]], tuple[str, dict]]:
    y_min = float(df["y"].min()) if args.y_min is None else args.y_min
    y_max = float(df["y"].max()) if args.y_max is None else args.y_max
    z_min = float(df["z"].min()) if args.z_min is None else args.z_min
    z_max = float(df["z"].max()) if args.z_max is None else args.z_max

    y_grid = make_regular_axis(y_min, y_max, args.dy)
    z_grid = make_regular_axis(z_min, z_max, args.dz)

    points = df[["y", "z"]].to_numpy(dtype=float)
    values = df["value"].to_numpy(dtype=float)
    linear_interp = LinearNDInterpolator(points, values, fill_value=np.nan)
    nearest_interp = NearestNDInterpolator(points, values) if args.fill_nearest else None

    mask_name, zmin_grid, zmax_grid, mask_info = resolve_mask(df, y_grid, args)
    pad = mask_info["pad"]

    table: list[tuple[float, list[tuple[float, float]]]] = []

    for i, yv in enumerate(y_grid):
        zlo = zmin_grid[i]
        zhi = zmax_grid[i]
        if not np.isfinite(zlo) or not np.isfinite(zhi):
            continue
        zlo -= pad
        zhi += pad
        mask = (z_grid >= zlo) & (z_grid <= zhi)
        if np.count_nonzero(mask) < 2:
            continue

        z_row = z_grid[mask]
        q = np.column_stack([np.full(z_row.shape, yv), z_row])
        vals = linear_interp(q)

        if nearest_interp is not None:
            nan_mask = np.isnan(vals)
            if np.any(nan_mask):
                vals[nan_mask] = nearest_interp(q[nan_mask])

        valid = ~np.isnan(vals)
        if np.count_nonzero(valid) < 2:
            continue

        row = [(float(z), float(v)) for z, v in zip(z_row[valid], vals[valid])]
        table.append((float(yv), row))

    if not table:
        raise RuntimeError("No table rows were generated.")

    return table, (mask_name, mask_info)


def write_table(path: str, table: list[tuple[float, list[tuple[float, float]]]], number_fmt: str) -> None:
    def fmt(x: float) -> str:
        return format(float(x), number_fmt)

    with open(path, "w", encoding="utf-8") as f:
        f.write("// -*- C++ -*-\n")
        f.write("(\n")
        for yv, row in table:
            f.write("    (\n")
            f.write(f"        {fmt(yv)}\n")
            f.write("        (\n")
            for zv, val in row:
                f.write(f"            ({fmt(zv)} {fmt(val)})\n")
            f.write("        )\n")
            f.write("    )\n")
        f.write(")\n")


def main() -> None:
    args = parse_args()
    df = read_input(args)
    table, mask_meta = build_table(df, args)
    write_table(args.output, table, args.fmt)

    counts = [len(row) for _, row in table]
    mask_name, info = mask_meta
    print(f"Wrote {len(table)} y-slices to: {Path(args.output).resolve()}")
    print(f"Mask method used: {mask_name}")
    if mask_name == "circle":
        print(
            f"Circle: center=({info['cy']:.8g}, {info['cz']:.8g}), R={info['radius']:.8g}, RMS={info['rms']:.3g}"
        )
    print(f"Min/avg/max points per y-slice: {min(counts)} / {sum(counts)/len(counts):.1f} / {max(counts)}")
    print(f"y-range: {table[0][0]} to {table[-1][0]}")


if __name__ == "__main__":
    main()
