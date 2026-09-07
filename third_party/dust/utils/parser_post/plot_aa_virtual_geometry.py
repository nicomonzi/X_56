#!/usr/bin/env python3
"""Quick visualization for aeroacoustic DAT geometry (base + virtual panels).

This script reads an aeroacoustic .dat file written by DUST postprocessing and
plots:
1) Base element centers (VL/surface centers)
2) Virtual upper/lower panel centers (if present in the new extended format)

It is intended as a lightweight sanity check for the new acoustic post output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


LEGACY_COLS = 15
EXTENDED_COLS = 35


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot geometry from DUST aeroacoustic .dat output"
    )
    parser.add_argument("dat_file", type=Path, help="Path to aeroacoustic .dat file")
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Use one point every N to reduce clutter (default: 1)",
    )
    parser.add_argument(
        "--normal-scale",
        type=float,
        default=0.2,
        help="Arrow length for normal vectors (default: 0.2)",
    )
    parser.add_argument(
        "--show-normals",
        action="store_true",
        help="Also draw normal vectors for sampled points",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="Save figure to this path instead of showing it",
    )
    parser.add_argument(
        "--elev",
        type=float,
        default=22.0,
        help="3D camera elevation angle in degrees (default: 22)",
    )
    parser.add_argument(
        "--azim",
        type=float,
        default=-55.0,
        help="3D camera azimuth angle in degrees (default: -55)",
    )
    return parser.parse_args()


def load_aa_data(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    with path.open("r", encoding="utf-8") as fobj:
        for line in fobj:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            vals = np.fromstring(stripped, sep=" ")
            if vals.size:
                rows.append(vals)

    if len(rows) < 2:
        raise ValueError(
            "File does not contain enough numeric rows. Expected metadata + data rows."
        )

    meta = rows[0]
    data_rows = rows[1:]
    ncols = data_rows[0].size

    for i_row, row in enumerate(data_rows, start=1):
        if row.size != ncols:
            raise ValueError(
                f"Inconsistent column count at data row {i_row}: "
                f"expected {ncols}, got {row.size}."
            )

    data = np.vstack(data_rows)
    return meta, data


def set_equal_3d(ax: plt.Axes, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> None:
    x_mid = 0.5 * (np.min(x) + np.max(x))
    y_mid = 0.5 * (np.min(y) + np.max(y))
    z_mid = 0.5 * (np.min(z) + np.max(z))
    max_range = max(np.ptp(x), np.ptp(y), np.ptp(z))
    if max_range <= 0.0:
        max_range = 1.0

    half = 0.5 * max_range
    ax.set_xlim(x_mid - half, x_mid + half)
    ax.set_ylim(y_mid - half, y_mid + half)
    ax.set_zlim(z_mid - half, z_mid + half)


def main() -> None:
    args = parse_args()

    if args.stride < 1:
        raise ValueError("--stride must be >= 1")

    meta, data = load_aa_data(args.dat_file)
    n_elem = data.shape[0]
    ncols = data.shape[1]

    if ncols < LEGACY_COLS:
        raise ValueError(
            f"Unsupported format: expected at least {LEGACY_COLS} columns, got {ncols}."
        )

    # Base element fields
    cen = data[:, 0:3]
    nor = data[:, 3:6]

    has_virtual = ncols >= EXTENDED_COLS

    if has_virtual:
        cen_u = data[:, 15:18]
        cen_l = data[:, 18:21]
        nor_u = data[:, 21:24]
        nor_l = data[:, 24:27]
        area_u = data[:, 27]
        area_l = data[:, 28]

        mask_u = area_u > 0.0
        mask_l = area_l > 0.0
    else:
        cen_u = cen_l = nor_u = nor_l = np.empty((0, 3))
        mask_u = mask_l = np.zeros(0, dtype=bool)

    sampled = np.arange(0, n_elem, args.stride)

    print(f"File: {args.dat_file}")
    print(f"Metadata values found: {meta.size}")
    if meta.size >= 1:
        print(f"Header num_el: {int(meta[0])}")
    print(f"Data rows: {n_elem}")
    print(f"Columns per row: {ncols}")
    print(f"New extended format detected ({EXTENDED_COLS} cols): {has_virtual}")
    if has_virtual:
        print(f"Rows with valid upper virtual panel: {np.count_nonzero(mask_u)}")
        print(f"Rows with valid lower virtual panel: {np.count_nonzero(mask_l)}")

    fig = plt.figure(figsize=(10, 9))
    ax3d = fig.add_subplot(1, 1, 1, projection="3d")

    # 3D view
    ax3d.scatter(
        cen[sampled, 0],
        cen[sampled, 1],
        cen[sampled, 2],
        s=14,
        c="0.15",
        alpha=0.82,
        depthshade=False,
        label="Base centers",
    )

    if has_virtual:
        idx_u = sampled[mask_u[sampled]]
        idx_l = sampled[mask_l[sampled]]

        ax3d.scatter(
            cen_u[idx_u, 0],
            cen_u[idx_u, 1],
            cen_u[idx_u, 2],
            s=20,
            c="tab:red",
            alpha=0.9,
            depthshade=False,
            label="Upper virtual",
        )
        ax3d.scatter(
            cen_l[idx_l, 0],
            cen_l[idx_l, 1],
            cen_l[idx_l, 2],
            s=20,
            c="tab:blue",
            alpha=0.9,
            depthshade=False,
            label="Lower virtual",
        )

    if args.show_normals:
        ax3d.quiver(
            cen[sampled, 0],
            cen[sampled, 1],
            cen[sampled, 2],
            nor[sampled, 0],
            nor[sampled, 1],
            nor[sampled, 2],
            length=args.normal_scale,
            normalize=True,
            color="0.35",
            linewidth=0.6,
            alpha=0.55,
        )

        if has_virtual:
            idx_u = sampled[mask_u[sampled]]
            idx_l = sampled[mask_l[sampled]]

            ax3d.quiver(
                cen_u[idx_u, 0],
                cen_u[idx_u, 1],
                cen_u[idx_u, 2],
                nor_u[idx_u, 0],
                nor_u[idx_u, 1],
                nor_u[idx_u, 2],
                length=args.normal_scale,
                normalize=True,
                color="tab:red",
                linewidth=0.6,
                alpha=0.45,
            )
            ax3d.quiver(
                cen_l[idx_l, 0],
                cen_l[idx_l, 1],
                cen_l[idx_l, 2],
                nor_l[idx_l, 0],
                nor_l[idx_l, 1],
                nor_l[idx_l, 2],
                length=args.normal_scale,
                normalize=True,
                color="tab:blue",
                linewidth=0.6,
                alpha=0.45,
            )

    ax3d.set_title("3D geometry check")
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.view_init(elev=args.elev, azim=args.azim)
    ax3d.grid(True, alpha=0.35)
    ax3d.set_facecolor("#f5f7fa")
    ax3d.xaxis.pane.set_facecolor((0.95, 0.96, 0.98, 1.0))
    ax3d.yaxis.pane.set_facecolor((0.95, 0.96, 0.98, 1.0))
    ax3d.zaxis.pane.set_facecolor((0.95, 0.96, 0.98, 1.0))

    all_x = [cen[:, 0]]
    all_y = [cen[:, 1]]
    all_z = [cen[:, 2]]
    if has_virtual:
        if np.any(mask_u):
            all_x.append(cen_u[mask_u, 0])
            all_y.append(cen_u[mask_u, 1])
            all_z.append(cen_u[mask_u, 2])
        if np.any(mask_l):
            all_x.append(cen_l[mask_l, 0])
            all_y.append(cen_l[mask_l, 1])
            all_z.append(cen_l[mask_l, 2])
    set_equal_3d(ax3d, np.concatenate(all_x), np.concatenate(all_y), np.concatenate(all_z))
    ax3d.legend(loc="upper left", frameon=True, framealpha=0.95)

    fig.suptitle(args.dat_file.name)
    fig.tight_layout(pad=1.2)

    if args.save is not None:
        fig.savefig(args.save, dpi=180)
        print(f"Saved plot to: {args.save}")
    else:
        backend = plt.get_backend().lower()
        if "agg" in backend:
            auto_out = args.dat_file.with_suffix(".png")
            fig.savefig(auto_out, dpi=180)
            print(
                "Non-interactive matplotlib backend detected "
                f"('{plt.get_backend()}'): saved plot to {auto_out}"
            )
        else:
            plt.show()


if __name__ == "__main__":
    main()
