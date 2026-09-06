#!/usr/bin/env python3
"""Plot station-local MBDyn--Nastran errors with upstream accumulation removed."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler


STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY))
from analyze_convergence import (  # noqa: E402
    RIGHT_SEMISPAN_NODES,
    read_mbdyn_semispan,
    read_nastran_displacements,
)


THESIS_COLORS = ["#8B0000", "#00008B", "#66B2FF", "#006400", "#CC5500"]
FONT_SIZE = 15
DISPLAYED_MODES = (1, 8, 25, 54)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.titlesize": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "axes.prop_cycle": cycler(color=THESIS_COLORS),
})


def save_figure(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(output / f"{stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--f06", type=Path,
        default=STUDY / "nastran/MAIN/sol101_gravity_5g.f06",
    )
    parser.add_argument("--results", type=Path, default=STUDY / "mbdyn/results")
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).resolve().parent / "plots"
    )
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    nastran = read_nastran_displacements(args.f06, RIGHT_SEMISPAN_NODES)
    nastran_field = np.asarray(
        [nastran[node] for node in RIGHT_SEMISPAN_NODES], dtype=float
    )
    cases = []
    for modes in DISPLAYED_MODES:
        result = args.results / f"gravity_5g_{modes:02d}_elastic_modes.nc"
        if not result.is_file():
            raise FileNotFoundError(f"Missing MBDyn result: {result}")
        b, displacement, rotation = read_mbdyn_semispan(
            result, RIGHT_SEMISPAN_NODES, args.tail_fraction
        )
        mbdyn_field = np.column_stack((displacement, rotation))
        db = np.diff(b)
        if np.any(db <= 0.0):
            raise RuntimeError("Semispan coordinates must be strictly increasing")

        # Translation curvature removes both accumulated displacement and
        # accumulated slope. Rotation gradient removes accumulated angle.
        b_mid = 0.5 * (b[:-1] + b[1:])
        db_mid = np.diff(b_mid)
        b_curvature = 0.5 * (b_mid[:-1] + b_mid[1:])
        mbdyn_translation_gradient = np.diff(mbdyn_field[:, :3], axis=0) / db[:, None]
        nastran_translation_gradient = np.diff(nastran_field[:, :3], axis=0) / db[:, None]
        mbdyn_translation_curvature = np.diff(
            mbdyn_translation_gradient, axis=0
        ) / db_mid[:, None]
        nastran_translation_curvature = np.diff(
            nastran_translation_gradient, axis=0
        ) / db_mid[:, None]
        translation_local_error = np.abs(
            mbdyn_translation_curvature - nastran_translation_curvature
        )

        mbdyn_rotation_gradient = np.diff(mbdyn_field[:, 3:], axis=0) / db[:, None]
        nastran_rotation_gradient = np.diff(nastran_field[:, 3:], axis=0) / db[:, None]
        rotation_local_error = np.abs(
            mbdyn_rotation_gradient - nastran_rotation_gradient
        )
        rotation_local_at_curvature = 0.5 * (
            rotation_local_error[:-1] + rotation_local_error[1:]
        )

        translation_reference = np.linalg.norm(
            nastran_translation_curvature, axis=1
        )
        rotation_reference_at_curvature = np.linalg.norm(
            0.5 * (
                nastran_rotation_gradient[:-1] + nastran_rotation_gradient[1:]
            ),
            axis=1,
        )
        translation_error = np.linalg.norm(translation_local_error, axis=1)
        rotation_error = np.linalg.norm(rotation_local_at_curvature, axis=1)
        translation_scale = max(float(np.max(translation_reference)), 1.0e-15)
        rotation_scale = max(
            float(np.max(rotation_reference_at_curvature)), 1.0e-15
        )
        local_global_percent = 100.0 * np.sqrt(
            0.5
            * (
                np.square(translation_error / translation_scale)
                + np.square(rotation_error / rotation_scale)
            )
        )
        cases.append(
            {
                "modes": modes,
                "b_curvature": b_curvature,
                "b_rotation": b_mid,
                "translation_local": translation_local_error,
                "rotation_local": rotation_local_error,
                "global_percent": local_global_percent,
            }
        )

    fields = [
        "elastic_modes", "quantity", "component", "station", "b_in",
        "local_error", "unit",
    ]
    csv_path = args.output / "gravity_5g_local_spanwise_errors.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for case in cases:
            for component, name in enumerate(("Ux", "Uy", "Uz")):
                for station, value in enumerate(case["translation_local"][:, component]):
                    writer.writerow({
                        "elastic_modes": case["modes"], "quantity": "translation_curvature",
                        "component": name, "station": station + 1,
                        "b_in": case["b_curvature"][station], "local_error": value,
                        "unit": "1/in",
                    })
            for component, name in enumerate(("Rx", "Ry", "Rz")):
                for station, value in enumerate(case["rotation_local"][:, component]):
                    writer.writerow({
                        "elastic_modes": case["modes"], "quantity": "rotation_gradient",
                        "component": name, "station": station + 1,
                        "b_in": case["b_rotation"][station], "local_error": value,
                        "unit": "rad/in",
                    })
            for station, value in enumerate(case["global_percent"]):
                writer.writerow({
                    "elastic_modes": case["modes"], "quantity": "normalized_global",
                    "component": "all", "station": station + 1,
                    "b_in": case["b_curvature"][station], "local_error": value,
                    "unit": "%",
                })

    titles = (
        "Displacement in x direction local curvature error",
        "Displacement in y direction local curvature error",
        "Displacement in z direction local curvature error",
        "Rotation about x axis local curvature error",
        "Rotation about y axis local curvature error",
        "Rotation about z axis local curvature error",
    )
    names = ("Ux", "Uy", "Uz", "Rx", "Ry", "Rz")
    symbols = ("u_x", "u_y", "u_z", r"\theta_x", r"\theta_y", r"\theta_z")
    markers = ("o", "s", "^", "D")
    for component, (title, name, symbol) in enumerate(zip(titles, names, symbols)):
        fig, axis = plt.subplots(figsize=(11.5, 7.0))
        for index, case in enumerate(cases):
            modes = int(case["modes"])
            if component < 3:
                plot_b = case["b_curvature"]
                plot_error = case["translation_local"][:, component]
            else:
                plot_b = case["b_rotation"]
                plot_error = case["rotation_local"][:, component - 3]
            axis.plot(
                plot_b,
                plot_error,
                color=THESIS_COLORS[index],
                marker=markers[index],
                markersize=4.5,
                linewidth=1.8,
                linestyle="-",
                label="1 elastic mode" if modes == 1 else f"{modes} elastic modes",
            )
        axis.set_title(title, fontweight="normal")
        axis.set_xlabel("b [in]")
        if component < 3:
            axis.set_ylabel(rf"Local $|\partial^2 {symbol}/\partial b^2|$ [1/in]")
        else:
            axis.set_ylabel(rf"Local $|\partial {symbol}/\partial b|$ [rad/in]")
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=2, frameon=False)
        save_figure(fig, args.output, f"gravity_5g_local_spanwise_error_{name}")

    fig, axis = plt.subplots(figsize=(11.5, 7.0))
    for index, case in enumerate(cases):
        modes = int(case["modes"])
        axis.plot(
            case["b_curvature"],
            case["global_percent"],
            color=THESIS_COLORS[index],
            marker=markers[index],
            markersize=4.5,
            linewidth=1.8,
            linestyle="-",
            label="1 elastic mode" if modes == 1 else f"{modes} elastic modes",
        )
    axis.set_title("Local normalized spanwise error", fontweight="normal")
    axis.set_xlabel("b [in]")
    axis.set_ylabel("Local normalized error [%]")
    axis.grid(True, alpha=0.3)
    axis.legend(ncol=2, frameon=False)
    save_figure(fig, args.output, "gravity_5g_local_normalized_spanwise_error")

    print(
        f"Used {len(RIGHT_SEMISPAN_NODES) - 1} rotation-gradient segments and "
        f"{len(RIGHT_SEMISPAN_NODES) - 2} translation-curvature stations."
    )
    print(f"Created 7 local-error plots in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
