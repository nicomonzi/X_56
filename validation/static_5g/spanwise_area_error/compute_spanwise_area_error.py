#!/usr/bin/env python3
"""Integrate a dimensionless MBDyn--Nastran error over the right semispan."""

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


def span_length_percentage(length_error: np.ndarray, b: np.ndarray) -> float:
    """Area-averaged length error as a percentage of semispan length."""
    length = float(b[-1] - b[0])
    return 100.0 * float(np.trapz(length_error, b) / length**2)


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

    first_result = args.results / "gravity_5g_01_elastic_modes.nc"
    b, _, _ = read_mbdyn_semispan(
        first_result, RIGHT_SEMISPAN_NODES, args.tail_fraction
    )
    length = float(b[-1] - b[0])

    summary_rows = []
    nodal_rows = []
    for modes in range(1, 55):
        result = args.results / f"gravity_5g_{modes:02d}_elastic_modes.nc"
        if not result.is_file():
            raise FileNotFoundError(f"Missing MBDyn result: {result}")
        case_b, displacement, rotation = read_mbdyn_semispan(
            result, RIGHT_SEMISPAN_NODES, args.tail_fraction
        )
        if not np.allclose(case_b, b, rtol=0.0, atol=1.0e-9):
            raise RuntimeError(f"Inconsistent span coordinates in {result.name}")

        translation_difference = displacement - nastran_field[:, :3]
        rotation_difference = rotation - nastran_field[:, 3:]
        translation_nodal_in = np.linalg.norm(translation_difference, axis=1)
        rotation_equivalent_nodal_in = b * np.linalg.norm(
            rotation_difference, axis=1
        )
        combined_nodal_in = np.sqrt(
            np.square(translation_nodal_in)
            + np.square(rotation_equivalent_nodal_in)
        )

        component_area = {}
        for component, name in enumerate(("Ux", "Uy", "Uz")):
            component_area[f"{name}_span_area_error_percent"] = span_length_percentage(
                np.abs(translation_difference[:, component]),
                b,
            )
        for component, name in enumerate(("Rx", "Ry", "Rz")):
            component_area[f"{name}_equivalent_span_area_error_percent"] = span_length_percentage(
                b * np.abs(rotation_difference[:, component]),
                b,
            )

        translation_area = span_length_percentage(translation_nodal_in, b)
        rotation_area = span_length_percentage(rotation_equivalent_nodal_in, b)
        combined_area = span_length_percentage(combined_nodal_in, b)
        summary_rows.append({
            "elastic_modes": modes,
            "spanwise_nodes": len(RIGHT_SEMISPAN_NODES),
            "semispan_length_in": length,
            "translation_span_area_error_percent": translation_area,
            "rotation_equivalent_span_area_error_percent": rotation_area,
            "combined_span_area_error_percent": combined_area,
            **component_area,
        })
        for index, node in enumerate(RIGHT_SEMISPAN_NODES):
            nodal_rows.append({
                "elastic_modes": modes,
                "node": node,
                "b_in": b[index],
                "translation_error_in": translation_nodal_in[index],
                "rotation_equivalent_error_in": rotation_equivalent_nodal_in[index],
                "combined_equivalent_error_in": combined_nodal_in[index],
                "combined_error_percent_of_span": (
                    100.0 * combined_nodal_in[index] / length
                ),
            })

    summary_path = args.output / "gravity_5g_spanwise_area_error.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    nodal_path = args.output / "gravity_5g_nodal_composite_error.csv"
    with nodal_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=nodal_rows[0].keys())
        writer.writeheader()
        writer.writerows(nodal_rows)

    modes = np.asarray([row["elastic_modes"] for row in summary_rows], dtype=int)
    combined = np.asarray(
        [row["combined_span_area_error_percent"] for row in summary_rows], dtype=float
    )
    fig, axis = plt.subplots(figsize=(9.2, 6.0))
    axis.plot(
        modes,
        combined,
        color=THESIS_COLORS[1],
        marker="o",
        markersize=5.5,
        linewidth=2.1,
        linestyle="-",
        label=f"{len(RIGHT_SEMISPAN_NODES)} spanwise nodes",
    )
    axis.set_title("Span-normalized integrated error", fontweight="normal")
    axis.set_xlabel("Number of elastic modes, N")
    axis.set_ylabel("Cumulative area error relative to span [%]")
    axis.set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54])
    axis.set_xlim(0, 56)
    error_range = float(np.max(combined) - np.min(combined))
    margin = max(0.08 * error_range, 0.005)
    axis.set_ylim(
        float(np.min(combined)) - margin,
        float(np.max(combined)) + margin,
    )
    axis.grid(False)
    axis.yaxis.grid(True, which="major", color="0.82", linewidth=0.7, alpha=0.65)
    axis.legend(frameon=False)
    axis.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(
        args.output / "gravity_5g_spanwise_integrated_normalized_error.png",
        dpi=450,
        bbox_inches="tight",
    )
    fig.savefig(
        args.output / "gravity_5g_spanwise_integrated_normalized_error.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    translation_combined = np.asarray(
        [row["translation_span_area_error_percent"] for row in summary_rows],
        dtype=float,
    )
    rotation_combined = np.asarray(
        [row["rotation_equivalent_span_area_error_percent"] for row in summary_rows],
        dtype=float,
    )
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 8.2), sharex=True)
    panel_data = (
        (
            translation_combined,
            THESIS_COLORS[1],
            "Equivalent translational error",
            r"$e_r$ [%]",
        ),
        (
            rotation_combined,
            THESIS_COLORS[0],
            "Equivalent rotational error",
            r"$e_r$ [%]",
        ),
    )
    for axis, (values, color, title, ylabel) in zip(axes, panel_data):
        axis.plot(
            modes,
            values,
            color=color,
            marker="o",
            markersize=5.5,
            linewidth=2.1,
            linestyle="-",
        )
        value_range = float(np.max(values) - np.min(values))
        panel_margin = max(0.08 * value_range, 0.003)
        axis.set_ylim(
            float(np.min(values)) - panel_margin,
            float(np.max(values)) + panel_margin,
        )
        axis.set_title(title, fontweight="normal")
        axis.set_ylabel(ylabel)
        axis.grid(False)
        axis.yaxis.grid(
            True, which="major", color="0.82", linewidth=0.7, alpha=0.65
        )
        axis.set_axisbelow(True)
    axes[-1].set_xlabel(r"$N$")
    axes[-1].set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54])
    axes[-1].set_xlim(0, 56)
    fig.tight_layout()
    fig.savefig(
        args.output / "gravity_5g_translation_rotation_span_area_error.png",
        dpi=450,
        bbox_inches="tight",
    )
    fig.savefig(
        args.output / "gravity_5g_translation_rotation_span_area_error.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    best = min(
        summary_rows, key=lambda row: row["combined_span_area_error_percent"]
    )
    print(f"Integrated {len(RIGHT_SEMISPAN_NODES)} nodes over b = {length:.3f} in.")
    print(
        "Minimum span-normalized combined area error: "
        f"{best['combined_span_area_error_percent']:.6g}% at "
        f"N={best['elastic_modes']} elastic modes"
    )
    translation_minimum = int(modes[int(np.argmin(translation_combined))])
    rotation_minimum = int(modes[int(np.argmin(rotation_combined))])
    print(
        f"Minimum translation area error: {np.min(translation_combined):.6g}% "
        f"at N={translation_minimum}"
    )
    print(
        f"Minimum rotation-equivalent area error: {np.min(rotation_combined):.6g}% "
        f"at N={rotation_minimum}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
