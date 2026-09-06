#!/usr/bin/env python3
"""Plot cumulative MBDyn--Nastran errors from root to each semispan station."""

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


THESIS_COLORS = [
    "#8B0000",
    "#00008B",
    "#66B2FF",
    "#006400",
    "#CC5500",
]
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

DISPLAYED_MODES = (1, 8, 25, 54)


def cumulative_l2(values: np.ndarray) -> np.ndarray:
    """Euclidean accumulation from the root through every current station."""
    return np.sqrt(np.cumsum(np.square(values), axis=0))


def modal_result_path(results: Path, modes: int) -> Path:
    return results / f"gravity_5g_{modes:02d}_elastic_modes.nc"


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
    parser.add_argument(
        "--results", type=Path, default=STUDY / "mbdyn/results"
    )
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
        path = modal_result_path(args.results, modes)
        if not path.is_file():
            raise FileNotFoundError(f"Missing MBDyn result: {path}")
        b, displacement, rotation = read_mbdyn_semispan(
            path, RIGHT_SEMISPAN_NODES, args.tail_fraction
        )
        difference = np.column_stack(
            (displacement - nastran_field[:, :3], rotation - nastran_field[:, 3:])
        )
        cumulative_components = cumulative_l2(difference)

        translation_error_energy = np.cumsum(
            np.sum(np.square(difference[:, :3]), axis=1)
        )
        rotation_error_energy = np.cumsum(
            np.sum(np.square(difference[:, 3:]), axis=1)
        )
        translation_reference_energy = np.cumsum(
            np.sum(np.square(nastran_field[:, :3]), axis=1)
        )
        rotation_reference_energy = np.cumsum(
            np.sum(np.square(nastran_field[:, 3:]), axis=1)
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            translation_relative = np.divide(
                translation_error_energy,
                translation_reference_energy,
                out=np.zeros_like(translation_error_energy),
                where=translation_reference_energy > 0.0,
            )
            rotation_relative = np.divide(
                rotation_error_energy,
                rotation_reference_energy,
                out=np.zeros_like(rotation_error_energy),
                where=rotation_reference_energy > 0.0,
            )
        global_percent = 100.0 * np.sqrt(
            0.5 * (translation_relative + rotation_relative)
        )
        cases.append(
            {
                "modes": modes,
                "b": b,
                "local_difference": difference,
                "cumulative_components": cumulative_components,
                "global_percent": global_percent,
            }
        )

    fields = ["elastic_modes", "node", "b_in"]
    for name, unit in zip(("Ux", "Uy", "Uz", "Rx", "Ry", "Rz"), ("in", "in", "in", "rad", "rad", "rad")):
        fields.extend((f"local_{name}_error_{unit}", f"cumulative_{name}_error_{unit}"))
    fields.append("cumulative_global_error_percent")
    with (args.output / "gravity_5g_cumulative_spanwise_errors.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for case in cases:
            for index, node in enumerate(RIGHT_SEMISPAN_NODES):
                row = {
                    "elastic_modes": case["modes"],
                    "node": node,
                    "b_in": case["b"][index],
                    "cumulative_global_error_percent": case["global_percent"][index],
                }
                for component, (name, unit) in enumerate(
                    zip(("Ux", "Uy", "Uz", "Rx", "Ry", "Rz"), ("in", "in", "in", "rad", "rad", "rad"))
                ):
                    row[f"local_{name}_error_{unit}"] = abs(
                        case["local_difference"][index, component]
                    )
                    row[f"cumulative_{name}_error_{unit}"] = case[
                        "cumulative_components"
                    ][index, component]
                writer.writerow(row)

    titles = (
        "Displacement in x direction cumulative error",
        "Displacement in y direction cumulative error",
        "Displacement in z direction cumulative error",
        "Rotation about x axis cumulative error",
        "Rotation about y axis cumulative error",
        "Rotation about z axis cumulative error",
    )
    names = ("Ux", "Uy", "Uz", "Rx", "Ry", "Rz")
    symbols = ("u_x", "u_y", "u_z", r"\theta_x", r"\theta_y", r"\theta_z")
    units = ("in", "in", "in", "rad", "rad", "rad")
    markers = ("o", "s", "^", "D")
    for component, (title, name, symbol, unit) in enumerate(
        zip(titles, names, symbols, units)
    ):
        fig, axis = plt.subplots(figsize=(11.5, 7.0))
        for index, case in enumerate(cases):
            modes = int(case["modes"])
            axis.plot(
                case["b"],
                case["cumulative_components"][:, component],
                color=THESIS_COLORS[index],
                marker=markers[index],
                markevery=2,
                markersize=4.5,
                linewidth=1.8,
                linestyle="-",
                label="1 elastic mode" if modes == 1 else f"{modes} elastic modes",
            )
        axis.set_title(title, fontweight="normal")
        axis.set_xlabel("b [in]")
        axis.set_ylabel(rf"$E^{{\mathrm{{cum}}}}_{{{symbol},N}}$ [{unit}]")
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=2, frameon=False)
        save_figure(fig, args.output, f"gravity_5g_cumulative_spanwise_error_{name}")

    fig, axis = plt.subplots(figsize=(11.5, 7.0))
    for index, case in enumerate(cases):
        modes = int(case["modes"])
        axis.plot(
            case["b"],
            case["global_percent"],
            color=THESIS_COLORS[index],
            marker=markers[index],
            markevery=2,
            markersize=4.5,
            linewidth=1.8,
            linestyle="-",
            label="1 elastic mode" if modes == 1 else f"{modes} elastic modes",
        )
    axis.set_title("Global cumulative spanwise error", fontweight="normal")
    axis.set_xlabel("b [in]")
    axis.set_ylabel("Normalized global cumulative error [%]")
    axis.grid(True, alpha=0.3)
    axis.legend(ncol=2, frameon=False)
    save_figure(fig, args.output, "gravity_5g_global_cumulative_spanwise_error")

    print(f"Used {len(RIGHT_SEMISPAN_NODES)} right-semispan nodes (winglet excluded).")
    print(f"Created {len(names) + 1} cumulative plots in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
