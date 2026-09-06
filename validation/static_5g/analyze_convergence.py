#!/usr/bin/env python3
"""Compare MBDyn modal truncations against Nastran 5g tip displacement."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from cycler import cycler

from plot_modal_truncation_convergence import generate_plot as generate_modal_plot


TIP_NODES = (990020, 991020)
RIGHT_SEMISPAN_NODES = (990001, *range(991002, 991021))
THESIS_COLORS = ["#8B0000", "#00008B", "#66B2FF", "#006400", "#CC5500"]
TITLE_SIZE = 15
LABEL_SIZE = 15
LEGEND_SIZE = 15
TICK_SIZE = 15
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.titlesize": TITLE_SIZE,
    "axes.labelsize": LABEL_SIZE,
    "legend.fontsize": LEGEND_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "axes.prop_cycle": cycler(color=THESIS_COLORS),
})
FLOAT = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[EDed][+-]?\d+)?"
F06_ROW = re.compile(
    rf"^\s*(\d+)\s+G\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})"
    rf"\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$"
)


def nastran_float(value: str) -> float:
    value = value.replace("D", "E").replace("d", "E")
    if "E" not in value.upper():
        value = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", value)
    return float(value)


def read_nastran_tips(path: Path) -> dict[int, np.ndarray]:
    values: dict[int, np.ndarray] = {}
    in_table = False
    fatal: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "FATAL MESSAGE" in line:
                fatal.append(line.strip())
            if "D I S P L A C E M E N T   V E C T O R" in line:
                in_table = True
                continue
            if not in_table:
                continue
            match = F06_ROW.match(line)
            if match and int(match.group(1)) in TIP_NODES:
                values[int(match.group(1))] = np.asarray(
                    [nastran_float(match.group(i)) for i in range(2, 5)]
                )
    if fatal:
        raise RuntimeError(f"Nastran fatal message: {fatal[0]}")
    missing = set(TIP_NODES) - values.keys()
    if missing:
        raise RuntimeError(f"Tip displacement missing from F06: {sorted(missing)}")
    return values


def read_nastran_displacements(
    path: Path, nodes: tuple[int, ...]
) -> dict[int, np.ndarray]:
    """Read requested translations and rotations from the last F06 table."""
    requested = set(nodes)
    values: dict[int, np.ndarray] = {}
    in_table = False
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "D I S P L A C E M E N T   V E C T O R" in line:
                in_table = True
                values = {}
                continue
            if not in_table:
                continue
            match = F06_ROW.match(line)
            if match and int(match.group(1)) in requested:
                values[int(match.group(1))] = np.asarray(
                    [nastran_float(match.group(i)) for i in range(2, 8)]
                )
    missing = requested - values.keys()
    if missing:
        raise RuntimeError(
            f"Semispan displacement missing from F06: {sorted(missing)}"
        )
    return values


def read_mbdyn_tips(path: Path, tail_fraction: float) -> tuple[dict[int, np.ndarray], float]:
    tips: dict[int, np.ndarray] = {}
    worst_range = 0.0
    with netCDF4.Dataset(path) as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        start = max(1, int((1.0 - tail_fraction) * len(time)))
        for node in TIP_NODES:
            position = np.asarray(dataset.variables[f"node.struct.{node}.X"][:], dtype=float)
            displacement = position[:, :3] - position[0, :3]
            tips[node] = np.mean(displacement[start:], axis=0)
            worst_range = max(worst_range, float(np.max(np.ptp(displacement[start:], axis=0))))
    return tips, worst_range


def read_mbdyn_semispan(
    path: Path, nodes: tuple[int, ...], tail_fraction: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return span coordinate, translation and rotation at semispan nodes."""
    coordinates = []
    displacements = []
    rotations = []
    with netCDF4.Dataset(path) as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        start = max(1, int((1.0 - tail_fraction) * len(time)))
        for node in nodes:
            key = f"node.struct.{node}.X"
            if key not in dataset.variables:
                raise RuntimeError(f"MBDyn node {node} missing from {path.name}")
            position = np.asarray(dataset.variables[key][:], dtype=float)
            rotation_key = f"node.struct.{node}.Phi"
            if rotation_key not in dataset.variables:
                raise RuntimeError(f"MBDyn rotation for node {node} missing from {path.name}")
            rotation = np.asarray(dataset.variables[rotation_key][:], dtype=float)
            coordinates.append(position[0, 1])
            displacements.append(np.mean(position[start:, :3] - position[0, :3], axis=0))
            rotations.append(np.mean(rotation[start:, :3] - rotation[0, :3], axis=0))
    coordinates = np.asarray(coordinates, dtype=float)
    coordinates -= coordinates[0]
    return (
        coordinates,
        np.asarray(displacements, dtype=float),
        np.asarray(rotations, dtype=float),
    )


def sweep_colors(number: int) -> list[str]:
    """Preserve the thesis colors first, then add distinguishable solid colors."""
    additions = [
        "#6A0DAD", "#008B8B", "#808000", "#C71585", "#4B4B4B",
        "#8B4513", "#2E8B57", "#4682B4", "#B8860B", "#7B68EE",
    ]
    colors = THESIS_COLORS + additions
    if number > len(colors):
        needed = number - len(colors)
        colors.extend(
            matplotlib.colors.to_hex(color)
            for color in plt.get_cmap("turbo")(np.linspace(0.05, 0.95, needed))
        )
    return colors[:number]


def stable_recommendation(
    rows: list[dict], threshold: float, error_key: str
) -> tuple[int | None, bool]:
    errors = [
        max(float(row[error_key]), float(row["tail_range_in"]))
        for row in rows
    ]
    for index, row in enumerate(rows):
        window = errors[index : index + 3]
        if all(error <= threshold for error in window):
            return int(row["elastic_modes"]), len(window) == 3
    return None, False


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "f06", nargs="?", type=Path,
        default=root / "nastran/MAIN/sol101_gravity_5g.f06"
    )
    parser.add_argument("--results", type=Path, default=root / "mbdyn/results")
    parser.add_argument("--output", type=Path, default=root / "plots")
    parser.add_argument("--absolute-threshold", type=float, default=0.01)
    parser.add_argument("--relative-threshold", type=float, default=0.01)
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    args = parser.parse_args()

    nastran = read_nastran_tips(args.f06)
    nastran_semispan = read_nastran_displacements(args.f06, RIGHT_SEMISPAN_NODES)
    reference = max(float(np.linalg.norm(value)) for value in nastran.values())
    threshold = max(args.absolute_threshold, args.relative_threshold * reference)
    rows: list[dict] = []
    semispan_cases: list[dict] = []
    pattern = re.compile(r"gravity_5g_(\d+)_elastic_modes\.nc$")
    for path in sorted(args.results.glob("gravity_5g_*_elastic_modes.nc")):
        match = pattern.search(path.name)
        if not match:
            continue
        count = int(match.group(1))
        mbdyn, tail_range = read_mbdyn_tips(path, args.tail_fraction)
        span_b, span_displacement, span_rotation = read_mbdyn_semispan(
            path, RIGHT_SEMISPAN_NODES, args.tail_fraction
        )
        semispan_cases.append(
            {
                "elastic_modes": count,
                "b_in": span_b,
                "displacement_in": span_displacement,
                "rotation_rad": span_rotation,
            }
        )
        errors = {node: mbdyn[node] - nastran[node] for node in TIP_NODES}
        rows.append(
            {
                "elastic_modes": count,
                "left_Ux_in": mbdyn[990020][0],
                "left_Uy_in": mbdyn[990020][1],
                "left_Uz_in": mbdyn[990020][2],
                "right_Ux_in": mbdyn[991020][0],
                "right_Uy_in": mbdyn[991020][1],
                "right_Uz_in": mbdyn[991020][2],
                "left_vector_error_in": np.linalg.norm(errors[990020]),
                "right_vector_error_in": np.linalg.norm(errors[991020]),
                "max_tip_vector_error_in": max(
                    np.linalg.norm(errors[990020]), np.linalg.norm(errors[991020])
                ),
                "max_tip_Uz_error_in": max(abs(errors[990020][2]), abs(errors[991020][2])),
                "tail_range_in": tail_range,
            }
        )
    rows.sort(key=lambda row: int(row["elastic_modes"]))
    semispan_cases.sort(key=lambda case: int(case["elastic_modes"]))
    if not rows:
        raise RuntimeError(f"No MBDyn NetCDF convergence results found in {args.results}")

    richest = rows[-1]
    for row in rows:
        internal_vector_errors = []
        internal_uz_errors = []
        for side in ("left", "right"):
            current = np.asarray(
                [row[f"{side}_{component}_in"] for component in ("Ux", "Uy", "Uz")],
                dtype=float,
            )
            reference_tip = np.asarray(
                [
                    richest[f"{side}_{component}_in"]
                    for component in ("Ux", "Uy", "Uz")
                ],
                dtype=float,
            )
            internal_vector_errors.append(float(np.linalg.norm(current - reference_tip)))
            internal_uz_errors.append(float(abs(current[2] - reference_tip[2])))
        row["internal_vector_error_vs_richest_in"] = max(internal_vector_errors)
        row["internal_Uz_error_vs_richest_in"] = max(internal_uz_errors)

    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "gravity_5g_convergence.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    nastran_span_array = np.asarray(
        [nastran_semispan[node] for node in RIGHT_SEMISPAN_NODES], dtype=float
    )
    span_table_fields = [
        "elastic_modes", "node", "b_in",
        "Ux_MBDyn_in", "Uy_MBDyn_in", "Uz_MBDyn_in",
        "Ux_Nastran_in", "Uy_Nastran_in", "Uz_Nastran_in",
        "abs_Ux_error_in", "abs_Uy_error_in", "abs_Uz_error_in",
        "Rx_MBDyn_rad", "Ry_MBDyn_rad", "Rz_MBDyn_rad",
        "Rx_Nastran_rad", "Ry_Nastran_rad", "Rz_Nastran_rad",
        "abs_Rx_error_rad", "abs_Ry_error_rad", "abs_Rz_error_rad",
    ]
    with (args.output / "gravity_5g_semispan_absolute_errors.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=span_table_fields)
        writer.writeheader()
        for case in semispan_cases:
            for index, node in enumerate(RIGHT_SEMISPAN_NODES):
                mbdyn_value = case["displacement_in"][index]
                mbdyn_rotation = case["rotation_rad"][index]
                nastran_value = nastran_span_array[index, :3]
                nastran_rotation = nastran_span_array[index, 3:]
                absolute_error = np.abs(mbdyn_value - nastran_value)
                absolute_rotation_error = np.abs(mbdyn_rotation - nastran_rotation)
                writer.writerow(
                    {
                        "elastic_modes": case["elastic_modes"],
                        "node": node,
                        "b_in": case["b_in"][index],
                        "Ux_MBDyn_in": mbdyn_value[0],
                        "Uy_MBDyn_in": mbdyn_value[1],
                        "Uz_MBDyn_in": mbdyn_value[2],
                        "Ux_Nastran_in": nastran_value[0],
                        "Uy_Nastran_in": nastran_value[1],
                        "Uz_Nastran_in": nastran_value[2],
                        "abs_Ux_error_in": absolute_error[0],
                        "abs_Uy_error_in": absolute_error[1],
                        "abs_Uz_error_in": absolute_error[2],
                        "Rx_MBDyn_rad": mbdyn_rotation[0],
                        "Ry_MBDyn_rad": mbdyn_rotation[1],
                        "Rz_MBDyn_rad": mbdyn_rotation[2],
                        "Rx_Nastran_rad": nastran_rotation[0],
                        "Ry_Nastran_rad": nastran_rotation[1],
                        "Rz_Nastran_rad": nastran_rotation[2],
                        "abs_Rx_error_rad": absolute_rotation_error[0],
                        "abs_Ry_error_rad": absolute_rotation_error[1],
                        "abs_Rz_error_rad": absolute_rotation_error[2],
                    }
                )

    counts = np.asarray([row["elastic_modes"] for row in rows])
    left_error = np.asarray([row["left_vector_error_in"] for row in rows])
    right_error = np.asarray([row["right_vector_error_in"] for row in rows])
    max_error = np.asarray([row["max_tip_vector_error_in"] for row in rows])
    uz_error = np.asarray([row["max_tip_Uz_error_in"] for row in rows])

    fig, axis = plt.subplots(figsize=(10, 6))
    axis.semilogy(
        counts, left_error, marker="o", linewidth=1.8, linestyle="-",
        label="Left tip vector error"
    )
    axis.semilogy(
        counts, right_error, marker="s", linewidth=1.8, linestyle="-",
        label="Right tip vector error"
    )
    axis.semilogy(
        counts, uz_error, marker="^", linewidth=1.8, linestyle="-",
        label="Maximum tip Uz error"
    )
    axis.axhline(
        threshold, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8,
        label=f"Threshold = {threshold:.4g} in"
    )
    axis.set_xlabel("Number of elastic modes")
    axis.set_ylabel("Absolute tip error [in]")
    axis.set_title("External validation against Nastran SOL 101 under 5g gravity")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(args.output / "gravity_5g_nastran_validation_error.png", dpi=180)
    plt.close(fig)

    # Publication-ready single-curve validation plot using Nastran as reference.
    fig, axis = plt.subplots(figsize=(9.2, 6.0))
    axis.plot(
        counts,
        max_error,
        color=THESIS_COLORS[1],
        marker="o",
        markersize=5.5,
        linewidth=2.1,
        linestyle="-",
        zorder=3,
    )
    axis.axhline(
        threshold,
        color="0.25",
        linewidth=1.2,
        linestyle="-",
        alpha=0.8,
        zorder=1,
    )
    axis.set_title("MBDyn–Nastran static convergence", fontweight="normal")
    axis.set_xlabel("Number of elastic modes, N")
    axis.set_ylabel("Maximum MBDyn–Nastran tip error [in]")
    axis.set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54])
    axis.set_xlim(0, 56)
    axis.set_ylim(1.1, 1.3)
    axis.grid(False)
    axis.yaxis.grid(True, which="major", color="0.82", linewidth=0.7, alpha=0.65)
    axis.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(
        args.output / "gravity_5g_nastran_static_convergence.png",
        dpi=450,
        bbox_inches="tight",
    )
    fig.savefig(
        args.output / "gravity_5g_nastran_static_convergence.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharex=True)
    for axis, node, label in zip(axes, TIP_NODES, ("Left tip", "Right tip")):
        for component, component_label in enumerate(("Ux", "Uy", "Uz")):
            values = np.asarray([
                row[f"{'left' if node == 990020 else 'right'}_{component_label}_in"]
                for row in rows
            ])
            axis.plot(
                counts, values, marker="o", linewidth=1.8, linestyle="-",
                label=f"MBDyn {component_label}"
            )
            axis.axhline(
                nastran[node][component], color="0.25", linestyle="-",
                linewidth=1.2, alpha=0.8,
                label=f"Nastran {component_label}"
            )
        axis.set_title(label)
        axis.set_xlabel("Number of elastic modes")
        axis.set_ylabel("Tip displacement [in]")
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=2, fontsize=LEGEND_SIZE)
    fig.suptitle("5g gravity tip-displacement convergence", fontsize=TITLE_SIZE)
    fig.tight_layout()
    fig.savefig(args.output / "gravity_5g_tip_displacements.png", dpi=180)
    plt.close(fig)

    displayed_counts = (1, 8, 25, 54)
    displayed_cases = [
        case for case in semispan_cases
        if int(case["elastic_modes"]) in displayed_counts
    ]
    missing_displayed = set(displayed_counts) - {
        int(case["elastic_modes"]) for case in displayed_cases
    }
    if missing_displayed:
        raise RuntimeError(
            f"Semispan plot cases are missing: {sorted(missing_displayed)}"
        )
    colors = THESIS_COLORS[:len(displayed_cases)]
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "*", "p")
    displacement_titles = (
        "Displacement in x direction absolute error",
        "Displacement in y direction absolute error",
        "Displacement in z direction absolute error",
    )
    for component_index, (component_name, plot_title) in enumerate(
        zip(("Ux", "Uy", "Uz"), displacement_titles)
    ):
        fig, axis = plt.subplots(figsize=(11.5, 7.0))
        for case_index, case in enumerate(displayed_cases):
            error = np.abs(
                case["displacement_in"][:, component_index]
                - nastran_span_array[:, component_index]
            )
            axis.plot(
                case["b_in"], error,
                color=colors[case_index],
                marker=markers[case_index % len(markers)],
                markevery=2,
                markersize=4.5,
                linewidth=1.8,
                linestyle="-",
                label=(
                    "1 elastic mode" if int(case["elastic_modes"]) == 1
                    else f"{case['elastic_modes']} elastic modes"
                ),
            )
        axis.set_xlabel("b [in]", fontsize=LABEL_SIZE)
        direction = ("x", "y", "z")[component_index]
        axis.set_ylabel(
            r"$|u_{" + direction + r",\mathrm{NASTRAN}}"
            r"-u_{" + direction + r",\mathrm{MBDyn},N}|$ [in]",
            fontsize=LABEL_SIZE,
        )
        axis.set_title(
            plot_title, fontsize=TITLE_SIZE, fontweight="normal"
        )
        axis.tick_params(axis="both", labelsize=TICK_SIZE)
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=2, fontsize=LEGEND_SIZE, frameon=False)
        fig.tight_layout()
        fig.savefig(
            args.output / f"gravity_5g_semispan_absolute_error_{component_name}.png",
            dpi=450,
        )
        plt.close(fig)

    rotation_titles = (
        "Rotation about x axis absolute error",
        "Rotation about y axis absolute error",
        "Rotation about z axis absolute error",
    )
    for component_index, (component_name, plot_title) in enumerate(
        zip(("Rx", "Ry", "Rz"), rotation_titles)
    ):
        fig, axis = plt.subplots(figsize=(11.5, 7.0))
        for case_index, case in enumerate(displayed_cases):
            error = np.abs(
                case["rotation_rad"][:, component_index]
                - nastran_span_array[:, component_index + 3]
            )
            axis.plot(
                case["b_in"], error,
                color=colors[case_index],
                marker=markers[case_index % len(markers)],
                markevery=2,
                markersize=4.5,
                linewidth=1.8,
                linestyle="-",
                label=(
                    "1 elastic mode" if int(case["elastic_modes"]) == 1
                    else f"{case['elastic_modes']} elastic modes"
                ),
            )
        axis.set_xlabel("b [in]", fontsize=LABEL_SIZE)
        direction = ("x", "y", "z")[component_index]
        axis.set_ylabel(
            r"$|\theta_{" + direction + r",\mathrm{NASTRAN}}"
            r"-\theta_{" + direction + r",\mathrm{MBDyn},N}|$ [rad]",
            fontsize=LABEL_SIZE,
        )
        axis.set_title(
            plot_title, fontsize=TITLE_SIZE, fontweight="normal"
        )
        axis.tick_params(axis="both", labelsize=TICK_SIZE)
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=2, fontsize=LEGEND_SIZE, frameon=False)
        fig.tight_layout()
        fig.savefig(
            args.output / f"gravity_5g_semispan_absolute_error_{component_name}.png",
            dpi=450,
        )
        plt.close(fig)

    # Six-panel thesis figure. The enlarged typography compensates for the
    # reduction of each panel when the complete figure is fitted to a page.
    combined_font_size = 24
    combined_ylabel_size = 30
    combined_title_size = 30
    fig, axes = plt.subplots(3, 2, figsize=(17.0, 19.0), sharex=True)
    panel_specs = (
        ("displacement_in", 0, r"$u_x$ absolute error", "u", "x", "in"),
        ("displacement_in", 1, r"$u_y$ absolute error", "u", "y", "in"),
        ("displacement_in", 2, r"$u_z$ absolute error", "u", "z", "in"),
        ("rotation_rad", 0, r"$\theta_x$ absolute error", r"\theta", "x", "rad"),
        ("rotation_rad", 1, r"$\theta_y$ absolute error", r"\theta", "y", "rad"),
        ("rotation_rad", 2, r"$\theta_z$ absolute error", r"\theta", "z", "rad"),
    )
    legend_handles = None
    legend_labels = None
    for axis, (field, component, title, symbol, direction, unit) in zip(
        axes.flat, panel_specs
    ):
        nastran_component = component if field == "displacement_in" else component + 3
        for case_index, case in enumerate(displayed_cases):
            modes = int(case["elastic_modes"])
            error = np.abs(
                case[field][:, component] - nastran_span_array[:, nastran_component]
            )
            axis.plot(
                case["b_in"],
                error,
                color=colors[case_index],
                marker=markers[case_index],
                markevery=2,
                markersize=7.0,
                linewidth=2.1,
                linestyle="-",
                label="1 elastic mode" if modes == 1 else f"{modes} elastic modes",
            )
        axis.set_title(title, fontsize=combined_title_size, fontweight="normal")
        axis.set_xlabel("b [in]", fontsize=combined_font_size)
        axis.set_ylabel(
            "$|"
            + symbol
            + "_{"
            + direction
            + r",\mathrm{NAST}}-"
            + symbol
            + "_{"
            + direction
            + r",\mathrm{MBDyn},N}|$ ["
            + unit
            + "]",
            fontsize=combined_ylabel_size,
        )
        axis.tick_params(axis="both", labelsize=combined_font_size)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=False)
        axis.yaxis.get_offset_text().set_fontsize(combined_font_size)
        axis.grid(True, alpha=0.3)
        if legend_handles is None:
            legend_handles, legend_labels = axis.get_legend_handles_labels()
    fig.legend(
        legend_handles,
        legend_labels,
        loc="lower center",
        ncol=4,
        fontsize=combined_font_size,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="0.25",
        bbox_to_anchor=(0.5, 0.005),
    )
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 1.0), pad=1.3, w_pad=2.0, h_pad=1.5)
    fig.savefig(
        args.output / "gravity_5g_semispan_absolute_errors_all.png",
        dpi=600,
        bbox_inches="tight",
    )
    fig.savefig(
        args.output / "gravity_5g_semispan_absolute_errors_all.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    # Dimensionless cumulative error over all six kinematic components and
    # every right-semispan node. Translation and rotation fields receive equal
    # weight after normalization by their respective maximum Nastran norm.
    node_count = len(RIGHT_SEMISPAN_NODES)
    translation_reference = float(
        np.max(np.linalg.norm(nastran_span_array[:, :3], axis=1))
    )
    rotation_reference = float(
        np.max(np.linalg.norm(nastran_span_array[:, 3:], axis=1))
    )
    if translation_reference <= 0.0 or rotation_reference <= 0.0:
        raise RuntimeError("Nastran semispan normalization reference is zero")
    cumulative_rows = []
    for case in semispan_cases:
        translation_difference = (
            case["displacement_in"] - nastran_span_array[:, :3]
        )
        rotation_difference = case["rotation_rad"] - nastran_span_array[:, 3:]
        translation_term = float(
            np.mean(np.sum(translation_difference**2, axis=1))
            / translation_reference**2
        )
        rotation_term = float(
            np.mean(np.sum(rotation_difference**2, axis=1))
            / rotation_reference**2
        )
        normalized_error_percent = 100.0 * np.sqrt(
            0.5 * (translation_term + rotation_term)
        )
        cumulative_rows.append(
            {
                "elastic_modes": int(case["elastic_modes"]),
                "spanwise_nodes": node_count,
                "normalized_cumulative_error_percent": normalized_error_percent,
                "translation_reference_in": translation_reference,
                "rotation_reference_rad": rotation_reference,
            }
        )

    cumulative_csv = args.output / "gravity_5g_normalized_cumulative_spanwise_error.csv"
    with cumulative_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=cumulative_rows[0].keys())
        writer.writeheader()
        writer.writerows(cumulative_rows)

    cumulative_modes = np.asarray(
        [row["elastic_modes"] for row in cumulative_rows], dtype=int
    )
    cumulative_error = np.asarray(
        [row["normalized_cumulative_error_percent"] for row in cumulative_rows],
        dtype=float,
    )
    fig, axis = plt.subplots(figsize=(9.2, 6.0))
    axis.plot(
        cumulative_modes,
        cumulative_error,
        color=THESIS_COLORS[0],
        marker="o",
        markersize=5.5,
        linewidth=2.1,
        linestyle="-",
        label=f"{node_count} spanwise nodes",
    )
    axis.set_title("Normalized cumulative spanwise error", fontweight="normal")
    axis.set_xlabel("Number of elastic modes, N")
    axis.set_ylabel("Normalized cumulative error [%]")
    axis.set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54])
    axis.set_xlim(0, 56)
    axis.set_ylim(bottom=0.0)
    axis.grid(False)
    axis.yaxis.grid(True, which="major", color="0.82", linewidth=0.7, alpha=0.65)
    axis.legend(frameon=False)
    axis.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(
        args.output / "gravity_5g_normalized_cumulative_spanwise_error.png",
        dpi=450,
        bbox_inches="tight",
    )
    fig.savefig(
        args.output / "gravity_5g_normalized_cumulative_spanwise_error.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    validation_recommendation, validation_confirmed = stable_recommendation(
        rows, threshold, "max_tip_vector_error_in"
    )
    internal_recommendation, internal_confirmed = stable_recommendation(
        rows, args.absolute_threshold, "internal_vector_error_vs_richest_in"
    )
    print(f"External Nastran threshold: {threshold:.6g} in")
    print(f"Best external tip error: {np.min(max_error):.6g} in")
    if validation_recommendation is None:
        print("External validation: no tested modal basis satisfies the threshold.")
    else:
        status = "confirmed" if validation_confirmed else "provisional"
        print(
            f"External validation minimum: {validation_recommendation} elastic modes "
            f"({status})."
        )
    if internal_recommendation is None:
        print("Internal modal convergence: no tested basis satisfies the threshold.")
    else:
        status = "confirmed" if internal_confirmed else "provisional"
        print(
            f"Internal modal convergence minimum: {internal_recommendation} elastic modes "
            f"({status})."
        )
    generate_modal_plot(
        args.output / "gravity_5g_convergence.csv",
        args.output / "gravity_5g_modal_truncation_error.png",
        args.output / "gravity_5g_modal_truncation_error.pdf",
        threshold=args.absolute_threshold,
        print_summary=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
