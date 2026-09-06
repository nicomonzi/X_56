#!/usr/bin/env python3
"""Plot the NASTRAN X-56 structural geometry and its fictitious beam."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MODEL_DIR = SCRIPT_DIR.parents[1] / "NASTRAN" / "NASTRAN40"

THESIS_COLORS = ["#8B0000", "#00008B", "#66B2FF", "#006400", "#CC5500"]
AIRCRAFT_COLOR = "#909090"
AIRCRAFT_EDGE = "#404040"
BEAM_COLOR = "#00B050"
BEAM_NODE_COLOR = "#005A24"
TEXT_SIZE = 17
TICK_SIZE = 15
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.titlesize": TEXT_SIZE,
    "axes.labelsize": TEXT_SIZE,
    "legend.fontsize": TEXT_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "axes.prop_cycle": cycler(color=THESIS_COLORS),
})


def nastran_float(value: str) -> float:
    """Read conventional and compact NASTRAN floating-point fields."""
    value = value.strip().replace("D", "E")
    if not value:
        return 0.0
    if "E" not in value.upper():
        for index in range(len(value) - 1, 0, -1):
            if value[index] in "+-" and value[index - 1].isdigit():
                value = value[:index] + "E" + value[index:]
                break
    return float(value)


def fields(line: str) -> list[str]:
    """Split a small-field NASTRAN record into its eight-character fields."""
    if "\t" in line:
        split_fields = line.split()
        return (split_fields + [""] * 10)[:10]
    padded = line.rstrip("\n").ljust(80)
    return [padded[index : index + 8].strip() for index in range(0, 80, 8)]


def read_grids(paths: list[Path]) -> dict[int, np.ndarray]:
    nodes: dict[int, np.ndarray] = {}
    for path in paths:
        with path.open(errors="replace") as stream:
            for line in stream:
                row = fields(line)
                if row[0] != "GRID":
                    continue
                node_id = int(row[1])
                coordinate_system = int(row[2] or 0)
                if coordinate_system != 0:
                    raise ValueError(
                        f"GRID {node_id} in {path} uses unsupported CP={coordinate_system}"
                    )
                nodes[node_id] = np.array(
                    [nastran_float(row[3]), nastran_float(row[4]), nastran_float(row[5])]
                )
    return nodes


def read_shells(paths: list[Path], nodes: dict[int, np.ndarray]) -> list[np.ndarray]:
    shell_sizes = {"CTRIA3": 3, "CTRIAR": 3, "CQUAD4": 4, "CQUADR": 4}
    shells: list[np.ndarray] = []
    for path in paths:
        with path.open(errors="replace") as stream:
            for line in stream:
                row = fields(line)
                node_count = shell_sizes.get(row[0])
                if node_count is None:
                    continue
                try:
                    node_ids = [int(value) for value in row[3 : 3 + node_count]]
                except ValueError:
                    continue
                if all(node_id in nodes for node_id in node_ids):
                    shells.append(np.vstack([nodes[node_id] for node_id in node_ids]))
    return shells


def read_beam(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes = read_grids([path])
    centre_id = 990001
    left_ids = [centre_id] + list(range(990002, 990025))
    right_ids = [centre_id] + list(range(991002, 991025))
    return (
        np.vstack([nodes[node_id] for node_id in left_ids]),
        np.vstack([nodes[node_id] for node_id in right_ids]),
        np.vstack([nodes[node_id] for node_id in left_ids[:0:-1] + right_ids]),
    )


def unique_edges(shells: list[np.ndarray]) -> np.ndarray:
    edges = []
    for shell in shells:
        for index in range(len(shell)):
            edges.append(np.vstack((shell[index], shell[(index + 1) % len(shell)])))
    return np.asarray(edges)


def padded_limits(values: np.ndarray, fraction: float = 0.04) -> tuple[float, float]:
    low, high = float(values.min()), float(values.max())
    margin = max((high - low) * fraction, 1.0)
    return low - margin, high + margin


def style_2d_axis(ax: plt.Axes, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=TEXT_SIZE, weight="normal")
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.grid(True, color="#D9E0E5", linewidth=0.55, alpha=0.75)
    ax.set_axisbelow(True)
    ax.set_facecolor("#FAFCFD")


def draw_orthogonal(
    ax: plt.Axes,
    shells: list[np.ndarray],
    edges: np.ndarray,
    beam: np.ndarray,
    axes: tuple[int, int],
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    first, second = axes
    projected_shells = [shell[:, [first, second]] for shell in shells]
    projected_edges = edges[:, :, [first, second]]
    ax.add_collection(
        PolyCollection(
            projected_shells,
            facecolors=AIRCRAFT_COLOR,
            edgecolors="none",
            alpha=0.55,
            rasterized=True,
        )
    )
    ax.add_collection(
        LineCollection(
            projected_edges,
            colors=AIRCRAFT_EDGE,
            linewidths=0.28,
            alpha=0.38,
            rasterized=True,
        )
    )
    ax.plot(
        beam[:, first],
        beam[:, second],
        color=BEAM_COLOR,
        linewidth=3.8,
        linestyle="-",
        marker="o",
        markersize=4.0,
        markerfacecolor=BEAM_NODE_COLOR,
        markeredgewidth=0,
        zorder=10,
    )
    all_points = np.vstack(shells)
    ax.set_xlim(*padded_limits(all_points[:, first]))
    ax.set_ylim(*padded_limits(all_points[:, second]))
    style_2d_axis(ax, xlabel, ylabel, title)


def save_figure(fig: plt.Figure, stem: Path, svg: bool = True) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    if svg:
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")


def create_plots(model_dir: Path, output_dir: Path) -> None:
    bulk_paths = sorted((model_dir / "BULK").glob("*.dat"))
    beam_path = model_dir / "MAIN" / "rbe3s.bdf"
    if not bulk_paths or not beam_path.is_file():
        raise FileNotFoundError(f"NASTRAN40 model not found below {model_dir}")

    nodes = read_grids(bulk_paths)
    shells = read_shells(bulk_paths, nodes)
    left_beam, right_beam, full_beam = read_beam(beam_path)
    edges = unique_edges(shells)
    output_dir.mkdir(parents=True, exist_ok=True)

    legend_handles = [
        Line2D([0], [0], color=AIRCRAFT_COLOR, linewidth=7, alpha=0.75,
               label="X-56 Geometry"),
        Line2D([0], [0], color=BEAM_COLOR, linewidth=3.8, marker="o",
               markersize=5, label="Equivalent Beam"),
    ]

    views = [
        ((1, 0), "Top View", "Y [in]", "X [in]", "top_view"),
        ((1, 2), "Front View", "Y [in]", "Z [in]", "front_view"),
        ((0, 2), "Side View", "X [in]", "Z [in]", "side_view"),
    ]

    fig, combined_axes = plt.subplots(
        3,
        1,
        figsize=(13.2, 13.0),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3.4, 1.0, 1.5]},
    )
    for ax, (projection, title, xlabel, ylabel, filename) in zip(combined_axes, views):
        draw_orthogonal(ax, shells, edges, full_beam, projection, title, xlabel, ylabel)
        if filename == "top_view":
            ax.invert_yaxis()
    combined_axes[1].sharex(combined_axes[0])
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.035),
        ncol=2,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        edgecolor="0.25",
        fontsize=TEXT_SIZE,
    )
    save_figure(fig, output_dir / "x56_beam_orthogonal_views", svg=False)
    plt.close(fig)

    for projection, title, xlabel, ylabel, filename in views:
        fig, ax = plt.subplots(figsize=(11.5, 7.2), constrained_layout=True)
        draw_orthogonal(ax, shells, edges, full_beam, projection, title, xlabel, ylabel)
        if filename == "top_view":
            ax.invert_yaxis()
        legend_options = {
            "handles": legend_handles,
            "frameon": True,
            "fancybox": False,
            "framealpha": 1.0,
            "edgecolor": "0.25",
        }
        if filename == "front_view":
            ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, -0.52),
                ncol=2,
                **legend_options,
            )
        else:
            ax.legend(loc="best", **legend_options)
        save_figure(fig, output_dir / f"x56_beam_{filename}")
        plt.close(fig)

    fig = plt.figure(figsize=(13.0, 8.2))
    ax = fig.add_axes((0.035, 0.075, 0.88, 0.84), projection="3d")
    ax.add_collection3d(
        Poly3DCollection(
            shells,
            facecolors=AIRCRAFT_COLOR,
            edgecolors="none",
            alpha=0.32,
            rasterized=True,
        )
    )
    ax.add_collection3d(
        Line3DCollection(
            edges,
            colors=AIRCRAFT_EDGE,
            linewidths=0.25,
            alpha=0.28,
            rasterized=True,
        )
    )
    for half in (left_beam, right_beam):
        ax.plot(
            half[:, 0],
            half[:, 1],
            half[:, 2],
            color=BEAM_COLOR,
            linewidth=4.2,
            linestyle="-",
            marker="o",
            markersize=4.2,
            markerfacecolor=BEAM_NODE_COLOR,
            markeredgewidth=0,
            zorder=20,
        )

    all_points = np.vstack(shells)
    ax.set_xlim(*padded_limits(all_points[:, 0]))
    ax.set_ylim(*padded_limits(all_points[:, 1]))
    ax.set_zlim(*padded_limits(all_points[:, 2]))
    spans = np.ptp(all_points, axis=0)
    ax.set_box_aspect(spans, zoom=0.98)
    ax.view_init(elev=24, azim=-58)
    ax.set_xlabel("X [in]", labelpad=10)
    ax.set_ylabel("Y [in]", labelpad=10)
    ax.set_zlabel("Z [in]", labelpad=8)
    ax.set_title(
        "3D View",
        fontsize=TEXT_SIZE,
        weight="normal",
        pad=8,
    )
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.legend(
        handles=legend_handles, loc="upper left",
        bbox_to_anchor=(0.0, 0.76), frameon=True,
        fancybox=False, framealpha=1.0, edgecolor="0.25"
    )
    ax.grid(True, linewidth=0.45, alpha=0.55)
    save_figure(fig, output_dir / "x56_beam_3d_view")
    plt.close(fig)

    print(f"Read {len(nodes)} structural nodes and {len(shells)} shell elements.")
    print(f"Read {len(full_beam)} fictitious-beam nodes.")
    print(f"Plots written to: {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=DEFAULT_MODEL_DIR,
        help="Path to NASTRAN40 (default: repository NASTRAN/NASTRAN40)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory receiving the figures (default: this script's directory)",
    )
    arguments = parser.parse_args()
    create_plots(arguments.model_dir.resolve(), arguments.output_dir.resolve())


if __name__ == "__main__":
    main()
