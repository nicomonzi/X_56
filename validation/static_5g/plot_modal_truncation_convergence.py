#!/usr/bin/env python3
"""Create the publication-ready MBDyn modal-truncation convergence plot."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler


THESIS_COLORS = [
    "#8B0000",  # dark red
    "#00008B",  # dark blue
    "#66B2FF",  # light blue
    "#006400",  # dark green
    "#CC5500",  # dark orange
]
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


def normalized(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def find_column(fieldnames: list[str], candidates: tuple[str, ...]) -> str | None:
    lookup = {normalized(name): name for name in fieldnames}
    for candidate in candidates:
        if normalized(candidate) in lookup:
            return lookup[normalized(candidate)]
    return None


def identify_columns(fieldnames: list[str]) -> dict[str, str | None]:
    columns: dict[str, str | None] = {
        "modes": find_column(
            fieldnames,
            ("elastic_modes", "number_of_elastic_modes", "active_elastic_modes", "modes"),
        ),
        "internal_max": find_column(
            fieldnames,
            (
                "internal_vector_error_vs_richest_in",
                "max_modal_truncation_error_in",
                "maximum_tip_displacement_error_in",
            ),
        ),
        "left_vector": find_column(
            fieldnames, ("left_vector_error_in", "left_tip_vector_error_in")
        ),
        "right_vector": find_column(
            fieldnames, ("right_vector_error_in", "right_tip_vector_error_in")
        ),
        "available_max": find_column(
            fieldnames, ("max_tip_vector_error_in", "maximum_tip_vector_error_in")
        ),
    }
    for side in ("left", "right"):
        for component in ("Ux", "Uy", "Uz"):
            columns[f"{side}_{component}"] = find_column(
                fieldnames,
                (f"{side}_{component}_in", f"{side}_tip_{component}_in"),
            )
    if columns["modes"] is None:
        raise ValueError(f"Cannot identify the elastic-mode column in: {fieldnames}")
    return columns


def load_truncation_data(csv_path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, str | None]]:
    with csv_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {csv_path}")
        columns = identify_columns(reader.fieldnames)
        rows = list(reader)
    if not rows:
        raise ValueError(f"CSV contains no data: {csv_path}")

    mode_column = columns["modes"]
    assert mode_column is not None
    modes = np.asarray([int(float(row[mode_column])) for row in rows], dtype=int)
    order = np.argsort(modes)
    modes = modes[order]
    rows = [rows[index] for index in order]

    displacement_columns = [
        columns[f"{side}_{component}"]
        for side in ("left", "right")
        for component in ("Ux", "Uy", "Uz")
    ]
    calculated_error: np.ndarray | None = None
    if all(column is not None for column in displacement_columns):
        reference_indices = np.flatnonzero(modes == 54)
        if len(reference_indices) != 1:
            raise ValueError("Exactly one 54-mode reference case is required")
        reference_index = int(reference_indices[0])
        left_columns = [columns[f"left_{component}"] for component in ("Ux", "Uy", "Uz")]
        right_columns = [columns[f"right_{component}"] for component in ("Ux", "Uy", "Uz")]
        left = np.asarray(
            [[float(row[column]) for column in left_columns] for row in rows], dtype=float
        )
        right = np.asarray(
            [[float(row[column]) for column in right_columns] for row in rows], dtype=float
        )
        left_error = np.linalg.norm(left - left[reference_index], axis=1)
        right_error = np.linalg.norm(right - right[reference_index], axis=1)
        calculated_error = np.maximum(left_error, right_error)

    internal_column = columns["internal_max"]
    if internal_column is not None:
        errors = np.asarray([float(row[internal_column]) for row in rows], dtype=float)
        if calculated_error is not None and not np.allclose(
            errors, calculated_error, rtol=2.0e-5, atol=1.0e-10
        ):
            difference = float(np.max(np.abs(errors - calculated_error)))
            raise ValueError(
                f"Stored and recomputed truncation errors disagree (max difference {difference:.6g} in)"
            )
    elif calculated_error is not None:
        errors = calculated_error
    else:
        raise ValueError(
            "The CSV has neither an internal maximum error nor the six tip-displacement columns"
        )
    return modes, errors, columns


def select_basis(
    modes: np.ndarray, errors: np.ndarray, threshold: float
) -> tuple[int | None, float | None, bool]:
    nonreference = modes != 54
    tested_modes = modes[nonreference]
    tested_errors = errors[nonreference]
    for index, (mode, error) in enumerate(zip(tested_modes, tested_errors)):
        if error > threshold:
            continue
        richer_errors = tested_errors[index + 1 :]
        confirmed = len(richer_errors) >= 2 and bool(np.all(richer_errors[:2] <= threshold))
        return int(mode), float(error), confirmed
    return None, None, False


def generate_plot(
    csv_path: Path,
    png_path: Path,
    pdf_path: Path,
    threshold: float = 0.01,
    print_summary: bool = True,
) -> tuple[int | None, float | None, bool]:
    modes, errors, columns = load_truncation_data(csv_path)
    selected_mode, selected_error, confirmed = select_basis(modes, errors, threshold)

    print(f"CSV mode column: {columns['modes']}")
    print(f"CSV left vector-error column: {columns['left_vector'] or 'not present'}")
    print(f"CSV right vector-error column: {columns['right_vector'] or 'not present'}")
    print(f"CSV available maximum-error column: {columns['available_max'] or 'not present'}")
    print(
        "Modal-truncation error source: "
        + (columns["internal_max"] or "recomputed from left/right tip displacement vectors")
    )

    plot_modes = modes
    plot_errors = errors

    fig, axis = plt.subplots(figsize=(9.2, 6.0))
    axis.plot(
        plot_modes,
        plot_errors,
        color=THESIS_COLORS[0],
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
        alpha=0.9,
        zorder=1,
    )
    if selected_mode is not None and selected_error is not None:
        axis.axvline(
            selected_mode,
            color="0.45",
            linewidth=1.2,
            linestyle="-",
            alpha=0.9,
            zorder=1,
        )
        axis.scatter(
            [selected_mode], [selected_error], s=82, marker="o",
            color=THESIS_COLORS[4], edgecolor=THESIS_COLORS[0],
            linewidth=1.2, zorder=5,
        )
    axis.set_xlabel("Number of elastic modes, N")
    axis.set_ylabel(r"$e_N$ [in]")
    if len(modes) > 20:
        axis.set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 54])
    else:
        axis.set_xticks(modes)
    axis.set_xlim(0, 56)
    axis.set_ylim(0.0, max(float(np.max(plot_errors)) * 1.08, threshold * 2.0))
    axis.grid(False)
    axis.yaxis.grid(True, which="major", color="0.82", linewidth=0.7, alpha=0.65)
    axis.set_axisbelow(True)
    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=450, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    if print_summary:
        if selected_mode is None or selected_error is None:
            print(f"WARNING: no tested basis satisfies the {threshold:.4g} in threshold.")
        else:
            print(f"First basis satisfying {threshold:.4g} in: {selected_mode} elastic modes")
            print(f"Corresponding error: {selected_error:.8f} in")
            print(f"At least two subsequent bases remain below threshold: {'yes' if confirmed else 'no'}")
    return selected_mode, selected_error, confirmed


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv", nargs="?", type=Path,
        default=root / "plots/gravity_5g_convergence.csv",
    )
    parser.add_argument(
        "--png", type=Path,
        default=root / "plots/gravity_5g_modal_truncation_error.png",
    )
    parser.add_argument(
        "--pdf", type=Path,
        default=root / "plots/gravity_5g_modal_truncation_error.pdf",
    )
    parser.add_argument("--threshold", type=float, default=0.01)
    args = parser.parse_args()
    generate_plot(args.csv, args.png, args.pdf, args.threshold)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
