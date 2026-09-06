#!/usr/bin/env python3
"""Pair lumped/coupled-mass modes by MAC and compare their frequencies."""

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
from scipy.optimize import linear_sum_assignment


SPAN_NODES = list(range(990001, 990024)) + list(range(991002, 991024))
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


def nastran_float(value: str) -> float:
    value = value.replace("D", "E").replace("d", "E")
    if "E" not in value.upper():
        value = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", value)
    return float(value)


def parse_f06(path: Path) -> tuple[dict[int, float], dict[int, np.ndarray]]:
    frequencies: dict[int, float] = {}
    node_shapes: dict[int, dict[int, np.ndarray]] = {}
    current_mode: int | None = None
    in_eigenvalues = False
    fatal: list[str] = []

    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            compact = re.sub(r"\s+", "", line).upper()
            if "FATAL MESSAGE" in line:
                fatal.append(line.strip())
            if "REALEIGENVALUES" in compact:
                in_eigenvalues = True
                current_mode = None
                continue
            if in_eigenvalues:
                fields = line.split()
                if len(fields) >= 7 and fields[0].isdigit() and fields[1].isdigit():
                    frequencies[int(fields[0])] = nastran_float(fields[4])
                    continue
                if "USER INFORMATION" in line or "EIGENVALUE=" in compact:
                    in_eigenvalues = False

            match = re.search(r"EIGENVECTORNO\.?([0-9]+)", compact)
            if match:
                current_mode = int(match.group(1))
                node_shapes.setdefault(current_mode, {})
                continue
            if current_mode is None:
                continue
            fields = line.split()
            if len(fields) < 8 or not fields[0].isdigit() or fields[1].upper() != "G":
                continue
            node_id = int(fields[0])
            if node_id in SPAN_NODES:
                node_shapes[current_mode][node_id] = np.asarray(
                    [nastran_float(value) for value in fields[2:8]], dtype=float
                )

    if fatal:
        raise RuntimeError(f"{path.name} contains a Nastran fatal message: {fatal[0]}")
    shapes = {
        mode: np.vstack([values[node] for node in SPAN_NODES])
        for mode, values in node_shapes.items()
        if all(node in values for node in SPAN_NODES)
    }
    if not shapes:
        raise RuntimeError(
            f"No complete printed eigenvectors in {path}; run the supplied deck with PRINT output."
        )
    return frequencies, shapes


def mac(left: np.ndarray, right: np.ndarray) -> float:
    a = left.ravel()
    b = right.ravel()
    denominator = float(np.vdot(a, a).real * np.vdot(b, b).real)
    return float(abs(np.vdot(a, b)) ** 2 / denominator) if denominator else 0.0


def global_metrics(shape: np.ndarray) -> tuple[float, float]:
    amplitude = np.linalg.norm(shape[:, :3], axis=1)
    maximum = float(np.max(amplitude))
    if maximum <= 0.0:
        return 0.0, 0.0
    participation = float(np.linalg.norm(amplitude) / (np.sqrt(len(amplitude)) * maximum))
    coverage = float(np.mean(amplitude >= 0.2 * maximum))
    return participation, coverage


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "lumped", nargs="?", type=Path,
        default=root / "MAIN/sol103_coupmass_lumped.f06"
    )
    parser.add_argument(
        "coupled", nargs="?", type=Path,
        default=root / "MAIN/sol103_coupmass_coupled.f06"
    )
    parser.add_argument("--output", type=Path, default=root / "results")
    parser.add_argument("--first-elastic-mode", type=int, default=7)
    args = parser.parse_args()

    freq_l, shapes_l = parse_f06(args.lumped)
    freq_c, shapes_c = parse_f06(args.coupled)
    modes_l = sorted(mode for mode in shapes_l if mode >= args.first_elastic_mode)
    modes_c = sorted(mode for mode in shapes_c if mode >= args.first_elastic_mode)
    mac_matrix = np.asarray(
        [[mac(shapes_l[left], shapes_c[right]) for right in modes_c] for left in modes_l]
    )
    frequency_gap = np.asarray(
        [
            [abs(freq_c[right] - freq_l[left]) / freq_l[left] for right in modes_c]
            for left in modes_l
        ]
    )
    # Prevent a visually similar local shape from being paired across a remote
    # frequency band. Inside the 15% window, MAC remains the dominant metric.
    matching_cost = 1.0 - mac_matrix + 0.5 * frequency_gap
    matching_cost[frequency_gap > 0.15] = 100.0 + frequency_gap[frequency_gap > 0.15]
    rows_index, columns_index = linear_sum_assignment(matching_cost)

    rows = []
    for i, j in sorted(zip(rows_index, columns_index), key=lambda pair: modes_l[pair[0]]):
        left, right = modes_l[i], modes_c[j]
        f_l, f_c = freq_l[left], freq_c[right]
        participation, coverage = global_metrics(shapes_l[left])
        paired_mac = mac_matrix[i, j]
        reliable = paired_mac >= 0.90
        rows.append(
            {
                "lumped_mode": left,
                "coupled_mode": right,
                "f_lumped_hz": f_l,
                "f_coupled_hz": f_c,
                "frequency_shift_percent": 100.0 * (f_c - f_l) / f_l,
                "MAC": paired_mac,
                "reliable_pair_MAC_ge_0p90": reliable,
                "span_participation_index": participation,
                "span_coverage_20_percent": coverage,
                "global_candidate": (
                    reliable and coverage >= 0.50 and participation >= 0.35
                ),
            }
        )

    args.output.mkdir(parents=True, exist_ok=True)
    csv_path = args.output / "coupmass_mode_matching.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    mode = np.asarray([row["lumped_mode"] for row in rows])
    shift = np.asarray([row["frequency_shift_percent"] for row in rows])
    mac_value = np.asarray([row["MAC"] for row in rows])
    coverage = 100.0 * np.asarray([row["span_coverage_20_percent"] for row in rows])
    reliable = mac_value >= 0.90
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)
    axes[0].plot(
        mode, np.where(reliable, shift, np.nan), marker="o", linewidth=1.8,
        linestyle="-", label="Reliable pairs (MAC >= 0.90)"
    )
    axes[0].plot(
        mode, np.where(~reliable, shift, np.nan), marker="x", linewidth=1.8,
        linestyle="-", label="Low-confidence pairs"
    )
    axes[0].axhline(0.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8)
    axes[0].set_ylabel("Frequency shift [%]")
    axes[0].legend()
    axes[1].plot(mode, mac_value, marker="o", linewidth=1.8, linestyle="-")
    axes[1].axhline(
        0.90, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8,
        label="MAC = 0.90"
    )
    axes[1].set_ylabel("MAC [-]")
    axes[1].legend()
    axes[2].plot(mode, coverage, marker="o", linewidth=1.8, linestyle="-")
    axes[2].axhline(
        50.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8,
        label="50% span coverage"
    )
    axes[2].set_ylabel("Span coverage [%]")
    axes[2].set_xlabel("Lumped-mass mode number")
    axes[2].legend()
    for axis in axes:
        axis.grid(True, alpha=0.3)
    fig.suptitle("COUPMASS sensitivity and global-mode indicators")
    fig.tight_layout()
    fig.savefig(args.output / "coupmass_comparison.png", dpi=180)
    plt.close(fig)
    print(f"Compared {len(rows)} elastic modes. Results: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
