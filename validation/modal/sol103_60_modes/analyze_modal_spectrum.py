#!/usr/bin/env python3
"""Read the SOL 103 F06 and plot the extracted modal spectrum."""

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
    return float(value.replace("D", "E").replace("d", "E"))


def parse_frequencies(path: Path) -> dict[int, float]:
    frequencies: dict[int, float] = {}
    in_table = False
    fatal: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "FATAL MESSAGE" in line:
                fatal.append(line.strip())
            if "REALEIGENVALUES" in re.sub(r"\s+", "", line).upper():
                in_table = True
                continue
            if not in_table:
                continue
            fields = line.split()
            if len(fields) >= 7 and fields[0].isdigit() and fields[1].isdigit():
                frequencies[int(fields[0])] = nastran_float(fields[4])
            elif "USER INFORMATION" in line:
                in_table = False
    if fatal:
        raise RuntimeError(f"Nastran fatal message: {fatal[0]}")
    if not frequencies:
        raise RuntimeError("No real-eigenvalue summary found in the F06")
    return frequencies


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "f06", nargs="?", type=Path, default=root / "MAIN/sol103_60_modes.f06"
    )
    parser.add_argument("--output", type=Path, default=root / "results")
    args = parser.parse_args()

    frequencies = parse_frequencies(args.f06)
    modes = np.asarray(sorted(frequencies))
    values = np.asarray([frequencies[mode] for mode in modes])
    elastic = modes >= 7
    spacing = np.diff(values[elastic])
    args.output.mkdir(parents=True, exist_ok=True)

    with (args.output / "modal_frequencies_60.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(("mode", "frequency_hz", "type"))
        for mode, frequency in zip(modes, values):
            writer.writerow((mode, frequency, "rigid" if mode <= 6 else "elastic"))

    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=False)
    axes[0].plot(
        modes[elastic], values[elastic], marker="o", markersize=4,
        linewidth=1.8, linestyle="-", label="Elastic modes"
    )
    axes[0].axhline(
        0.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8,
        label="Zero frequency"
    )
    axes[0].set_xlabel("Nastran mode number")
    axes[0].set_ylabel("Natural frequency [Hz]")
    axes[0].set_title("SOL 103 modal spectrum: 54 elastic modes")
    axes[0].legend()

    axes[1].plot(
        modes[elastic][1:], spacing, marker="s", markersize=4,
        linewidth=1.8, linestyle="-", label="Adjacent-mode spacing"
    )
    axes[1].axhline(
        0.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8
    )
    axes[1].set_xlabel("Upper mode number of each pair")
    axes[1].set_ylabel("Frequency spacing [Hz]")
    axes[1].set_title("Modal-density indicator")
    axes[1].legend()
    for axis in axes:
        axis.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output / "modal_spectrum_60.png", dpi=180)
    plt.close(fig)

    print(f"Modes found: {len(modes)} total, {int(np.sum(elastic))} elastic")
    print(f"Elastic frequency range: {values[elastic][0]:.6g} to {values[elastic][-1]:.6g} Hz")
    print(f"Smallest adjacent elastic spacing: {np.min(spacing):.6g} Hz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
