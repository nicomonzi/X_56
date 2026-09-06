#!/usr/bin/env python3
"""Plot SOL 103 frequencies for lumped and selected coupled mass matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

from compare_coupmass import parse_f06


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


def ordered_frequencies(path: Path) -> tuple[np.ndarray, np.ndarray]:
    frequencies, _ = parse_f06(path)
    modes = np.asarray(sorted(frequencies), dtype=int)
    values = np.asarray([frequencies[mode] for mode in modes], dtype=float)
    return modes, values


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lumped",
        type=Path,
        default=root / "MAIN/sol103_coupmass_lumped.f06",
    )
    parser.add_argument(
        "--coupled",
        type=Path,
        default=root / "MAIN/sol103_coupmass_coupled.f06",
    )
    parser.add_argument("--output", type=Path, default=root / "results")
    args = parser.parse_args()

    lumped_modes, lumped_frequency = ordered_frequencies(args.lumped)
    coupled_modes, coupled_frequency = ordered_frequencies(args.coupled)
    elastic_lumped = lumped_modes > 6
    elastic_coupled = coupled_modes > 6
    lumped_modes = lumped_modes[elastic_lumped]
    lumped_frequency = lumped_frequency[elastic_lumped]
    coupled_modes = coupled_modes[elastic_coupled]
    coupled_frequency = coupled_frequency[elastic_coupled]

    common_modes = np.intersect1d(lumped_modes, coupled_modes)
    lumped_by_mode = dict(zip(lumped_modes, lumped_frequency))
    coupled_by_mode = dict(zip(coupled_modes, coupled_frequency))
    relative_difference = np.asarray([
        abs(coupled_by_mode[mode] - lumped_by_mode[mode])
        / lumped_by_mode[mode]
        for mode in common_modes
    ])
    separation_threshold = 0.01
    separated = np.flatnonzero(relative_difference > separation_threshold)
    separation_mode = int(common_modes[separated[0]]) if len(separated) else None

    fig, axis = plt.subplots(figsize=(9.2, 6.0))
    axis.plot(
        lumped_modes,
        lumped_frequency,
        color=THESIS_COLORS[0],
        marker="o",
        markersize=4.8,
        linewidth=1.8,
        linestyle="-",
        label="Lumped mass (COUPMASS = -1)",
    )
    axis.plot(
        coupled_modes,
        coupled_frequency,
        color=THESIS_COLORS[1],
        marker="s",
        markersize=4.5,
        linewidth=1.8,
        linestyle="-",
        label="Coupled mass (COUPMASS = 1)",
    )
    axis.set_title("Natural frequency comparison", fontweight="normal")
    axis.set_xlabel("Mode number")
    axis.set_ylabel("Natural frequency [Hz]")
    axis.set_xlim(7, max(int(lumped_modes[-1]), int(coupled_modes[-1])))
    axis.set_xticks([7, 10, 20, 30, 40, 50, 60])
    axis.grid(False)
    axis.yaxis.grid(True, color="0.82", linewidth=0.7, alpha=0.65)
    axis.set_axisbelow(True)
    axis.legend(frameon=True, fancybox=False, edgecolor="0.25")

    if separation_mode is not None:
        axis.axvline(
            separation_mode,
            color="0.25",
            linewidth=1.2,
            linestyle="-",
            alpha=0.75,
            zorder=1,
        )
        zoom = axis.inset_axes([0.53, 0.12, 0.43, 0.42])
        zoom.plot(
            lumped_modes,
            lumped_frequency,
            color=THESIS_COLORS[0],
            marker="o",
            markersize=4.8,
            linewidth=1.8,
            linestyle="-",
        )
        zoom.plot(
            coupled_modes,
            coupled_frequency,
            color=THESIS_COLORS[1],
            marker="s",
            markersize=4.5,
            linewidth=1.8,
            linestyle="-",
        )
        zoom.axvline(
            separation_mode,
            color="0.25",
            linewidth=1.2,
            linestyle="-",
            alpha=0.75,
        )
        zoom_start = max(7, separation_mode - 3)
        zoom_end = min(int(common_modes[-1]), separation_mode + 4)
        zoom_mask_l = (lumped_modes >= zoom_start) & (lumped_modes <= zoom_end)
        zoom_mask_c = (coupled_modes >= zoom_start) & (coupled_modes <= zoom_end)
        zoom_values = np.concatenate(
            [lumped_frequency[zoom_mask_l], coupled_frequency[zoom_mask_c]]
        )
        zoom_margin = 0.08 * float(np.ptp(zoom_values))
        zoom.set_xlim(zoom_start - 0.3, zoom_end + 0.3)
        zoom.set_ylim(
            float(np.min(zoom_values)) - zoom_margin,
            float(np.max(zoom_values)) + zoom_margin,
        )
        zoom.set_xticks(np.arange(zoom_start, zoom_end + 1))
        zoom.tick_params(axis="both", labelsize=10)
        zoom.grid(False)
        zoom.yaxis.grid(True, color="0.85", linewidth=0.6, alpha=0.6)
        zoom.set_title(
            f"First separation above 1%: mode {separation_mode}",
            fontsize=11,
            fontweight="normal",
        )
    fig.tight_layout()

    args.output.mkdir(parents=True, exist_ok=True)
    png = args.output / "coupmass_frequency_comparison.png"
    pdf = args.output / "coupmass_frequency_comparison.pdf"
    fig.savefig(png, dpi=450, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Lumped modes plotted: {len(lumped_modes)}")
    print(f"COUPMASS = 1 modes plotted: {len(coupled_modes)}")
    if separation_mode is not None:
        index = int(np.where(common_modes == separation_mode)[0][0])
        print(
            f"First same-index frequency separation above 1%: mode {separation_mode} "
            f"({100.0 * relative_difference[index]:.3f}%)"
        )
    print(f"PNG: {png}")
    print(f"PDF: {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
