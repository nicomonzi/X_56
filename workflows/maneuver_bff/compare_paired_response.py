#!/usr/bin/env python3
"""Compare steady and maneuver BFF growth after excited-shadow subtraction."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bff_paired_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.signal import butter, hilbert, sosfiltfilt


def read_mode7(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with Dataset(path) as data:
        time = np.asarray(data["time"][:]).squeeze()
        modal = np.asarray(data["elem.joint.5.a"][:])
    return time, modal[:, 0]


def paired_delta(shadow: Path, excited: Path) -> tuple[np.ndarray, np.ndarray]:
    ts, ys = read_mode7(shadow)
    te, ye = read_mode7(excited)
    if ts.shape != te.shape or not np.allclose(ts, te, atol=1e-9, rtol=0.0):
        raise ValueError(f"assi temporali diversi: {shadow} / {excited}")
    return ts, ye - ys


def growth_metrics(
    time: np.ndarray, signal: np.ndarray, start: float, end: float, frequency: float,
) -> tuple[dict, np.ndarray, np.ndarray]:
    use = (time >= start) & (time <= end)
    t = time[use]
    y = signal[use]
    if len(t) < 50:
        raise ValueError("finestra di confronto troppo corta")
    dt = float(np.median(np.diff(t)))
    nyquist = 0.5 / dt
    low = max(0.2, 0.55 * frequency)
    high = min(1.65 * frequency, 0.90 * nyquist)
    if not low < high:
        raise ValueError("banda BFF incompatibile con il time step")
    sos = butter(4, (low, high), btype="bandpass", fs=1.0 / dt, output="sos")
    filtered = sosfiltfilt(sos, y)
    envelope = np.abs(hilbert(filtered))
    floor = max(0.05 * float(np.max(envelope)), 1e-12)
    period = 1.0 / frequency
    # Exclude Hilbert/filter edge transients from the logarithmic slope.
    fit = (envelope >= floor) & (t >= t[0] + period) & (t <= t[-1] - period)
    if np.count_nonzero(fit) < 20:
        raise ValueError("perturbazione o finestra insufficiente per stimare sigma")
    coefficients = np.polyfit(t[fit] - t[fit][0], np.log(np.maximum(envelope[fit], floor)), 1)
    sigma = float(coefficients[0])
    early = t <= t[0] + period
    late = t >= t[-1] - period
    early_energy = float(np.mean(envelope[early] ** 2))
    late_energy = float(np.mean(envelope[late] ** 2))
    return ({
        "sigma_per_s": sigma,
        "early_cycle_energy_proxy": early_energy,
        "late_cycle_energy_proxy": late_energy,
        "energy_growth_ratio": late_energy / max(early_energy, 1e-24),
    }, t, envelope)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steady-shadow", type=Path, required=True)
    parser.add_argument("--steady-excited", type=Path, required=True)
    parser.add_argument("--maneuver-shadow", type=Path, required=True)
    parser.add_argument("--maneuver-excited", type=Path, required=True)
    parser.add_argument("--start", type=float, required=True)
    parser.add_argument("--end", type=float, required=True)
    parser.add_argument("--frequency", type=float, default=2.0)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ts, steady = paired_delta(args.steady_shadow, args.steady_excited)
    tm, maneuver = paired_delta(args.maneuver_shadow, args.maneuver_excited)
    if ts.shape != tm.shape or not np.allclose(ts, tm, atol=1e-9, rtol=0.0):
        raise ValueError("steady e maneuver devono usare lo stesso asse temporale")
    steady_metrics, time, env_steady = growth_metrics(ts, steady, args.start, args.end, args.frequency)
    maneuver_metrics, _, env_maneuver = growth_metrics(tm, maneuver, args.start, args.end, args.frequency)
    delta_sigma = maneuver_metrics["sigma_per_s"] - steady_metrics["sigma_per_s"]
    relative_energy_growth = (
        maneuver_metrics["energy_growth_ratio"] / max(steady_metrics["energy_growth_ratio"], 1e-24)
    )
    sigma_resolution = 0.05
    verdict = (
        "amplified" if delta_sigma > sigma_resolution
        else "suppressed" if delta_sigma < -sigma_resolution
        else "indistinguishable_with_current_window"
    )
    summary = {
        "comparison_definition": "mode7 excited-minus-shadow, identical perturbation and time grid",
        "window_s": [args.start, args.end],
        "filter_center_frequency_hz": args.frequency,
        "steady": steady_metrics,
        "maneuver": maneuver_metrics,
        "delta_sigma_per_s": delta_sigma,
        "maneuver_to_steady_energy_growth_ratio": relative_energy_growth,
        "sigma_resolution_per_s": sigma_resolution,
        "verdict": verdict,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "paired_growth_summary.json").write_text(json.dumps(summary, indent=2))
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    use = (ts >= args.start) & (ts <= args.end)
    axes[0].plot(ts[use], steady[use], label="steady excited-shadow")
    axes[0].plot(tm[use], maneuver[use], label="maneuver excited-shadow", alpha=0.85)
    axes[0].set_ylabel("delta q7")
    axes[0].legend()
    axes[1].semilogy(time, np.maximum(env_steady, 1e-12), label="steady envelope")
    axes[1].semilogy(time, np.maximum(env_maneuver, 1e-12), label="maneuver envelope")
    axes[1].set(xlabel="time [s]", ylabel="BFF envelope")
    axes[1].legend()
    for axis in axes:
        axis.grid(True, alpha=0.3)
    fig.suptitle(f"{verdict}: delta sigma = {delta_sigma:+.4f} 1/s")
    fig.tight_layout()
    fig.savefig(args.output / "paired_growth_comparison.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
