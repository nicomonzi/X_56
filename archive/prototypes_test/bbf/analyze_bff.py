#!/usr/bin/env python3
"""Post-process one NASA X-56 BFF run without launching MBDyn.

Usage:
    python3 analyze_bff.py output/test_bff.nc
    python3 analyze_bff.py output/test_bff

The script reads the named NetCDF variables actually written by MBDyn, saves
PNG figures beside the run, and prints quantitative frequency/envelope/trim
metrics.  It deliberately does not label a case as stable or fluttering.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_bff_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

try:
    from scipy.signal import find_peaks
except ImportError:  # Small fallback; SciPy is convenient but not mandatory.
    find_peaks = None


INCH_TO_M = 0.0254
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167
RAD_TO_DEG = 180.0 / math.pi

BASE_NODE = 990000
MODAL_JOINT = 5
FINAL_JOINT = 1
TRIM_JOINT = 23
FLAPS = {
    1004: "BFL",
    1008: "WFL1",
    1011: "WFL2",
    1014: "WFL3",
    1017: "WFL4",
    2004: "BFR",
    2008: "WFR1",
    2011: "WFR2",
    2014: "WFR3",
    2017: "WFR4",
}
EXCITED = {1014, 1017, 2014, 2017}
FLIGHT_CONTROLLED = {1004, 2004}
FLIGHT_CONTROL_PIDS = {
    9301: "Altitude/Vz PI",
    9303: "Pitch P",
    9305: "Filtered pitch-rate damper",
    9307: "Heave-error integrator",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot and quantify one MBDyn BFF run")
    parser.add_argument("result", nargs="?", default="output/test_bff.nc")
    return parser.parse_args()


def result_path(value: str) -> Path:
    path = Path(value)
    if path.suffix != ".nc":
        path = path.with_suffix(".nc")
    if not path.is_file():
        raise FileNotFoundError(f"MBDyn NetCDF result not found: {path}")
    return path


def main_constant(name: str, fallback: float) -> float:
    text = (Path(__file__).resolve().parent / "main_x56_bff.mbd").read_text()
    match = re.search(
        rf"set\s*:\s*const\s+real\s+{re.escape(name)}\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*;",
        text,
        flags=re.IGNORECASE,
    )
    return fallback if match is None else float(match.group(1).replace("D", "E"))


def variable(nc: Dataset, name: str) -> np.ndarray:
    if name not in nc.variables:
        available = "\n  ".join(sorted(nc.variables))
        raise KeyError(f"Required MBDyn variable '{name}' is absent. Available:\n  {available}")
    return np.ma.filled(nc.variables[name][:], np.nan).astype(float)


def optional_variable(nc: Dataset, name: str) -> np.ndarray | None:
    if name not in nc.variables:
        return None
    return np.ma.filled(nc.variables[name][:], np.nan).astype(float)


def detrend_linear(time: np.ndarray, signal: np.ndarray) -> np.ndarray:
    good = np.isfinite(signal)
    if np.count_nonzero(good) < 2:
        return signal - np.nanmean(signal)
    coefficients = np.polyfit(time[good], signal[good], 1)
    return signal - np.polyval(coefficients, time)


def remove_slow_component(time: np.ndarray, signal: np.ndarray, window_seconds: float = 0.8):
    """Remove quasi-static drift before estimating the elastic-mode envelope."""
    dt = float(np.median(np.diff(time)))
    count = max(3, int(round(window_seconds / dt)))
    if count % 2 == 0:
        count += 1
    pad = count // 2
    extended = np.pad(signal, pad, mode="edge")
    trend = np.convolve(extended, np.ones(count) / count, mode="valid")
    return signal - trend


def spectrum(time: np.ndarray, signal: np.ndarray, min_frequency: float = 0.2):
    good = np.isfinite(time) & np.isfinite(signal)
    time = time[good]
    signal = signal[good]
    if len(time) < 8:
        return np.array([]), np.array([]), math.nan
    dt = float(np.median(np.diff(time)))
    y = detrend_linear(time, signal)
    window = np.hanning(len(y))
    amplitude = 2.0 * np.abs(np.fft.rfft(y * window)) / max(np.sum(window), 1.0)
    frequency = np.fft.rfftfreq(len(y), dt)
    eligible = frequency >= min_frequency
    dominant = frequency[eligible][np.argmax(amplitude[eligible])] if np.any(eligible) else math.nan
    return frequency, amplitude, float(dominant)


def envelope_fit(time: np.ndarray, signal: np.ndarray, dominant_frequency: float):
    """Fit log amplitudes of local |signal| peaks; result is growth rate [1/s]."""
    if len(time) < 3:
        return math.nan, np.array([], dtype=int), np.full_like(signal, np.nan)
    y = np.abs(detrend_linear(time, signal))
    dt = float(np.median(np.diff(time)))
    if not np.isfinite(dominant_frequency):
        dominant_frequency = 0.1
    distance = max(1, int(0.65 / (max(dominant_frequency, 0.1) * dt)))
    if find_peaks is not None:
        indices, _ = find_peaks(y, distance=distance)
    else:
        indices = np.flatnonzero((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])) + 1
        if len(indices):
            kept = [int(indices[0])]
            for index in indices[1:]:
                if index - kept[-1] >= distance:
                    kept.append(int(index))
            indices = np.asarray(kept)
    floor = max(np.nanmax(y) * 1.0e-6, np.finfo(float).eps)
    indices = indices[np.isfinite(y[indices]) & (y[indices] > floor)]
    if len(indices) < 3:
        return math.nan, indices, np.full_like(y, np.nan)
    slope, intercept = np.polyfit(time[indices], np.log(y[indices]), 1)
    fitted = np.exp(intercept + slope * time)
    return float(slope), indices, fitted


def normalized(signal: np.ndarray) -> np.ndarray:
    scale = np.nanstd(signal)
    return (signal - np.nanmean(signal)) / scale if scale > 0.0 else np.zeros_like(signal)


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def run() -> None:
    path = result_path(arguments().result)
    plot_dir = path.parent / f"{path.stem}_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    release_time = main_constant("RELEASE_TIME", 3.0)
    excitation_start = main_constant("EXCITATION_START", 4.0)
    excitation_frequency = main_constant("EXCITATION_FREQUENCY", 3.21713438486)
    excitation_end = excitation_start + 1.0 / excitation_frequency
    velocity = main_constant("V_INF", math.nan)

    with Dataset(path) as nc:
        time = variable(nc, "time")
        base_x = variable(nc, f"node.struct.{BASE_NODE}.X")
        base_phi = variable(nc, f"node.struct.{BASE_NODE}.Phi")
        base_v = variable(nc, f"node.struct.{BASE_NODE}.XP")
        base_omega = variable(nc, f"node.struct.{BASE_NODE}.Omega")
        modal_q = variable(nc, f"elem.joint.{MODAL_JOINT}.a")
        modal_qp = variable(nc, f"elem.joint.{MODAL_JOINT}.aPrime")
        final_force = variable(nc, f"elem.joint.{FINAL_JOINT}.f")
        final_moment = variable(nc, f"elem.joint.{FINAL_JOINT}.m")
        trim_force = variable(nc, f"elem.joint.{TRIM_JOINT}.f")
        trim_moment = variable(nc, f"elem.joint.{TRIM_JOINT}.m")
        flap = {
            label: variable(nc, f"elem.joint.{label}.Phi")[:, 1] * RAD_TO_DEG
            for label in FLAPS
        }
        flight_control = {
            label: optional_variable(nc, f"elem.loadable.{label}.output")
            for label in FLIGHT_CONTROL_PIDS
        }
        tip_l = variable(nc, "node.struct.990020.X")[:, 2]
        tip_r = variable(nc, "node.struct.991020.X")[:, 2]

    # Z is reported relative to the initial CG position: Z(0) = 0 by definition.
    z = (base_x[:, 2] - base_x[0, 2]) * INCH_TO_M
    vz = base_v[:, 2] * INCH_TO_M
    az = np.gradient(vz, time)
    pitch = base_phi[:, 1] * RAD_TO_DEG
    pitch_rate = base_omega[:, 1] * RAD_TO_DEG
    pitch_acceleration = np.gradient(pitch_rate, time)
    bending = modal_q[:, 0]
    bending_rate = modal_qp[:, 0]
    tip_symmetric = ((tip_l + tip_r) / 2.0 - (tip_l[0] + tip_r[0]) / 2.0) * INCH_TO_M

    fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    rigid = [(z, "Z relative [m]"), (vz, "Vertical speed [m/s]"),
             (az, "Vertical acceleration [m/s²]"), (pitch, "Pitch [deg]"),
             (pitch_rate, "Pitch rate [deg/s]"),
             (pitch_acceleration, "Pitch acceleration [deg/s²]")]
    for axis, (signal, label) in zip(axes.flat, rigid):
        axis.plot(time, signal, lw=1.0)
        axis.axvline(release_time, color="0.5", ls="--", lw=0.8)
        axis.axvspan(excitation_start, excitation_end, color="tab:orange", alpha=0.2)
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.25)
    axes[-1, 0].set_xlabel("Time [s]")
    axes[-1, 1].set_xlabel("Time [s]")
    fig.suptitle(f"Rigid-body longitudinal response — V = {velocity:g} m/s")
    save_figure(fig, plot_dir / "rigid_body.png")

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(time, bending, label="FEM mode 7 / modal coordinate 1")
    axes[1].plot(time, bending_rate)
    axes[2].plot(time, tip_symmetric)
    for axis in axes:
        axis.axvspan(excitation_start, excitation_end, color="tab:orange", alpha=0.2)
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel("Modal q₁ [-]")
    axes[1].set_ylabel("Modal q̇₁ [1/s]")
    axes[2].set_ylabel("Mean tip ΔZ [m]")
    axes[2].set_xlabel("Time [s]")
    fig.suptitle("First symmetric bending response")
    save_figure(fig, plot_dir / "bending.png")

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for label, name in FLAPS.items():
        axis = axes[0] if label < 2000 else axes[1]
        axis.plot(time, flap[label], label=name, lw=1.0)
    for axis, side in zip(axes, ("Left", "Right")):
        axis.axvspan(excitation_start, excitation_end, color="tab:orange", alpha=0.2)
        axis.set_ylabel(f"{side} deflection [deg]")
        axis.grid(True, alpha=0.25)
        axis.legend(ncol=5, fontsize=8)
    axes[1].set_xlabel("Time [s]")
    fig.suptitle("Control-surface deflections")
    save_figure(fig, plot_dir / "flaps.png")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axes[0, 0].plot(time, trim_force[:, 2] * LBF_TO_N, label="Fz")
    axes[0, 0].set_ylabel("Trim-joint Fz [N]")
    axes[0, 1].plot(time, trim_moment[:, 1] * LBFIN_TO_NM, label="My")
    axes[0, 1].set_ylabel("Trim-joint My [N·m]")
    axes[1, 0].plot(time, final_force[:, 0] * LBF_TO_N, label="Fx")
    axes[1, 0].plot(time, final_force[:, 1] * LBF_TO_N, label="Fy")
    axes[1, 0].set_ylabel("Final-pin force [N]")
    axes[1, 1].plot(time, final_moment[:, 0] * LBFIN_TO_NM, label="Mx")
    axes[1, 1].plot(time, final_moment[:, 2] * LBFIN_TO_NM, label="Mz")
    axes[1, 1].set_ylabel("Final-pin moment [N·m]")
    for axis in axes.flat:
        axis.axvline(release_time, color="0.5", ls="--", lw=0.8)
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=8)
        axis.set_xlabel("Time [s]")
    fig.suptitle("Constraint reactions (trim-joint values are masked after release)")
    save_figure(fig, plot_dir / "reactions.png")

    post = time >= excitation_end + 0.15
    if np.count_nonzero(post) < 8:
        # Useful diagnostic fallback for deliberately shortened probe runs.
        post = time >= max(release_time, time[-1] - 0.5 * (time[-1] - time[0]))
    analysis_time = time[post]
    elastic_bending = remove_slow_component(analysis_time, bending[post])
    signals = {"Pitch": pitch[post], "Heave": z[post], "Bending q₁": elastic_bending}
    spectra = {
        name: spectrum(analysis_time, signal, 1.0 if name == "Bending q₁" else 0.2)
        for name, signal in signals.items()
    }
    modal_frequency = spectra["Bending q₁"][2]
    growth_rate, peaks, fitted_envelope = envelope_fit(
        analysis_time, signals["Bending q₁"], modal_frequency
    )

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    for name, signal in signals.items():
        axes[0].plot(analysis_time, normalized(detrend_linear(analysis_time, signal)), label=name)
    bending_detrended = signals["Bending q₁"]
    axes[1].plot(analysis_time, np.abs(bending_detrended), color="0.6", label="|detrended q₁|")
    if len(peaks):
        axes[1].plot(analysis_time[peaks], np.abs(bending_detrended)[peaks], "o", ms=4, label="peaks")
    if np.any(np.isfinite(fitted_envelope)):
        axes[1].plot(analysis_time, fitted_envelope, "r--", label="exp. fit")
    axes[0].set_ylabel("Normalized amplitude")
    axes[1].set_ylabel("Modal envelope [-]")
    axes[1].set_xlabel("Time [s]")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend()
    fig.suptitle("Pitch–heave–bending coupling and post-pulse envelope")
    save_figure(fig, plot_dir / "coupling_envelope.png")

    fig, axis = plt.subplots(figsize=(11, 6))
    for name, (frequency, amplitude, _) in spectra.items():
        if len(frequency):
            scale = np.nanmax(amplitude)
            axis.plot(frequency, amplitude / scale if scale > 0 else amplitude, label=name)
    axis.set_xlim(0.0, min(15.0, 0.5 / np.median(np.diff(time))))
    axis.set_xlabel("Frequency [Hz]")
    axis.set_ylabel("Normalized FFT amplitude")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.suptitle("Post-perturbation spectra")
    save_figure(fig, plot_dir / "spectra.png")

    trim_window = (time >= max(0.0, release_time - 1.0)) & (time < release_time)
    fz_trim = trim_force[trim_window, 2] * LBF_TO_N
    my_trim = trim_moment[trim_window, 1] * LBFIN_TO_NM
    free_window = time >= release_time
    settled_window = time >= max(release_time, time[-1] - 3.0)
    x_drift = (base_x[free_window, :2] - base_x[np.flatnonzero(free_window)[0], :2]) * INCH_TO_M
    roll_yaw = (base_phi[free_window][:, [0, 2]] - base_phi[np.flatnonzero(free_window)[0], [0, 2]]) * RAD_TO_DEG
    after_pulse = time >= excitation_end + 0.05
    if not np.any(after_pulse):
        after_pulse = time >= time[-1] - min(0.5, 0.1 * (time[-1] - time[0]))
    trim_only_variation = max(
        np.ptp(flap[label][after_pulse])
        for label in FLAPS
        if label not in EXCITED and label not in FLIGHT_CONTROLLED
    )
    body_flap_control_variation = max(
        np.ptp(flap[label][after_pulse]) for label in FLIGHT_CONTROLLED
    )
    excited_final_spread = max(
        np.ptp(flap[label][time >= time[-1] - 0.5]) for label in EXCITED
    )

    lines = [
        f"Result: {path}",
        f"Velocity: {velocity:.6g} m/s",
        f"Trim residual Fz mean/std: {np.nanmean(fz_trim):.6g} / {np.nanstd(fz_trim):.6g} N",
        f"Trim residual My mean/std: {np.nanmean(my_trim):.6g} / {np.nanstd(my_trim):.6g} N m",
        f"Post-pulse dominant pitch frequency: {spectra['Pitch'][2]:.6g} Hz",
        f"Post-pulse dominant heave frequency: {spectra['Heave'][2]:.6g} Hz",
        f"Post-pulse dominant bending frequency: {modal_frequency:.6g} Hz",
        f"Bending-envelope exponential rate: {growth_rate:.6g} 1/s",
        f"Maximum constrained X/Y drift: {np.nanmax(np.abs(x_drift)):.6g} m",
        f"Maximum constrained roll/yaw drift: {np.nanmax(np.abs(roll_yaw)):.6g} deg",
        f"Free heave range after release: {np.ptp(z[free_window]):.6g} m",
        f"Mean rigid heave after pulse: {np.nanmean(z[post]):.6g} m",
        f"Mean rigid heave in final 3 s: {np.nanmean(z[settled_window]):.6g} m",
        f"Free pitch range after release: {np.ptp(pitch[free_window]):.6g} deg",
        f"Trim-only wing-flap post-pulse variation: {trim_only_variation:.6g} deg",
        f"Controlled body-flap post-pulse variation: {body_flap_control_variation:.6g} deg",
        f"Excited flap variation in final 0.5 s: {excited_final_spread:.6g} deg",
        f"Figures: {plot_dir}",
    ]
    for label, name in FLIGHT_CONTROL_PIDS.items():
        signal = flight_control[label]
        if signal is not None:
            signal_deg = signal * RAD_TO_DEG
            lines.insert(
                -1,
                f"{name} command min/max: "
                f"{np.nanmin(signal_deg):.6g} / {np.nanmax(signal_deg):.6g} deg",
            )
    report = "\n".join(lines) + "\n"
    print(report, end="")
    (plot_dir / "analysis_summary.txt").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    run()
