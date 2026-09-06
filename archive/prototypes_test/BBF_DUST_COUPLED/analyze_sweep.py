#!/usr/bin/env python3
"""Analyze a coupled MBDyn/DUST sweep and create technical English plots."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_bbf_sweep_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.signal import butter, detrend, find_peaks, periodogram, sosfiltfilt
from scipy.stats import linregress


INCH_TO_M = 0.0254
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167
RAD_TO_DEG = 180.0 / math.pi
BASE_NODE = 990000
TRIM_JOINT = 23
MODAL_JOINT = 5

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

CONTROLLERS = {
    9301: ("Vertical-speed pitch command", 4.0),
    9307: ("Altitude PI pitch command", 6.0),
    9303: ("Pitch PI surface command", 8.0),
    9305: ("Pitch-rate damping command", 4.0),
}


@dataclass
class ModeEstimate:
    frequency_hz: float
    growth_rate: float
    damping_ratio: float
    fit_r2: float
    filtered: np.ndarray
    peak_indices: np.ndarray
    fitted_envelope: np.ndarray
    spectrum_frequency: np.ndarray
    spectrum_psd: np.ndarray


@dataclass
class Run:
    velocity: float
    path: Path
    time: np.ndarray
    trim_end: float
    excitation_start: float
    bff_start: float
    bff_end: float
    heave: np.ndarray
    vertical_speed: np.ndarray
    pitch: np.ndarray
    aoa: np.ndarray
    pitch_rate: np.ndarray
    modal_q: np.ndarray
    trim_force: np.ndarray
    trim_moment: np.ndarray
    trim_pitch_output: np.ndarray
    trim_surface_output: np.ndarray
    flaps: dict[int, np.ndarray]
    controllers: dict[int, np.ndarray]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze coupled MBDyn/DUST sweep results without running either solver."
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
        help="Directory containing V_*_mps/case.nc folders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Analysis destination; default: RESULTS/ANALYSIS.",
    )
    return parser.parse_args()


def nc_data(nc: Dataset, name: str) -> np.ndarray:
    if name not in nc.variables:
        raise KeyError(f"Missing NetCDF variable: {name}")
    return np.ma.filled(nc.variables[name][:], np.nan).astype(float)


def log_constant(log_path: Path, name: str, fallback: float) -> float:
    if not log_path.is_file():
        return fallback
    text = log_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        rf"(?m)^\s*const\s+real\s+{re.escape(name)}\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*$",
        text,
    )
    return fallback if match is None else float(match.group(1).replace("D", "E"))


def read_run(nc_path: Path) -> Run:
    log_path = nc_path.with_suffix(".log")
    velocity = log_constant(log_path, "V_INF", math.nan)
    trim_end = log_constant(log_path, "TRIM_END", 12.0)
    excitation_start = log_constant(log_path, "EXCITATION_START", 13.5)
    bff_start = log_constant(log_path, "BFF_WINDOW_START", 13.75)
    bff_end = log_constant(log_path, "BFF_WINDOW_END", 28.75)

    with Dataset(nc_path) as nc:
        time = nc_data(nc, "time")
        position = nc_data(nc, f"node.struct.{BASE_NODE}.X")
        velocity_xyz = nc_data(nc, f"node.struct.{BASE_NODE}.XP")
        attitude = nc_data(nc, f"node.struct.{BASE_NODE}.Phi")
        omega = nc_data(nc, f"node.struct.{BASE_NODE}.Omega")
        modal_q = nc_data(nc, f"elem.joint.{MODAL_JOINT}.a")[:, 0]
        trim_force = (
            nc_data(nc, f"elem.joint.{TRIM_JOINT}.f")[:, 2] * LBF_TO_N
        )
        trim_moment = (
            nc_data(nc, f"elem.joint.{TRIM_JOINT}.m")[:, 1] * LBFIN_TO_NM
        )
        trim_pitch_output = (
            nc_data(nc, "elem.loadable.9101.output").squeeze() * RAD_TO_DEG
        )
        trim_surface_output = (
            nc_data(nc, "elem.loadable.9102.output").squeeze() * RAD_TO_DEG
        )
        flaps = {
            label: nc_data(nc, f"elem.joint.{label}.Phi")[:, 1] * RAD_TO_DEG
            for label in FLAPS
        }
        controllers = {
            label: nc_data(nc, f"elem.loadable.{label}.output").squeeze()
            * RAD_TO_DEG
            for label in CONTROLLERS
        }

    heave = (position[:, 2] - position[0, 2]) * INCH_TO_M
    vertical_speed = velocity_xyz[:, 2] * INCH_TO_M
    pitch = attitude[:, 1] * RAD_TO_DEG
    pitch_rate = omega[:, 1] * RAD_TO_DEG
    aoa = pitch - np.degrees(np.arctan2(vertical_speed, velocity))

    return Run(
        velocity=velocity,
        path=nc_path,
        time=time,
        trim_end=trim_end,
        excitation_start=excitation_start,
        bff_start=bff_start,
        bff_end=bff_end,
        heave=heave,
        vertical_speed=vertical_speed,
        pitch=pitch,
        aoa=aoa,
        pitch_rate=pitch_rate,
        modal_q=modal_q,
        trim_force=trim_force,
        trim_moment=trim_moment,
        trim_pitch_output=trim_pitch_output,
        trim_surface_output=trim_surface_output,
        flaps=flaps,
        controllers=controllers,
    )


def mode_estimate(
    time: np.ndarray,
    signal: np.ndarray,
    frequency_band: tuple[float, float],
) -> ModeEstimate:
    """Estimate modal frequency and exponential growth from a band-passed signal."""
    dt = float(np.median(np.diff(time)))
    fs = 1.0 / dt
    cleaned = detrend(np.nan_to_num(signal, nan=np.nanmedian(signal)))
    frequency, psd = periodogram(cleaned, fs=fs, window="hann", scaling="density")
    valid = (frequency >= frequency_band[0]) & (frequency <= frequency_band[1])
    if not np.any(valid):
        raise RuntimeError("No spectral samples inside the requested modal band")
    peak_frequency = float(frequency[valid][np.argmax(psd[valid])])

    half_width = max(0.30, 0.22 * peak_frequency)
    low = max(frequency_band[0], peak_frequency - half_width)
    high = min(frequency_band[1], peak_frequency + half_width)
    if high - low < 0.20:
        low, high = frequency_band
    sos = butter(4, (low, high), btype="bandpass", fs=fs, output="sos")
    filtered = sosfiltfilt(sos, cleaned)

    minimum_spacing = max(1, int(0.60 * fs / max(peak_frequency, 0.1)))
    amplitude = np.abs(filtered)
    peaks, _ = find_peaks(amplitude, distance=minimum_spacing)
    edge = max(2, int(0.40 * fs))
    peaks = peaks[(peaks >= edge) & (peaks < len(time) - edge)]
    if len(peaks):
        floor = max(float(np.nanmax(amplitude)) * 1.0e-5, 1.0e-14)
        peaks = peaks[amplitude[peaks] > floor]

    fit = np.full_like(time, np.nan, dtype=float)
    growth_rate = damping_ratio = fit_r2 = math.nan
    if len(peaks) >= 5:
        regression = linregress(time[peaks], np.log(amplitude[peaks]))
        growth_rate = float(regression.slope)
        fit_r2 = float(regression.rvalue**2)
        fit = np.exp(regression.intercept + regression.slope * time)
        omega = 2.0 * math.pi * peak_frequency
        damping_ratio = float(
            -growth_rate / math.sqrt(omega * omega + growth_rate * growth_rate)
        )

    return ModeEstimate(
        frequency_hz=peak_frequency,
        growth_rate=growth_rate,
        damping_ratio=damping_ratio,
        fit_r2=fit_r2,
        filtered=filtered,
        peak_indices=peaks,
        fitted_envelope=fit,
        spectrum_frequency=frequency,
        spectrum_psd=psd,
    )


def save_figure(destination: Path) -> None:
    plt.tight_layout()
    plt.savefig(destination, dpi=160)
    plt.close()


def time_plot(
    run: Run,
    destination: Path,
    signal: np.ndarray,
    ylabel: str,
    title: str,
    zero_line: bool = False,
) -> None:
    plt.figure(figsize=(10, 5.5))
    plt.plot(run.time, signal, linewidth=1.25)
    if zero_line:
        plt.axhline(0.0, color="0.25", linestyle="--", linewidth=0.8)
    plt.axvline(run.trim_end, color="black", linestyle="--", linewidth=0.9)
    plt.axvspan(
        run.excitation_start, run.bff_start, color="tab:orange", alpha=0.20,
        label="Control-surface doublet",
    )
    plt.axvspan(
        run.bff_start, run.bff_end, color="tab:red", alpha=0.08,
        label="BFF observation window",
    )
    plt.xlabel("Time [s]")
    plt.ylabel(ylabel)
    plt.title(f"{title} — $V_\\infty$ = {run.velocity:g} m/s")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    save_figure(destination)


def per_run_plots(
    run: Run,
    destination: Path,
    pitch_mode: ModeEstimate,
    bending_mode: ModeEstimate,
    linear_valid: bool,
) -> None:
    destination.mkdir(parents=True, exist_ok=True)

    time_plot(
        run, destination / "01_trim_vertical_force.png", run.trim_force,
        "Vertical trim reaction [N]", "Automatic-trim vertical-force residual",
        zero_line=True,
    )
    time_plot(
        run, destination / "02_trim_pitching_moment.png", run.trim_moment,
        "Pitching-moment trim reaction [N m]",
        "Automatic-trim pitching-moment residual", zero_line=True,
    )
    time_plot(
        run, destination / "03_trim_pitch_solution.png", run.trim_pitch_output,
        "Trim pitch command [deg]", "Automatic trim pitch solution",
    )
    time_plot(
        run, destination / "04_trim_surface_solution.png",
        run.trim_surface_output, "Symmetric trim deflection [deg]",
        "Automatic trim control-surface solution",
    )
    time_plot(
        run, destination / "05_rigid_body_heave.png", run.heave,
        "Rigid-body heave relative to initial position [m]",
        "Altitude-hold response", zero_line=True,
    )
    time_plot(
        run, destination / "06_vertical_speed.png", run.vertical_speed,
        "Vertical speed [m/s]", "Rigid-body vertical-speed response",
        zero_line=True,
    )
    time_plot(
        run, destination / "07_pitch_attitude.png", run.pitch,
        "Pitch attitude [deg]", "Rigid-body pitch response",
    )
    time_plot(
        run, destination / "08_angle_of_attack.png", run.aoa,
        "Rigid-body angle of attack [deg]", "Angle-of-attack response",
    )
    time_plot(
        run, destination / "09_pitch_rate.png", run.pitch_rate,
        "Pitch rate [deg/s]", "Rigid-body pitch-rate response", zero_line=True,
    )
    time_plot(
        run, destination / "10_first_symmetric_bending.png", run.modal_q,
        "Modal coordinate $q_1$ [-]", "First symmetric wing-bending response",
    )

    window = (run.time >= run.bff_start + 0.02) & (run.time <= run.bff_end)
    test_time = run.time[window]
    plt.figure(figsize=(10, 5.5))
    plt.plot(
        test_time, np.abs(bending_mode.filtered), color="0.55",
        label="Absolute band-passed 1SWB response",
    )
    if len(bending_mode.peak_indices):
        indices = bending_mode.peak_indices
        plt.plot(
            test_time[indices], np.abs(bending_mode.filtered[indices]), "o",
            markersize=3.5, label="Envelope samples",
        )
    if np.any(np.isfinite(bending_mode.fitted_envelope)):
        plt.plot(
            test_time, bending_mode.fitted_envelope, "r--",
            label=f"Exponential fit, $R^2$={bending_mode.fit_r2:.3f}",
        )
    status = "linear response" if linear_valid else "nonlinear response: fit not valid"
    plt.xlabel("Time [s]")
    plt.ylabel("Band-passed modal amplitude [-]")
    plt.title(
        f"1SWB modal envelope — $V_\\infty$ = {run.velocity:g} m/s "
        f"({status})"
    )
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    save_figure(destination / "11_bending_envelope.png")

    plt.figure(figsize=(10, 5.5))
    pitch_psd = pitch_mode.spectrum_psd / max(
        np.nanmax(pitch_mode.spectrum_psd), 1.0e-30
    )
    bending_psd = bending_mode.spectrum_psd / max(
        np.nanmax(bending_mode.spectrum_psd), 1.0e-30
    )
    plt.semilogy(
        pitch_mode.spectrum_frequency, pitch_psd, label="Pitch-attitude PSD",
    )
    plt.semilogy(
        bending_mode.spectrum_frequency, bending_psd, label="1SWB-coordinate PSD",
    )
    plt.xlim(0.0, 6.0)
    plt.ylim(1.0e-10, 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized power spectral density [-]")
    plt.title(f"Modal frequency content — $V_\\infty$ = {run.velocity:g} m/s")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(loc="best")
    save_figure(destination / "12_modal_frequency_spectrum.png")

    plt.figure(figsize=(10, 5.5))
    for label, name in FLAPS.items():
        plt.plot(run.time, run.flaps[label], linewidth=1.0, label=name)
    plt.axvspan(run.excitation_start, run.bff_start, color="tab:orange", alpha=0.2)
    plt.axvspan(run.bff_start, run.bff_end, color="tab:red", alpha=0.08)
    plt.xlabel("Time [s]")
    plt.ylabel("Control-surface deflection [deg]")
    plt.title(
        f"Symmetric control-surface activity — $V_\\infty$ = "
        f"{run.velocity:g} m/s"
    )
    plt.grid(True, alpha=0.3)
    plt.legend(ncol=5, fontsize=8, loc="best")
    save_figure(destination / "13_control_surface_activity.png")

    plt.figure(figsize=(10, 5.5))
    for label, (name, _) in CONTROLLERS.items():
        plt.plot(run.time, run.controllers[label], linewidth=1.0, label=name)
    plt.axvspan(run.excitation_start, run.bff_start, color="tab:orange", alpha=0.2)
    plt.axvspan(run.bff_start, run.bff_end, color="tab:red", alpha=0.08)
    plt.xlabel("Time [s]")
    plt.ylabel("Controller output [deg]")
    plt.title(
        f"Longitudinal flight-controller activity — $V_\\infty$ = "
        f"{run.velocity:g} m/s"
    )
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, loc="best")
    save_figure(destination / "14_controller_activity.png")


def one_comparison_plot(
    destination: Path,
    velocity: np.ndarray,
    curves: list[tuple[np.ndarray, str, str]],
    ylabel: str,
    title: str,
    zero_line: bool = False,
) -> None:
    plt.figure(figsize=(10, 5.5))
    for values, label, style in curves:
        plt.plot(velocity, values, style, linewidth=1.3, markersize=5, label=label)
    if zero_line:
        plt.axhline(0.0, color="black", linestyle="--", linewidth=0.9)
    plt.xlabel("Physical freestream velocity [m/s]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    save_figure(destination)


def flutter_crossings(
    velocity: np.ndarray,
    damping: np.ndarray,
    eligible: np.ndarray,
) -> list[float]:
    crossings: list[float] = []
    for index in range(len(velocity) - 1):
        if not (eligible[index] and eligible[index + 1]):
            continue
        if abs(velocity[index + 1] - velocity[index] - 1.0) > 1.0e-9:
            continue
        first, second = damping[index], damping[index + 1]
        if first > 0.0 and second <= 0.0:
            fraction = first / (first - second)
            crossings.append(float(velocity[index] + fraction))
    return crossings


def main() -> None:
    args = arguments()
    results = args.results.resolve()
    output = (args.output or results / "ANALYSIS").resolve()
    output.mkdir(parents=True, exist_ok=True)

    paths = sorted(results.glob("V_*_mps/case.nc"))
    if not paths:
        raise FileNotFoundError(f"No V_###_mps/case.nc results found in {results}")

    rows: list[dict[str, object]] = []
    for nc_path in paths:
        run = read_run(nc_path)
        trim_window = (run.time >= run.trim_end - 1.0) & (run.time < run.trim_end)
        bff_window = (
            (run.time >= run.bff_start + 0.02) & (run.time <= run.bff_end)
        )
        recovery_window = run.time >= max(run.bff_end, run.time[-1] - 1.0)
        if np.count_nonzero(bff_window) < 100:
            raise RuntimeError(f"Incomplete BFF window in {nc_path}")

        trim_pitch = float(np.mean(run.trim_pitch_output[trim_window]))
        trim_surface = float(np.mean(run.trim_surface_output[trim_window]))
        pitch_deviation = float(
            np.max(np.abs(run.pitch[bff_window] - trim_pitch))
        )
        heave_range = float(np.ptp(run.heave[bff_window]))
        linear_valid = bool(pitch_deviation <= 10.0 and heave_range <= 5.0)
        complete = bool(run.time[-1] >= run.bff_end)

        test_time = run.time[bff_window]
        pitch_mode = mode_estimate(
            test_time, run.pitch[bff_window], (0.30, 2.00)
        )
        bending_mode = mode_estimate(
            test_time, run.modal_q[bff_window], (1.00, 5.00)
        )

        flap_peak_dynamic: dict[str, float] = {}
        pretest = (
            (run.time >= run.trim_end + 1.0)
            & (run.time < run.excitation_start)
        )
        for label, name in FLAPS.items():
            baseline = float(np.mean(run.flaps[label][pretest]))
            flap_peak_dynamic[name] = float(
                np.max(np.abs(run.flaps[label][bff_window] - baseline))
            )

        saturation_fraction: dict[str, float] = {}
        for label, (name, limit) in CONTROLLERS.items():
            saturation_fraction[name] = float(
                np.mean(np.abs(run.controllers[label][bff_window]) >= 0.995 * limit)
            )

        fit_eligible = bool(
            complete
            and linear_valid
            and np.isfinite(bending_mode.damping_ratio)
            and bending_mode.fit_r2 >= 0.50
        )
        if not linear_valid:
            verdict = "NONLINEAR_RIGID_BODY_DEPARTURE"
        elif not np.isfinite(bending_mode.growth_rate):
            verdict = "INSUFFICIENT_MODAL_PEAKS"
        elif bending_mode.fit_r2 < 0.50:
            verdict = "LOW_CONFIDENCE_MODAL_FIT"
        elif bending_mode.growth_rate > 0.0:
            verdict = "LINEAR_MODAL_GROWTH"
        else:
            verdict = "LINEAR_MODAL_DECAY"

        row: dict[str, object] = {
            "velocity_mps": run.velocity,
            "complete": complete,
            "linear_valid": linear_valid,
            "fit_eligible": fit_eligible,
            "verdict": verdict,
            "trim_pitch_deg": trim_pitch,
            "trim_surface_deg": trim_surface,
            "trim_fz_mean_n": float(np.mean(run.trim_force[trim_window])),
            "trim_fz_std_n": float(np.std(run.trim_force[trim_window])),
            "trim_my_mean_nm": float(np.mean(run.trim_moment[trim_window])),
            "trim_my_std_nm": float(np.std(run.trim_moment[trim_window])),
            "bff_heave_range_m": heave_range,
            "recovery_mean_heave_m": float(np.mean(run.heave[recovery_window])),
            "bff_pitch_deviation_deg": pitch_deviation,
            "bff_aoa_range_deg": float(np.ptp(run.aoa[bff_window])),
            "bff_q1_peak_to_peak": float(np.ptp(run.modal_q[bff_window])),
            "pitch_frequency_hz": pitch_mode.frequency_hz,
            "pitch_growth_rate_1ps": pitch_mode.growth_rate,
            "pitch_damping_ratio": pitch_mode.damping_ratio,
            "pitch_fit_r2": pitch_mode.fit_r2,
            "bending_frequency_hz": bending_mode.frequency_hz,
            "bending_growth_rate_1ps": bending_mode.growth_rate,
            "bending_damping_ratio": bending_mode.damping_ratio,
            "bending_fit_r2": bending_mode.fit_r2,
            "frequency_separation_hz": abs(
                bending_mode.frequency_hz - pitch_mode.frequency_hz
            ),
            "peak_inboard_surface_command_deg": max(
                flap_peak_dynamic["BFL"],
                flap_peak_dynamic["BFR"],
                flap_peak_dynamic["WFL1"],
                flap_peak_dynamic["WFR1"],
            ),
            "peak_excitation_surface_command_deg": max(
                flap_peak_dynamic["WFL4"], flap_peak_dynamic["WFR4"]
            ),
            "maximum_controller_saturation_fraction": max(
                saturation_fraction.values()
            ),
        }
        rows.append(row)

        per_run_destination = output / f"V_{int(round(run.velocity)):03d}_mps"
        per_run_plots(
            run, per_run_destination, pitch_mode, bending_mode, linear_valid
        )
        with (per_run_destination / "run_assessment.txt").open(
            "w", encoding="utf-8"
        ) as stream:
            stream.write(f"Velocity: {run.velocity:g} m/s\n")
            stream.write(f"Assessment: {verdict}\n")
            stream.write(f"Linear-response validity: {linear_valid}\n")
            stream.write(
                f"1SWB frequency: {bending_mode.frequency_hz:.6g} Hz\n"
            )
            stream.write(
                f"1SWB growth rate: {bending_mode.growth_rate:.6g} 1/s\n"
            )
            stream.write(
                f"1SWB damping ratio: {bending_mode.damping_ratio:.6g}\n"
            )
            stream.write(f"1SWB fit R^2: {bending_mode.fit_r2:.6g}\n")

        print(
            f"Analyzed {run.velocity:g} m/s: {verdict}, "
            f"f_1SWB={bending_mode.frequency_hz:.3f} Hz, "
            f"sigma={bending_mode.growth_rate:.4g} 1/s",
            flush=True,
        )

    summary_path = output / "sweep_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    velocity = np.asarray([float(row["velocity_mps"]) for row in rows])
    linear = np.asarray([bool(row["linear_valid"]) for row in rows])
    eligible = np.asarray([bool(row["fit_eligible"]) for row in rows])

    def values(name: str, valid_only: bool = False) -> np.ndarray:
        result = np.asarray([float(row[name]) for row in rows], dtype=float)
        if valid_only:
            result[~linear] = np.nan
        return result

    comparative = output / "COMPARATIVE"
    comparative.mkdir(exist_ok=True)

    bending_damping = values("bending_damping_ratio", valid_only=True)
    pitch_damping = values("pitch_damping_ratio", valid_only=True)
    one_comparison_plot(
        comparative / "01_vg_diagram.png",
        velocity,
        [
            (100.0 * bending_damping, "1SWB damping ratio", "o-"),
            (100.0 * pitch_damping, "Short-period damping ratio", "s-"),
        ],
        "Modal damping ratio [%]",
        "V–g diagram: closed-loop longitudinal modal damping",
        zero_line=True,
    )
    one_comparison_plot(
        comparative / "02_vf_diagram.png",
        velocity,
        [
            (
                values("bending_frequency_hz", valid_only=True),
                "First symmetric wing-bending branch",
                "o-",
            ),
            (
                values("pitch_frequency_hz", valid_only=True),
                "Rigid-body short-period spectral branch",
                "s-",
            ),
        ],
        "Modal frequency [Hz]",
        "V–f diagram: closed-loop longitudinal modal frequencies",
    )
    one_comparison_plot(
        comparative / "03_modal_growth_rate.png",
        velocity,
        [
            (
                values("bending_growth_rate_1ps", valid_only=True),
                "1SWB exponential rate",
                "o-",
            ),
            (
                values("pitch_growth_rate_1ps", valid_only=True),
                "Short-period exponential rate",
                "s-",
            ),
        ],
        "Exponential rate [1/s]",
        "Modal growth rate versus physical freestream velocity",
        zero_line=True,
    )
    one_comparison_plot(
        comparative / "04_trim_pitch.png",
        velocity,
        [(values("trim_pitch_deg"), "Trim pitch attitude", "o-")],
        "Trim pitch attitude [deg]",
        "Automatic trim pitch versus physical freestream velocity",
    )
    one_comparison_plot(
        comparative / "05_trim_surface_deflection.png",
        velocity,
        [(values("trim_surface_deg"), "Symmetric trim deflection", "o-")],
        "Symmetric trim deflection [deg]",
        "Automatic trim control-surface deflection versus physical freestream velocity",
    )
    one_comparison_plot(
        comparative / "06_altitude_hold_performance.png",
        velocity,
        [
            (values("bff_heave_range_m"), "BFF-window heave range", "o-"),
            (
                np.abs(values("recovery_mean_heave_m")),
                "Absolute recovery mean heave",
                "s-",
            ),
        ],
        "Rigid-body heave metric [m]",
        "Altitude-hold performance versus physical freestream velocity",
    )
    one_comparison_plot(
        comparative / "07_angle_of_attack_excursion.png",
        velocity,
        [(values("bff_aoa_range_deg"), "BFF-window AoA range", "o-")],
        "Angle-of-attack range [deg]",
        "Rigid-body angle-of-attack excursion versus physical freestream velocity",
    )
    one_comparison_plot(
        comparative / "08_control_surface_activity.png",
        velocity,
        [
            (
                values("peak_inboard_surface_command_deg"),
                "Peak inboard control command",
                "o-",
            ),
            (
                values("peak_excitation_surface_command_deg"),
                "Peak WF4 excitation command",
                "s-",
            ),
        ],
        "Peak dynamic deflection [deg]",
        "Control-surface activity versus physical freestream velocity",
    )
    one_comparison_plot(
        comparative / "09_modal_frequency_separation.png",
        velocity,
        [
            (
                values("frequency_separation_hz", valid_only=True),
                "Short-period to 1SWB separation",
                "o-",
            )
        ],
        "Frequency separation [Hz]",
        "Longitudinal modal-frequency separation versus physical freestream velocity",
    )
    one_comparison_plot(
        comparative / "10_modal_fit_quality.png",
        velocity,
        [(values("bending_fit_r2"), "1SWB exponential-fit quality", "o-")],
        "Coefficient of determination $R^2$ [-]",
        "1SWB damping-estimation quality",
    )

    crossings = flutter_crossings(velocity, bending_damping, eligible)
    nonlinear_velocities = velocity[~linear]
    assessment_path = output / "flutter_assessment.txt"
    with assessment_path.open("w", encoding="utf-8") as stream:
        stream.write("FLUTTER-SPEED ASSESSMENT\n")
        stream.write("========================\n\n")
        stream.write(
            "Criterion: zero crossing of the 1--5 Hz 1SWB damping ratio, using "
            "only complete, linear-response cases with R^2 >= 0.50 and "
            "consecutive 1 m/s samples.\n\n"
        )
        if len(crossings) == 1:
            stream.write(
                f"Estimated flutter speed: {crossings[0]:.4g} m/s "
                "(linear interpolation of the V-g zero crossing).\n"
            )
        elif len(crossings) > 1:
            stream.write(
                "Flutter speed is inconclusive: multiple admissible zero "
                f"crossings were found at {crossings} m/s.\n"
            )
        else:
            stream.write(
                "Flutter speed was not identified from a valid monotonic V-g "
                "zero crossing inside the tested range.\n"
            )
        if len(nonlinear_velocities):
            stream.write(
                "\nCases rejected because of nonlinear rigid-body departure "
                "(heave range > 5 m or pitch deviation > 10 deg): "
                + ", ".join(f"{value:g}" for value in nonlinear_velocities)
                + " m/s.\n"
            )
        stream.write(
            "\nA nonlinear loss of altitude or pitch is not classified as "
            "flutter. Inspect the per-run modal envelope, V-g diagram and "
            "fit-quality diagram before accepting a flutter boundary.\n"
        )

    print(f"\nAnalysis written to: {output}")
    print(f"Summary table: {summary_path}")
    print(f"Flutter assessment: {assessment_path}")


if __name__ == "__main__":
    main()
