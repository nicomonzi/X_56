#!/usr/bin/env python3
"""Build V-g and V-f diagrams from free responses using Hilbert envelopes."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_bbf_hilbert_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.signal import butter, detrend, hilbert, periodogram, sosfiltfilt


DEFAULT_RESULTS = Path(
    "/mnt/c/Users/Utente/Desktop/RESULTS_BBF_DT002_MODES_7_23_DENSE"
)
BASE_NODE = 990000
LEFT_TIP = 990020
RIGHT_TIP = 991020
MODAL_JOINT = 5
INCH_TO_M = 0.0254

# Il ramo BBF osservato nelle simulazioni a 17 modi e' vicino a 1.3 Hz.
SEARCH_BAND_HZ = (0.60, 2.00)
TRACK_HALF_WIDTH_HZ = 0.30
FREE_RESPONSE_DELAY_S = 0.50
FILTER_EDGE_S = 0.40

CONTROLLER_LIMITS_DEG = {
    9301: 4.0,
    9307: 6.0,
    9303: 8.0,
    9305: 4.0,
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create separate Hilbert-based V-g, V-f and sigma diagrams from "
            "all complete 17-mode simulations."
        )
    )
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Default: RESULTS/ANALYSIS_HILBERT.",
    )
    return parser.parse_args()


def nc_array(dataset: Dataset, name: str) -> np.ndarray:
    if name not in dataset.variables:
        raise KeyError(f"Missing NetCDF variable: {name}")
    return np.ma.filled(dataset.variables[name][:], np.nan).astype(float)


def log_constant(path: Path, name: str, fallback: float) -> float:
    if not path.is_file():
        return fallback
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        rf"(?m)^\s*const\s+real\s+{re.escape(name)}\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*$",
        text,
    )
    return fallback if match is None else float(match.group(1).replace("D", "E"))


def common_frequency(
    time: np.ndarray,
    pitch: np.ndarray,
    bending: np.ndarray,
) -> float:
    """Locate the spectral component shared by pitch and symmetric bending."""
    sample_rate = 1.0 / float(np.median(np.diff(time)))
    frequency, pitch_psd = periodogram(
        detrend(pitch), fs=sample_rate, window="hann"
    )
    _, bending_psd = periodogram(
        detrend(bending), fs=sample_rate, window="hann"
    )
    selected = (
        (frequency >= SEARCH_BAND_HZ[0])
        & (frequency <= SEARCH_BAND_HZ[1])
    )
    if not np.any(selected):
        return math.nan
    pitch_scale = max(float(np.max(pitch_psd[selected])), 1.0e-30)
    bending_scale = max(float(np.max(bending_psd[selected])), 1.0e-30)
    shared_power = (
        pitch_psd[selected] / pitch_scale
        + bending_psd[selected] / bending_scale
    )
    return float(frequency[selected][np.argmax(shared_power)])


def hilbert_parameters(
    time: np.ndarray,
    signal: np.ndarray,
    center_frequency_hz: float,
) -> tuple[float, float, float, int]:
    """Return sigma, damped frequency and g without a regression."""
    dt = float(np.median(np.diff(time)))
    sample_rate = 1.0 / dt
    low = max(SEARCH_BAND_HZ[0], center_frequency_hz - TRACK_HALF_WIDTH_HZ)
    high = min(SEARCH_BAND_HZ[1], center_frequency_hz + TRACK_HALF_WIDTH_HZ)
    if high - low < 0.20:
        return math.nan, math.nan, math.nan, 0

    filtered = sosfiltfilt(
        butter(4, (low, high), btype="bandpass", fs=sample_rate, output="sos"),
        detrend(signal),
    )
    analytic = hilbert(filtered)
    envelope = np.abs(analytic)
    phase = np.unwrap(np.angle(analytic))

    usable = (
        (time >= time[0] + FILTER_EDGE_S)
        & (time <= time[-1] - FILTER_EDGE_S)
    )
    # Conserva anche i casi prossimi al flutter: una soglia basata su un
    # percentile eliminerebbe le parti iniziali/finali di un inviluppo che
    # cresce o decade. Il 2% del massimo rimuove soltanto il fondo numerico.
    amplitude_floor = max(
        0.02 * float(np.max(envelope[usable])),
        1.0e-14,
    )
    usable &= envelope >= amplitude_floor

    # Differenze su circa due cicli: robuste al ripple dell'inviluppo e senza
    # regressione globale o coefficiente R^2.
    lag = max(1, int(round(2.0 * sample_rate / center_frequency_hz)))
    if len(time) <= lag:
        return math.nan, math.nan, math.nan, 0
    pair_valid = usable[:-lag] & usable[lag:]
    if np.count_nonzero(pair_valid) < 8:
        return math.nan, math.nan, math.nan, int(np.count_nonzero(pair_valid))

    elapsed = time[lag:] - time[:-lag]
    sigma_samples = (
        np.log(envelope[lag:]) - np.log(envelope[:-lag])
    ) / elapsed
    omega_samples = (phase[lag:] - phase[:-lag]) / elapsed
    sigma = float(np.median(sigma_samples[pair_valid]))
    omega_d = float(np.median(omega_samples[pair_valid]))
    frequency_hz = omega_d / (2.0 * math.pi)
    omega_n = math.sqrt(omega_d * omega_d + sigma * sigma)

    # Convenzione richiesta: g=-2*zeta=2*sigma/omega_n.
    g_value = 2.0 * sigma / omega_n
    return sigma, frequency_hz, g_value, int(np.count_nonzero(pair_valid))


def first_saturation_time(
    time: np.ndarray,
    controllers: dict[int, np.ndarray],
    start: float,
    stop: float,
) -> float:
    saturated = np.zeros(len(time), dtype=bool)
    observation = (time >= start) & (time <= stop)
    for label, limit in CONTROLLER_LIMITS_DEG.items():
        saturated |= np.abs(controllers[label]) >= 0.995 * limit
    indices = np.flatnonzero(observation & saturated)
    return math.inf if len(indices) == 0 else float(time[indices[0]])


def analyze_case(path: Path) -> dict[str, float | int | str | bool]:
    log = path.with_suffix(".log")
    velocity = log_constant(log, "V_INF", math.nan)
    trim_end = log_constant(log, "TRIM_END", 15.0)
    excitation_off = log_constant(log, "BFF_WINDOW_START", 16.75)
    observation_end = log_constant(log, "BFF_WINDOW_END", 31.75)

    with Dataset(path) as dataset:
        time = nc_array(dataset, "time")
        base_position = nc_array(dataset, f"node.struct.{BASE_NODE}.X")
        base_attitude = nc_array(dataset, f"node.struct.{BASE_NODE}.Phi")
        left_position = nc_array(dataset, f"node.struct.{LEFT_TIP}.X")
        right_position = nc_array(dataset, f"node.struct.{RIGHT_TIP}.X")
        modal_coordinates = nc_array(
            dataset, f"elem.joint.{MODAL_JOINT}.a"
        )
        if modal_coordinates.ndim == 1:
            modal_coordinates = modal_coordinates[:, np.newaxis]
        controllers = {
            label: nc_array(
                dataset, f"elem.loadable.{label}.output"
            ).squeeze() * 180.0 / math.pi
            for label in CONTROLLER_LIMITS_DEG
        }

    imported_modes = int(modal_coordinates.shape[1])
    complete = bool(len(time) and time[-1] >= observation_end)
    if imported_modes != 17 or not complete:
        return {
            "velocity_mps": velocity,
            "result": str(path),
            "imported_modes": imported_modes,
            "complete": complete,
            "accepted": False,
            "reason": "not a complete 17-mode result",
        }

    trim = (time >= trim_end - 1.0) & (time < trim_end)
    theta = base_attitude[:, 1]
    theta_trim = float(np.mean(theta[trim]))
    z_cg_trim = float(np.mean(base_position[trim, 2]))
    z_left_trim = float(np.mean(left_position[trim, 2]))
    z_right_trim = float(np.mean(right_position[trim, 2]))
    x_cg_trim = float(np.mean(base_position[trim, 0]))
    x_left = float(np.mean(left_position[trim, 0])) - x_cg_trim
    x_right = float(np.mean(right_position[trim, 0])) - x_cg_trim

    theta_relative = theta - theta_trim
    z_cg_relative = base_position[:, 2] - z_cg_trim
    rigid_left_z = z_cg_relative - x_left * theta_relative
    rigid_right_z = z_cg_relative - x_right * theta_relative
    elastic_left = left_position[:, 2] - z_left_trim - rigid_left_z
    elastic_right = right_position[:, 2] - z_right_trim - rigid_right_z
    symmetric_elastic = 0.5 * (elastic_left + elastic_right) * INCH_TO_M

    analysis_start = excitation_off + FREE_RESPONSE_DELAY_S
    saturation_time = first_saturation_time(
        time, controllers, analysis_start, observation_end
    )
    analysis_end = min(observation_end, saturation_time - 0.20)
    free = (time >= analysis_start) & (time <= analysis_end)
    if np.count_nonzero(free) < 150:
        return {
            "velocity_mps": velocity,
            "result": str(path),
            "imported_modes": imported_modes,
            "complete": complete,
            "accepted": False,
            "reason": "insufficient unsaturated free response",
        }

    free_time = time[free]
    free_pitch = theta_relative[free]
    free_bending = symmetric_elastic[free]
    center_frequency = common_frequency(
        free_time, free_pitch, free_bending
    )
    pitch_sigma, pitch_frequency, pitch_g, pitch_samples = hilbert_parameters(
        free_time, free_pitch, center_frequency
    )
    bend_sigma, bend_frequency, bend_g, bend_samples = hilbert_parameters(
        free_time, free_bending, center_frequency
    )
    finite = np.all(
        np.isfinite(
            [
                center_frequency,
                pitch_sigma,
                pitch_frequency,
                pitch_g,
                bend_sigma,
                bend_frequency,
                bend_g,
            ]
        )
    )
    frequency_agreement = abs(pitch_frequency - bend_frequency)
    accepted = bool(finite and frequency_agreement <= 0.15)
    if not finite:
        reason = "insufficient Hilbert-envelope samples"
    elif frequency_agreement > 0.15:
        reason = "pitch/bending frequency mismatch"
    else:
        reason = "accepted"
    return {
        "velocity_mps": velocity,
        "result": str(path),
        "imported_modes": imported_modes,
        "complete": complete,
        "accepted": accepted,
        "reason": reason,
        "free_start_s": analysis_start,
        "free_end_s": analysis_end,
        "first_saturation_s": saturation_time,
        "shared_spectral_frequency_hz": center_frequency,
        "pitch_sigma_1ps": pitch_sigma,
        "pitch_frequency_hz": pitch_frequency,
        "pitch_g": pitch_g,
        "pitch_hilbert_pairs": pitch_samples,
        "symmetric_bending_sigma_1ps": bend_sigma,
        "symmetric_bending_frequency_hz": bend_frequency,
        "symmetric_bending_g": bend_g,
        "symmetric_bending_hilbert_pairs": bend_samples,
        "frequency_difference_hz": frequency_agreement,
        "combined_g": 0.5 * (pitch_g + bend_g),
        "combined_sigma_1ps": 0.5 * (pitch_sigma + bend_sigma),
        "combined_frequency_hz": 0.5 * (pitch_frequency + bend_frequency),
    }


def zero_crossings(
    velocity: np.ndarray, values: np.ndarray
) -> list[float]:
    crossings: list[float] = []
    for index in range(len(velocity) - 1):
        first, second = values[index], values[index + 1]
        if not (np.isfinite(first) and np.isfinite(second)):
            continue
        if first < 0.0 <= second:
            fraction = -first / (second - first)
            crossings.append(
                float(
                    velocity[index]
                    + fraction * (velocity[index + 1] - velocity[index])
                )
            )
    return crossings


def save_plot(
    path: Path,
    velocity: np.ndarray,
    curves: tuple[tuple[np.ndarray, str, str], ...],
    ylabel: str,
    title: str,
    flutter_speed: float,
    zero_line: bool = False,
) -> None:
    plt.figure(figsize=(10, 5.5))
    for values, label, style in curves:
        plt.plot(velocity, values, style, linewidth=1.3, markersize=5, label=label)
    if zero_line:
        plt.axhline(0.0, color="black", linestyle="--", linewidth=0.9)
    if np.isfinite(flutter_speed):
        plt.axvline(
            flutter_speed,
            color="tab:red",
            linestyle="--",
            linewidth=1.2,
            label=f"Estimated flutter speed: {flutter_speed:.3f} m/s",
        )
    plt.xlabel("Equivalent airspeed [m/s]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def save_root_locus(
    path: Path,
    velocity: np.ndarray,
    pitch_sigma: np.ndarray,
    pitch_frequency: np.ndarray,
    bending_sigma: np.ndarray,
    bending_frequency: np.ndarray,
    flutter_speed: float,
) -> None:
    """Plot the response-identified complex roots in classical root-locus form."""
    pitch_real_hz = pitch_sigma / (2.0 * math.pi)
    bending_real_hz = bending_sigma / (2.0 * math.pi)

    figure, axis = plt.subplots(figsize=(10.5, 6.2))
    colour_scale = plt.Normalize(float(np.min(velocity)), float(np.max(velocity)))
    colour_map = plt.get_cmap("viridis")

    # Le linee mostrano la traiettoria; il colore dei punti identifica V.
    axis.plot(
        pitch_real_hz,
        pitch_frequency,
        color="tab:blue",
        linewidth=1.0,
        alpha=0.65,
        label="Root identified from pitch",
    )
    axis.plot(
        bending_real_hz,
        bending_frequency,
        color="tab:orange",
        linewidth=1.0,
        alpha=0.65,
        label="Root identified from symmetric elastic bending",
    )
    axis.scatter(
        pitch_real_hz,
        pitch_frequency,
        c=velocity,
        cmap=colour_map,
        norm=colour_scale,
        marker="s",
        s=38,
        edgecolors="tab:blue",
        linewidths=0.5,
        zorder=3,
    )
    points = axis.scatter(
        bending_real_hz,
        bending_frequency,
        c=velocity,
        cmap=colour_map,
        norm=colour_scale,
        marker="o",
        s=38,
        edgecolors="tab:orange",
        linewidths=0.5,
        zorder=3,
    )

    # Una freccia sul ramo combinato indica la direzione di V crescente.
    combined_real_hz = 0.5 * (pitch_real_hz + bending_real_hz)
    combined_frequency = 0.5 * (pitch_frequency + bending_frequency)
    arrow_index = max(0, len(velocity) - 3)
    if len(velocity) >= 2:
        axis.annotate(
            "",
            xy=(
                combined_real_hz[arrow_index + 1],
                combined_frequency[arrow_index + 1],
            ),
            xytext=(
                combined_real_hz[arrow_index],
                combined_frequency[arrow_index],
            ),
            arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.6},
        )
        axis.text(
            combined_real_hz[arrow_index + 1],
            combined_frequency[arrow_index + 1],
            "  Increasing airspeed",
            fontsize=9,
            va="center",
        )

    axis.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
    if np.isfinite(flutter_speed):
        flutter_index = int(np.argmin(np.abs(velocity - flutter_speed)))
        axis.annotate(
            f"Flutter boundary\n$V_f$ = {flutter_speed:.3f} m/s",
            xy=(
                combined_real_hz[flutter_index],
                combined_frequency[flutter_index],
            ),
            xytext=(18, -45),
            textcoords="offset points",
            arrowprops={"arrowstyle": "->", "color": "tab:red"},
            color="tab:red",
            fontsize=9,
        )

    colour_bar = figure.colorbar(points, ax=axis, pad=0.02)
    colour_bar.set_label("Equivalent airspeed [m/s]")
    axis.set_xlabel(r"Real part, $\mathrm{Re}(\lambda)/(2\pi)$ [Hz]")
    axis.set_ylabel(r"Imaginary part, $\mathrm{Im}(\lambda)/(2\pi)$ [Hz]")
    axis.set_title(
        "Response-identified root locus of the coupled pitch–bending mode"
    )
    axis.grid(True, alpha=0.3)
    axis.legend(loc="best")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    args = arguments()
    results = args.results.expanduser().resolve()
    output = (args.output or results / "ANALYSIS_HILBERT").resolve()
    output.mkdir(parents=True, exist_ok=True)

    paths = sorted(results.glob("V_*_mps/case.nc"))
    if not paths:
        raise FileNotFoundError(f"No simulation results found in {results}")
    rows = [analyze_case(path) for path in paths]
    rows.sort(key=lambda row: float(row["velocity_mps"]))

    summary = output / "hilbert_summary.csv"
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    accepted_rows = [row for row in rows if bool(row["accepted"])]
    if len(accepted_rows) < 2:
        raise RuntimeError("Fewer than two valid Hilbert results are available")
    velocity = np.asarray(
        [float(row["velocity_mps"]) for row in accepted_rows]
    )

    def series(name: str) -> np.ndarray:
        return np.asarray([float(row[name]) for row in accepted_rows])

    pitch_g = series("pitch_g")
    bending_g = series("symmetric_bending_g")
    combined_g = series("combined_g")
    crossings = zero_crossings(velocity, combined_g)
    flutter_speed = crossings[0] if len(crossings) == 1 else math.nan

    save_plot(
        output / "01_vg_hilbert.png",
        velocity,
        (
            (pitch_g, "Pitch response", "s-"),
            (bending_g, "Symmetric elastic bending", "o-"),
            (combined_g, "Combined coupled mode", "^-"),
        ),
        "Flutter parameter $g=2\\sigma/\\omega_n$ [-]",
        "V–g diagram from the free-response Hilbert envelope",
        flutter_speed,
        zero_line=True,
    )
    save_plot(
        output / "02_vf_hilbert.png",
        velocity,
        (
            (series("pitch_frequency_hz"), "Pitch response", "s-"),
            (
                series("symmetric_bending_frequency_hz"),
                "Symmetric elastic bending",
                "o-",
            ),
        ),
        "Damped frequency [Hz]",
        "V–f diagram of the coupled pitch–bending mode",
        flutter_speed,
    )
    save_plot(
        output / "03_sigma_hilbert.png",
        velocity,
        (
            (series("pitch_sigma_1ps"), "Pitch response", "s-"),
            (
                series("symmetric_bending_sigma_1ps"),
                "Symmetric elastic bending",
                "o-",
            ),
            (series("combined_sigma_1ps"), "Combined coupled mode", "^-"),
        ),
        "Exponential rate $\\sigma$ [1/s]",
        "Free-response growth rate from the Hilbert envelope",
        flutter_speed,
        zero_line=True,
    )
    save_root_locus(
        output / "04_root_locus_hilbert.png",
        velocity,
        series("pitch_sigma_1ps"),
        series("pitch_frequency_hz"),
        series("symmetric_bending_sigma_1ps"),
        series("symmetric_bending_frequency_hz"),
        flutter_speed,
    )

    assessment = output / "flutter_assessment_hilbert.txt"
    with assessment.open("w", encoding="utf-8") as stream:
        stream.write("HILBERT FREE-RESPONSE FLUTTER ASSESSMENT\n")
        stream.write("========================================\n\n")
        stream.write(
            "Method: band-pass filtering of the common pitch/symmetric-elastic-"
            "bending component, followed by Hilbert envelope and phase "
            "differences. No regression or R-squared criterion is used.\n\n"
        )
        if len(crossings) == 1:
            stream.write(
                f"Estimated flutter speed: {flutter_speed:.6g} m/s.\n"
            )
        elif len(crossings) > 1:
            stream.write(
                "Multiple g=0 crossings were found: "
                + ", ".join(f"{value:.6g}" for value in crossings)
                + " m/s.\n"
            )
        else:
            stream.write("No g=0 crossing was found in the accepted cases.\n")
        stream.write(
            f"Accepted simulations: {len(accepted_rows)} of {len(rows)}.\n"
        )

    print(f"Summary: {summary}")
    print(f"V-g diagram: {output / '01_vg_hilbert.png'}")
    print(f"V-f diagram: {output / '02_vf_hilbert.png'}")
    print(f"Sigma diagram: {output / '03_sigma_hilbert.png'}")
    print(f"Root locus: {output / '04_root_locus_hilbert.png'}")
    print(f"Assessment: {assessment}")


if __name__ == "__main__":
    main()
