"""Build V-g and V-f diagrams from the free response of the X-56 model.

The flap command in INCLUDE/control_surface_joints_flutter.mbd is

    sine, 1.0, 12.566, +/-0.12217, 2.5, 0.0

so the 2 Hz, 2.5-cycle command finishes at approximately 2.25 s.  Samples
before that instant are deliberately excluded: damping is identified only
from the subsequent free response.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.signal import butter, detrend, find_peaks, sosfiltfilt
from scipy.stats import theilslopes


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output"

# Parameters of the flap drive in control_surface_joints_flutter.mbd.
FLAP_START = 1.0
FLAP_ANGULAR_FREQUENCY = 12.566  # rad/s
FLAP_CYCLES = 2.5
FLAP_END = FLAP_START + FLAP_CYCLES * 2.0 * np.pi / FLAP_ANGULAR_FREQUENCY

LEFT_TIP = "990020"
RIGHT_TIP = "991020"

# The retained FEM modes start at structural mode 7.  Its first four vacuum
# frequencies are about 3.22, 5.30, 8.71 and 11.16 Hz.  Mode 7 is not reported:
# this flap input excites it too weakly, so its band is dominated by a residual
# at the 2 Hz flap-drive frequency.  Separate bands prevent the other damping
# estimates from mixing two modes in one logarithmic decrement.
MODE_BANDS = {
    "Mode 8 (antisymmetric tip Z)": (4.2, 6.8, "antisymmetric"),
    "Mode 9 (symmetric tip Z)": (6.8, 9.7, "symmetric"),
    "Mode 10 (symmetric tip Z)": (9.7, 12.5, "symmetric"),
}

RUN_PATTERN = re.compile(r"out_V(?P<velocity>\d+(?:\.\d+)?)\.nc$")
MIN_PEAKS = 6
MIN_R_SQUARED = 0.75
FILTER_ORDER = 4
EDGE_TRIM = 0.20  # remove the short filtfilt boundary transient [s]


@dataclass
class ModalEstimate:
    velocity: float
    mode: str
    damping_g: float
    frequency_hz: float
    growth_rate: float
    r_squared: float
    peaks_used: int
    peak_amplitude: float
    reliable: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate modal damping from the post-flap free response."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="directory containing out_V*.nc files (default: V-g/output)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show the figure interactively after saving it",
    )
    return parser.parse_args()


def discover_runs(output_dir: Path) -> list[tuple[float, Path]]:
    """Return all available simulations, sorted by numeric velocity."""
    runs: list[tuple[float, Path]] = []
    for path in output_dir.glob("out_V*.nc"):
        match = RUN_PATTERN.match(path.name)
        if match:
            runs.append((float(match.group("velocity")), path))
    return sorted(runs, key=lambda item: item[0])


def bandpass(signal: np.ndarray, band: tuple[float, float], fs: float) -> np.ndarray:
    """Zero-phase Butterworth filtering using numerically stable SOS form."""
    low, high = band
    nyquist = 0.5 * fs
    if not (0.0 < low < high < nyquist):
        raise ValueError(f"invalid filter band {band} Hz for fs={fs:.3f} Hz")
    sos = butter(
        FILTER_ORDER,
        [low / nyquist, high / nyquist],
        btype="band",
        output="sos",
    )
    return sosfiltfilt(sos, signal)


def modal_damping(
    time: np.ndarray,
    signal: np.ndarray,
    band: tuple[float, float],
) -> tuple[float, float, float, float, int, float]:
    """Estimate g=2*zeta and frequency from a band-limited free decay.

    A robust line is fitted to log(peak amplitude) versus time.  Its slope is
    the modal growth rate sigma: sigma < 0 is stable and therefore gives
    positive damping; sigma > 0 gives negative damping (flutter convention).
    """
    time = np.asarray(time, dtype=float)
    signal = np.asarray(signal, dtype=float)
    finite = np.isfinite(time) & np.isfinite(signal)
    time, signal = time[finite], signal[finite]
    if time.size < 50:
        return (np.nan,) * 4 + (0, np.nan)

    dt = float(np.median(np.diff(time)))
    fs = 1.0 / dt
    response = detrend(signal, type="linear")
    filtered = bandpass(response, band, fs)

    usable = (time >= time[0] + EDGE_TRIM) & (time <= time[-1] - EDGE_TRIM)
    t_fit, s_fit = time[usable], filtered[usable]
    if t_fit.size < 50:
        return (np.nan,) * 4 + (0, np.nan)

    # Positive peaks give one timestamp per complete oscillation.  The minimum
    # distance rejects small ripples without imposing the expected frequency.
    min_distance = max(1, int(0.70 * fs / band[1]))
    prominence = max(np.ptp(s_fit) * 1.0e-5, np.finfo(float).eps)
    peaks, _ = find_peaks(s_fit, distance=min_distance, prominence=prominence)
    if peaks.size < MIN_PEAKS:
        return (np.nan,) * 4 + (int(peaks.size), float(np.max(np.abs(s_fit))))

    peak_times = t_fit[peaks]
    amplitudes = np.abs(s_fit[peaks])

    # Once a stable response reaches its numerical floor, late peaks bias the
    # slope toward zero.  Keep the coherent portion of a decaying response;
    # growing responses retain all peaks so negative damping is measurable.
    first_level = float(np.median(amplitudes[: min(3, amplitudes.size)]))
    last_level = float(np.median(amplitudes[-min(3, amplitudes.size) :]))
    if last_level < first_level:
        tail = filtered[-max(20, filtered.size // 8) :]
        noise_level = 1.4826 * float(np.median(np.abs(tail - np.median(tail))))
        threshold = max(2.5 * noise_level, 0.01 * float(np.max(amplitudes)))
        coherent = amplitudes >= threshold
        # Use the first continuous coherent block; isolated late noise peaks
        # are not part of the free decay.
        below = np.flatnonzero(~coherent)
        stop = int(below[0]) if below.size else amplitudes.size
        if stop >= MIN_PEAKS:
            peak_times = peak_times[:stop]
            amplitudes = amplitudes[:stop]

    if amplitudes.size < MIN_PEAKS or np.any(amplitudes <= 0.0):
        return (np.nan,) * 4 + (int(amplitudes.size), float(np.max(amplitudes)))

    log_amplitudes = np.log(amplitudes)
    growth_rate = float(theilslopes(log_amplitudes, peak_times).slope)
    intercept = float(np.median(log_amplitudes - growth_rate * peak_times))
    fitted = growth_rate * peak_times + intercept
    residual = float(np.sum((log_amplitudes - fitted) ** 2))
    total = float(np.sum((log_amplitudes - np.mean(log_amplitudes)) ** 2))
    r_squared = 1.0 - residual / total if total > 0.0 else np.nan

    # Robustly fit time versus cycle number.  This uses every peak and avoids
    # quantising the reported frequency to an integer number of time steps.
    cycle_number = np.arange(peak_times.size, dtype=float)
    period = float(theilslopes(peak_times, cycle_number).slope)
    if period <= 0.0:
        return (np.nan,) * 4 + (int(amplitudes.size), float(np.max(amplitudes)))
    frequency = 1.0 / period
    omega_d = 2.0 * np.pi * frequency
    zeta = -growth_rate / np.sqrt(omega_d**2 + growth_rate**2)
    damping_g = 2.0 * zeta
    return (
        float(damping_g),
        float(frequency),
        growth_rate,
        float(r_squared),
        int(amplitudes.size),
        float(np.max(amplitudes)),
    )


def analyse_run(velocity: float, path: Path) -> list[ModalEstimate]:
    """Read one NetCDF file and identify the configured flexible modes."""
    estimates: list[ModalEstimate] = []
    with Dataset(path, "r") as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        left_z = np.asarray(
            dataset.variables[f"node.struct.{LEFT_TIP}.X"][:, 2], dtype=float
        )
        right_z = np.asarray(
            dataset.variables[f"node.struct.{RIGHT_TIP}.X"][:, 2], dtype=float
        )

    response_mask = time >= FLAP_END
    response_time = time[response_mask]
    symmetric_z = 0.5 * (left_z[response_mask] + right_z[response_mask])
    antisymmetric_z = 0.5 * (left_z[response_mask] - right_z[response_mask])
    channels = {"symmetric": symmetric_z, "antisymmetric": antisymmetric_z}

    for mode, (low, high, channel) in MODE_BANDS.items():
        g, frequency, sigma, r2, peaks, amplitude = modal_damping(
            response_time, channels[channel], (low, high)
        )
        estimates.append(
            ModalEstimate(
                velocity=velocity,
                mode=mode,
                damping_g=g,
                frequency_hz=frequency,
                growth_rate=sigma,
                r_squared=r2,
                peaks_used=peaks,
                peak_amplitude=amplitude,
                reliable=bool(
                    np.isfinite(g)
                    and np.isfinite(frequency)
                    and np.isfinite(r2)
                    and peaks >= MIN_PEAKS
                    and r2 >= MIN_R_SQUARED
                ),
            )
        )
    return estimates


def save_csv(estimates: list[ModalEstimate], path: Path) -> None:
    fields = list(ModalEstimate.__dataclass_fields__)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for estimate in estimates:
            writer.writerow(estimate.__dict__)


def zero_crossings(velocity: np.ndarray, damping: np.ndarray) -> list[float]:
    """Linearly interpolate finite V-g zero crossings."""
    crossings: list[float] = []
    adjacent_points = zip(
        velocity[:-1], velocity[1:], damping[:-1], damping[1:]
    )
    for v0, v1, g0, g1 in adjacent_points:
        if not np.all(np.isfinite([v0, v1, g0, g1])):
            continue
        if g0 == 0.0:
            crossings.append(float(v0))
        elif g0 * g1 < 0.0:
            crossings.append(float(v0 - g0 * (v1 - v0) / (g1 - g0)))
    return crossings


def plot_results(estimates: list[ModalEstimate], output_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    colors = plt.get_cmap("tab10").colors

    for color, mode in zip(colors, MODE_BANDS):
        rows = [item for item in estimates if item.mode == mode]
        velocity = np.array([item.velocity for item in rows])
        damping = np.array([item.damping_g for item in rows])
        frequency = np.array([item.frequency_hz for item in rows])
        reliable = np.array([item.reliable for item in rows])
        valid = np.isfinite(damping) & np.isfinite(frequency) & reliable
        plot_damping = np.where(valid, damping, np.nan)
        plot_frequency = np.where(valid, frequency, np.nan)
        axes[0].plot(
            velocity, plot_damping, "o-", ms=4, color=color, label=mode
        )
        axes[1].plot(
            velocity, plot_frequency, "s-", ms=4, color=color, label=mode
        )

        for crossing in zero_crossings(velocity, plot_damping):
            axes[0].axvline(crossing, color=color, linestyle=":", alpha=0.35)

    axes[0].axhline(0.0, color="black", linestyle="--", linewidth=1.0)
    axes[0].set_title(f"V-g diagram: free response from t = {FLAP_END:.2f} s")
    axes[0].set_ylabel("Damping g = 2ζ  (negative = unstable)")
    axes[1].set_title("V-f diagram: identified flexible modes")
    axes[1].set_ylabel("Damped frequency [Hz]")
    for axis in axes:
        axis.set_xlabel("Velocity [m/s]")
        axis.grid(True, alpha=0.35)
        axis.legend(fontsize=8)

    fig.tight_layout()
    figure_path = output_dir / "Vg_Vf_free_response.png"
    fig.savefig(figure_path, dpi=180)
    return figure_path


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    runs = discover_runs(output_dir)
    if not runs:
        raise SystemExit(f"No out_V*.nc simulations found in {output_dir}")

    print(
        f"Found {len(runs)} simulations: {runs[0][0]:g} to {runs[-1][0]:g} m/s."
    )
    print(
        f"Flap command: {FLAP_START:.2f} to {FLAP_END:.3f} s; "
        f"only t >= {FLAP_END:.3f} s is analysed."
    )

    estimates: list[ModalEstimate] = []
    for velocity, path in runs:
        try:
            estimates.extend(analyse_run(velocity, path))
        except (KeyError, OSError, ValueError) as error:
            print(f"Skipping V={velocity:g} m/s ({path.name}): {error}")

    csv_path = output_dir / "Vg_free_response.csv"
    save_csv(estimates, csv_path)
    figure_path = plot_results(estimates, output_dir)

    print(f"Saved {csv_path}")
    print(f"Saved {figure_path}")
    for mode in MODE_BANDS:
        rows = [item for item in estimates if item.mode == mode]
        reliable = np.array([item.reliable for item in rows])
        velocity = np.array([item.velocity for item in rows])
        damping = np.array(
            [item.damping_g if item.reliable else np.nan for item in rows]
        )
        crossings = zero_crossings(velocity, damping)
        print(f"{mode}: {int(reliable.sum())}/{len(rows)} reliable fits")
        if crossings:
            values = ", ".join(f"{value:.2f}" for value in crossings)
            print(f"{mode}: estimated g=0 crossing(s) at {values} m/s")

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
