#!/usr/bin/env python3
"""Visualizza trim, volo e accoppiamento short-period/1SWB di un run MBDyn."""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_bbf_auto_trim_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

INCH_TO_M = 0.0254
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167
RAD_TO_DEG = 180.0 / math.pi
BASE_NODE = 990000
TRIM_JOINT = 23
MODAL_JOINT = 5
FLAPS = {
    1004: "BFL", 1008: "WFL1", 1011: "WFL2", 1014: "WFL3", 1017: "WFL4",
    2004: "BFR", 2008: "WFR1", 2011: "WFR2", 2014: "WFR3", 2017: "WFR4",
}
PIDS = {
    9101: "trim pitch", 9102: "trim elevator", 9301: "Vz/altitude",
    9303: "pitch", 9305: "pitch-rate", 9307: "heave integral",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot one BBF_AUTO_TRIM NetCDF result")
    parser.add_argument("result", nargs="?", default="output/case.nc")
    return parser.parse_args()


def source_constant(source: str, name: str, fallback: float) -> float:
    """Estrae una costante gia' valutata dal testo MBDyn."""
    match = re.search(
        rf"(?:set\s*:\s*)?const\s+real\s+{re.escape(name)}\s*=\s*"
        # Il sorgente termina con ';', mentre il symbol table del .log no.
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*;?",
        source,
        flags=re.IGNORECASE,
    )
    return fallback if match is None else float(match.group(1).replace("D", "E"))


def result_constant(result: Path, name: str, fallback: float) -> float:
    """Usa il log del run: un risultato conserva cosi' velocita' e tempi propri."""
    log_path = result.with_suffix(".log")
    if log_path.is_file():
        value = source_constant(log_path.read_text(encoding="utf-8"), name, math.nan)
        if np.isfinite(value):
            return value
    main = (Path(__file__).parent / "main_bbf.mbd").read_text(encoding="utf-8")
    return source_constant(main, name, fallback)


def data(nc: Dataset, name: str) -> np.ndarray:
    if name not in nc.variables:
        raise KeyError(f"Variabile NetCDF MBDyn assente: {name}")
    return np.ma.filled(nc.variables[name][:], np.nan).astype(float)


def detrend(time: np.ndarray, signal: np.ndarray) -> np.ndarray:
    good = np.isfinite(signal)
    if np.count_nonzero(good) < 2:
        return signal - np.nanmean(signal)
    return signal - np.polyval(np.polyfit(time[good], signal[good], 1), time)


def remove_slow_component(time: np.ndarray, signal: np.ndarray, seconds: float = 0.8):
    """Separa il 1SWB dalla deformazione quasi-statica e dal moto rigido."""
    count = max(3, int(round(seconds / np.median(np.diff(time)))))
    if count % 2 == 0:
        count += 1
    pad = count // 2
    trend = np.convolve(
        np.pad(signal, pad, mode="edge"), np.ones(count) / count, mode="valid"
    )
    return signal - trend


def spectrum(time: np.ndarray, signal: np.ndarray, fmin: float = 0.1):
    signal = detrend(time, signal)
    window = np.hanning(len(signal))
    amplitude = 2.0 * np.abs(np.fft.rfft(signal * window)) / max(np.sum(window), 1.0)
    frequency = np.fft.rfftfreq(len(signal), np.median(np.diff(time)))
    valid = frequency >= fmin
    peak = frequency[valid][np.argmax(amplitude[valid])] if np.any(valid) else math.nan
    return frequency, amplitude, float(peak)


def growth(time: np.ndarray, signal: np.ndarray, frequency: float):
    """Stima il tasso esponenziale dai massimi locali dell'ampiezza detrendizzata."""
    amplitude = np.abs(detrend(time, signal))
    candidates = np.flatnonzero(
        (amplitude[1:-1] > amplitude[:-2]) & (amplitude[1:-1] >= amplitude[2:])
    ) + 1
    dt = float(np.median(np.diff(time)))
    spacing = max(1, int(0.60 / (max(frequency, 0.5) * dt)))
    peaks: list[int] = []
    for index in candidates:
        if not peaks or index - peaks[-1] >= spacing:
            peaks.append(int(index))
        elif amplitude[index] > amplitude[peaks[-1]]:
            peaks[-1] = int(index)
    indices = np.asarray(peaks, dtype=int)
    floor = max(np.nanmax(amplitude) * 1.e-5, np.finfo(float).eps)
    indices = indices[amplitude[indices] > floor]
    if len(indices) < 3:
        return math.nan, indices, np.full_like(amplitude, np.nan)
    rate, intercept = np.polyfit(time[indices], np.log(amplitude[indices]), 1)
    return float(rate), indices, np.exp(intercept + rate * time)


def damping_ratio(rate: float, frequency: float) -> float:
    if not np.isfinite(rate) or not np.isfinite(frequency):
        return math.nan
    omega = 2.0 * math.pi * frequency
    return -rate / math.sqrt(omega * omega + rate * rate)


def save(fig: plt.Figure, destination: Path) -> None:
    fig.tight_layout()
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def run() -> None:
    path = Path(arguments().result)
    if path.suffix != ".nc":
        path = path.with_suffix(".nc")
    if not path.is_file():
        raise FileNotFoundError(path)
    figures = path.parent / f"{path.stem}_plots"
    figures.mkdir(parents=True, exist_ok=True)

    velocity = result_constant(path, "V_INF", math.nan)
    trim_end = result_constant(path, "TRIM_END", 12.0)
    excitation_start = result_constant(path, "EXCITATION_START", trim_end + 6.0)
    # Compatibilita' con i primi run di sviluppo che usavano OPEN_LOOP_*.
    old_start = result_constant(path, "OPEN_LOOP_START", excitation_start + 0.25)
    old_end = result_constant(path, "OPEN_LOOP_END", old_start + 1.50)
    bff_start = result_constant(path, "BFF_WINDOW_START", old_start)
    bff_end = result_constant(path, "BFF_WINDOW_END", old_end)

    with Dataset(path) as nc:
        time = data(nc, "time")
        position = data(nc, f"node.struct.{BASE_NODE}.X")
        attitude = data(nc, f"node.struct.{BASE_NODE}.Phi")
        velocity_xyz = data(nc, f"node.struct.{BASE_NODE}.XP")
        omega = data(nc, f"node.struct.{BASE_NODE}.Omega")
        modal_q = data(nc, f"elem.joint.{MODAL_JOINT}.a")[:, 0]
        trim_force = data(nc, f"elem.joint.{TRIM_JOINT}.f")[:, 2] * LBF_TO_N
        trim_moment = data(nc, f"elem.joint.{TRIM_JOINT}.m")[:, 1] * LBFIN_TO_NM
        flap = {
            label: data(nc, f"elem.joint.{label}.Phi")[:, 1] * RAD_TO_DEG
            for label in FLAPS
        }
        pid = {
            label: data(nc, f"elem.loadable.{label}.output").squeeze() * RAD_TO_DEG
            for label in PIDS
        }

    heave = (position[:, 2] - position[0, 2]) * INCH_TO_M
    vertical_speed = velocity_xyz[:, 2] * INCH_TO_M
    pitch = attitude[:, 1] * RAD_TO_DEG
    pitch_rate = omega[:, 1] * RAD_TO_DEG
    # AoA rigido: assetto meno angolo di traiettoria verticale del CG.
    aoa = pitch - np.degrees(np.arctan2(vertical_speed, velocity))
    xy_drift = np.linalg.norm((position[:, :2] - position[0, :2]) * INCH_TO_M, axis=1)
    roll_yaw_drift = np.linalg.norm(attitude[:, [0, 2]] * RAD_TO_DEG, axis=1)
    trim_window = (time >= trim_end - 1.0) & (time < trim_end)
    settle_window = (time >= trim_end + 0.5) & (time < excitation_start)
    bff_window = (time >= bff_start + 0.02) & (time <= bff_end)
    recovery_window = time >= max(bff_end, time[-1] - 2.0)

    if np.count_nonzero(bff_window) < 20:
        raise RuntimeError("Finestra BFF troppo corta o assente nel risultato")

    # Figura 1: convergenza del trim calcolato nello stesso run.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axes[0, 0].plot(time, trim_force)
    axes[0, 0].set_ylabel("Trim Fz [N]")
    axes[0, 1].plot(time, trim_moment)
    axes[0, 1].set_ylabel("Trim My [N m]")
    axes[1, 0].plot(time, pid[9101])
    axes[1, 0].set_ylabel("Pitch di trim [deg]")
    axes[1, 1].plot(time, pid[9102])
    axes[1, 1].set_ylabel("Elevator di trim [deg]")
    for axis in axes.flat:
        axis.axvline(trim_end, color="k", ls="--", lw=0.8)
        axis.grid(True, alpha=0.3)
        axis.set_xlabel("Tempo [s]")
    fig.suptitle(f"Trim automatico nello stesso run — V={velocity:g} m/s")
    save(fig, figures / "trim.png")

    # Figura 2: il controllo deve mantenere quota e assetto prima e dopo il test.
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for axis, signal, label in zip(
        axes.flat,
        (heave, vertical_speed, pitch, aoa, pitch_rate, modal_q),
        (
            "Heave rigido [m]", "Vz [m/s]", "Pitch [deg]",
            "AoA rigido [deg]", "Pitch-rate [deg/s]", "q1 / 1SWB [-]",
        ),
    ):
        axis.plot(time, signal)
        axis.axvline(trim_end, color="k", ls="--", lw=0.8)
        axis.axvspan(excitation_start, bff_start, color="tab:orange", alpha=0.25)
        axis.axvspan(bff_start, bff_end, color="tab:red", alpha=0.12)
        axis.set_ylabel(label)
        axis.set_xlabel("Tempo [s]")
        axis.grid(True, alpha=0.3)
    fig.suptitle("Volo longitudinale: doublet arancione, osservazione BFF rossa")
    save(fig, figures / "flight.png")

    # Figura 3: confronto diretto short-period/1SWB con il SAS di volo attivo.
    test_time = time[bff_window]
    pitch_test = detrend(test_time, pitch[bff_window])
    bending_test = remove_slow_component(test_time, modal_q[bff_window])
    fp, ap, pitch_frequency = spectrum(test_time, pitch_test)
    # Il primo bending e' cercato sopra 1 Hz, fuori da heave e short-period.
    fb, ab, bending_frequency = spectrum(test_time, bending_test, fmin=1.0)
    pitch_growth, pitch_peaks, pitch_fit = growth(test_time, pitch_test, pitch_frequency)
    bending_growth, bending_peaks, bending_fit = growth(test_time, bending_test, bending_frequency)

    # Oltre questi limiti la risposta non e' piu' una perturbazione del trim:
    # frequenze e damping lineari non hanno significato e vengono marcati NaN.
    trim_pitch = float(np.mean(pid[9101][trim_window]))
    linear_response = bool(
        np.nanmax(np.abs(pitch[bff_window] - trim_pitch)) <= 10.0
        and np.nanmax(np.abs(heave[bff_window])) <= 5.0
    )
    if not linear_response:
        pitch_frequency = pitch_growth = math.nan
        bending_frequency = bending_growth = math.nan
        pitch_fit[:] = np.nan
        bending_fit[:] = np.nan

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].plot(test_time, pitch_test / max(np.std(pitch_test), 1.e-12), label="pitch")
    axes[0, 0].plot(test_time, bending_test / max(np.std(bending_test), 1.e-12), label="q1")
    axes[0, 0].set_ylabel("Ampiezza normalizzata")
    axes[0, 0].legend()
    axes[0, 1].plot(fp, ap / max(np.max(ap), 1.e-12), label="pitch")
    axes[0, 1].plot(fb, ab / max(np.max(ab), 1.e-12), label="q1")
    axes[0, 1].set_xlim(0.0, 8.0)
    axes[0, 1].set_ylabel("FFT normalizzata")
    axes[0, 1].legend()
    axes[1, 0].plot(test_time, np.abs(pitch_test), color="0.6")
    axes[1, 0].plot(test_time[pitch_peaks], np.abs(pitch_test)[pitch_peaks], "o")
    axes[1, 0].plot(test_time, pitch_fit, "r--")
    axes[1, 0].set_ylabel("Inviluppo pitch [deg]")
    axes[1, 1].plot(test_time, np.abs(bending_test), color="0.6")
    axes[1, 1].plot(test_time[bending_peaks], np.abs(bending_test)[bending_peaks], "o")
    axes[1, 1].plot(test_time, bending_fit, "r--")
    axes[1, 1].set_ylabel("Inviluppo q1 [-]")
    for axis in axes.flat:
        axis.set_xlabel("Tempo [s]" if axis in axes[1] else "Frequenza [Hz]" if axis is axes[0, 1] else "Tempo [s]")
        axis.grid(True, alpha=0.3)
    fig.suptitle("Interazione short-period / primo bending simmetrico")
    save(fig, figures / "bff.png")

    # Figura 4: WF4 riceve il doublet; WF1 e body flap mantengono il volo.
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for label, name in FLAPS.items():
        axes[0 if label < 2000 else 1].plot(time, flap[label], label=name)
    for axis in axes:
        axis.axvspan(excitation_start, bff_start, color="tab:orange", alpha=0.25)
        axis.axvspan(bff_start, bff_end, color="tab:red", alpha=0.12)
        axis.set_ylabel("Deflessione [deg]")
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=5, fontsize=8)
    axes[1].set_xlabel("Tempo [s]")
    fig.suptitle("Superfici di controllo")
    save(fig, figures / "controls.png")

    report = [
        f"Result: {path}",
        f"Velocity: {velocity:.6g} m/s",
        f"Automatic trim pitch/elevator: {np.mean(pid[9101][trim_window]):.6g} / "
        f"{np.mean(pid[9102][trim_window]):.6g} deg",
        f"Trim residual Fz mean/std: {np.mean(trim_force[trim_window]):.6g} / "
        f"{np.std(trim_force[trim_window]):.6g} N",
        f"Trim residual My mean/std: {np.mean(trim_moment[trim_window]):.6g} / "
        f"{np.std(trim_moment[trim_window]):.6g} N m",
        f"Pre-test mean rigid heave: {np.mean(heave[settle_window]):.6g} m",
        f"Recovery mean rigid heave: {np.mean(heave[recovery_window]):.6g} m",
        f"BFF-window rigid AoA range: {np.ptp(aoa[bff_window]):.6g} deg",
        f"BFF-window rigid heave range: {np.ptp(heave[bff_window]):.6g} m",
        f"Maximum constrained X/Y drift: {np.nanmax(xy_drift):.6g} m",
        f"Maximum constrained roll/yaw drift: {np.nanmax(roll_yaw_drift):.6g} deg",
        f"Linear-range response: {'yes' if linear_response else 'no; damping estimate invalid'}",
        f"BFF-window pitch frequency/rate/damping: {pitch_frequency:.6g} Hz / "
        f"{pitch_growth:.6g} 1/s / {damping_ratio(pitch_growth, pitch_frequency):.6g}",
        f"BFF-window 1SWB frequency/rate/damping: {bending_frequency:.6g} Hz / "
        f"{bending_growth:.6g} 1/s / {damping_ratio(bending_growth, bending_frequency):.6g}",
        f"Pitch-1SWB frequency separation: {abs(pitch_frequency-bending_frequency):.6g} Hz",
        f"Figures: {figures}",
    ]
    print("\n".join(report))


if __name__ == "__main__":
    run()
