#!/usr/bin/env python3
"""Postprocess true open-loop X-56 longitudinal free responses."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bff_open_loop_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.optimize import least_squares
from scipy.signal import butter, detrend, sosfiltfilt
from scipy.spatial.transform import Rotation

from modal_identification import DETRENDS, WINDOWS_S, detrend_variant, identify_multimodal
from nastran_flutter_reference import DEFAULT_F06, interpolate_at_velocity, parse_flutter_point, zero_crossing

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path("/mnt/c/Users/Utente/Desktop/BFF_open_loop")
IN_TO_M = 0.0254
RAD_TO_DEG = 180.0 / math.pi
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167
WEIGHT_LBF = 419.63493154586183
S_REF_IN2 = 8064.0
NASTRAN_CZ_ALPHA_PER_DEG = 0.1083271
GAUSS_W = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
PID_LABELS = (9101, 9102, 9103, 9104, 9201, 9202, 9301, 9302, 9303)
PID_NAMES = ("altitude", "vertical_velocity", "pitch", "pitch_rate", "roll", "roll_rate", "yaw", "yaw_rate", "lateral_velocity")
SURFACE_NODES = {
    "BFL": (990004, 880004), "BFR": (991004, 881004),
    "WF1L": (990008, 880008), "WF1R": (991008, 881008),
    "WF2L": (990011, 880011), "WF2R": (991011, 881011),
    "WF3L": (990014, 880014), "WF3R": (991014, 881014),
    "WF4L": (990017, 880017), "WF4R": (991017, 881017),
}


def arr(data: Dataset, name: str) -> np.ndarray:
    return np.asarray(data[name][:]).squeeze()


def model_constants(text: str) -> dict[str, float]:
    values = {"pi": math.pi, "deg2rad": math.pi / 180.0, "m2in": 1.0 / IN_TO_M}
    for match in re.finditer(r"(?m)^\s*set:\s*const\s+(?:real|integer)\s+(\w+)\s*=\s*([^;]+);", text):
        name, expression = match.groups()
        expression = expression.split("#", 1)[0].strip()
        try:
            values[name] = float(eval(expression, {"__builtins__": {}, "cos": math.cos, "sin": math.sin}, values))
        except (NameError, SyntaxError, TypeError, ValueError, ZeroDivisionError):
            pass
    return values


def filter_first_order(signal: np.ndarray, a1: float, b0: float, b1: float) -> np.ndarray:
    y = np.zeros_like(signal, dtype=float)
    for i, value in enumerate(signal):
        y[i] = a1 * (y[i - 1] if i else 0.0) + b0 * value + b1 * (signal[i - 1] if i else 0.0)
    return y


def sample_and_hold(live: np.ndarray, trigger: np.ndarray, initial: float) -> np.ndarray:
    """Replay MBDyn SHDriveCaller: update after convergence of a triggered step."""
    held = np.empty_like(live)
    stored = initial
    for i in range(len(live)):
        held[i] = stored
        if trigger[i]:
            stored = live[i]
    return held


def detrended_harmonic_amplitude(time: np.ndarray, signal: np.ndarray, frequency: float) -> float:
    """Least-squares sinusoidal amplitude after removing offset and linear drift."""
    t = np.asarray(time, float) - float(time[0])
    y = np.asarray(signal, float)
    design = np.column_stack((
        np.ones_like(t), t,
        np.sin(2.0 * math.pi * frequency * t),
        np.cos(2.0 * math.pi * frequency * t),
    ))
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    return float(math.hypot(coefficients[2], coefficients[3]))


def surface_kinematics(data: Dataset) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for name, (parent, control) in SURFACE_NODES.items():
        r_parent = Rotation.from_rotvec(arr(data, f"node.struct.{parent}.Phi"))
        r_control = Rotation.from_rotvec(arr(data, f"node.struct.{control}.Phi"))
        result[name] = np.linalg.norm((r_parent.inv() * r_control).as_rotvec(), axis=1) * RAD_TO_DEG
    return result


@lru_cache(maxsize=1)
def aerodynamic_spans() -> dict[int, float]:
    text = (ROOT / "INCLUDE/aerobody.mbd").read_text()
    spans: dict[int, float] = {}
    for block in re.split(r"(?=\s*aerodynamic body:)", text):
        match = re.search(r"aerodynamic body:\s*(\d+)", block)
        if not match:
            continue
        lines = [line.strip() for line in block.splitlines()]
        references = [i for i, line in enumerate(lines) if line.startswith("reference, node")]
        expression = lines[references[1] + 1].split(",", 1)[0]
        spans[int(match.group(1))] = float(eval(expression, {"__builtins__": {}}, {}))
    if len(spans) != 58:
        raise RuntimeError(f"attesi 58 elementi aerodinamici, trovati {len(spans)}")
    return spans


def point_force_moment_residual(data: Dataset, index: int, cg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    force = np.zeros(3)
    moment = np.zeros(3)
    for label, span in aerodynamic_spans().items():
        for point, weight in enumerate(GAUSS_W):
            scale = 0.5 * span * weight
            f = np.asarray(data[f"elem.aerodynamic.{label}.F_{point}"][index]) * scale
            m = np.asarray(data[f"elem.aerodynamic.{label}.M_{point}"][index]) * scale
            x = np.asarray(data[f"elem.aerodynamic.{label}.X_{point}"][index])
            force += f
            moment += m + np.cross(x - cg, f)
    force[2] -= WEIGHT_LBF
    return force * LBF_TO_N, moment * LBFIN_TO_NM


@lru_cache(maxsize=1)
def identify_swb1_from_fem() -> dict[str, object]:
    """Find the lowest-frequency symmetric vertical wing mode, not a fixed index."""
    text = (ROOT / "INCLUDE/mbdyn_modal.fem").read_text(errors="ignore")

    def section(start: str, end: str) -> str:
        return text.split(start, 1)[1].split(end, 1)[0]

    nodes = [int(x) for x in re.findall(r"[-+]?\d+", section("** RECORD GROUP 2, FINITE ELEMENT NODE LIST", "** RECORD GROUP 3"))]
    index = {node: i for i, node in enumerate(nodes)}
    blocks = re.split(r"\*\*\s+NORMAL MODE SHAPE #\s*(\d+)", section("** RECORD GROUP 8, MODE SHAPES", "** RECORD GROUP 9"))
    shapes: dict[int, np.ndarray] = {}
    for i in range(1, len(blocks), 2):
        mode = int(blocks[i])
        if 7 <= mode <= 12:
            shapes[mode] = np.fromstring(blocks[i + 1].replace("**", ""), sep=" ").reshape(len(nodes), 6)
    pairs = [(990000 + i, 991000 + i) for i in range(2, 24) if 990000 + i in index and 991000 + i in index]
    stiffness = np.fromstring(section("** RECORD GROUP 10, MODAL STIFFNESS MATRIX", "** RECORD GROUP 11").replace("**", ""), sep=" ").reshape(60, 60)
    candidates = []
    for mode, shape in shapes.items():
        left = np.array([shape[index[a], 2] for a, _ in pairs])
        right = np.array([shape[index[b], 2] for _, b in pairs])
        denom = np.linalg.norm(left) + np.linalg.norm(right) + 1e-30
        symmetric_residual = float(np.linalg.norm(left - right) / denom)
        frequency = float(math.sqrt(max(stiffness[mode - 1, mode - 1], 0.0)) / (2.0 * math.pi))
        candidates.append((frequency, mode, symmetric_residual))
    symmetric = [item for item in candidates if item[2] < 0.20]
    if not symmetric:
        raise RuntimeError("nessun modo wing-bending simmetrico trovato nella base FEM 7--12")
    frequency, mode, residual = min(symmetric)
    return {"fem_mode": mode, "modal_column": mode - 7, "dry_frequency_hz": frequency, "symmetry_residual": residual, "candidates": candidates}


def fit_damped_sine(time: np.ndarray, signal: np.ndarray, low: float = 0.7, high: float = 4.0) -> dict[str, float | bool | None]:
    t = np.asarray(time, float) - float(time[0])
    y = np.asarray(signal, float)
    if len(t) < 80 or not np.all(np.isfinite(y)) or np.ptp(y) < 1e-10:
        return {"frequency_hz": None, "sigma_per_s": None, "damping_ratio": None, "phase_deg": None, "r2": 0.0, "reliable": False}
    base = np.polyval(np.polyfit(t, y, 1), t)
    yd = y - base
    dt = float(np.median(np.diff(t)))
    frequencies = np.fft.rfftfreq(len(t), dt)
    spectrum = np.abs(np.fft.rfft(yd * np.hanning(len(t))))
    use = (frequencies >= low) & (frequencies <= high)
    f0 = float(frequencies[use][np.argmax(spectrum[use])])
    amp = max(float(np.std(yd) * math.sqrt(2.0)), 1e-10)
    x0 = np.array([base[0], (base[-1] - base[0]) / max(t[-1], dt), amp, 0.0, 0.0, f0])
    scale = max(float(np.std(yd)), 1e-10)

    def model(x: np.ndarray) -> np.ndarray:
        c0, c1, a, b, sigma, frequency = x
        return c0 + c1 * t + np.exp(np.clip(sigma * t, -50.0, 50.0)) * (
            a * np.sin(2.0 * np.pi * frequency * t) + b * np.cos(2.0 * np.pi * frequency * t)
        )

    fit = least_squares(
        lambda x: (model(x) - y) / scale,
        x0,
        bounds=([-np.inf, -np.inf, -np.inf, -np.inf, -8.0, low], [np.inf, np.inf, np.inf, np.inf, 8.0, high]),
        loss="soft_l1",
        max_nfev=5000,
    )
    c0, c1, a, b, sigma, frequency = fit.x
    prediction = model(fit.x)
    ss_res = float(np.sum((y - prediction) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-30)
    omega_n = math.hypot(sigma, 2.0 * math.pi * frequency)
    return {
        "frequency_hz": float(frequency), "sigma_per_s": float(sigma),
        "damping_ratio": float(-sigma / omega_n),
        "phase_deg": float(math.degrees(math.atan2(b, a))), "r2": float(r2),
        "reliable": bool(fit.success and r2 >= 0.65 and frequency * t[-1] >= 3.0),
    }


def reconstruct_commands(
    time: np.ndarray, pid: dict[int, np.ndarray], c: dict[str, float]
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    a1, b0, b1 = c["ACT_A1"], c["ACT_B0"], c["ACT_B1"]
    limit = c["SURFACE_CORRECTION_LIMIT"]
    long_act = filter_first_order(np.clip(pid[9103] + pid[9104], -limit, limit), a1, b0, b1)
    lat_act = filter_first_order(np.clip(pid[9201] + pid[9202], -limit, limit), a1, b0, b1)
    dir_act = np.zeros_like(time)
    alt_drive = filter_first_order(pid[9101], c["ALT_LP_A1"], c["ALT_LP_B0"], c["ALT_LP_B1"])
    lift_raw = c["LIFT_HOLD_GAIN"] * (alt_drive + pid[9102])
    lift_act = filter_first_order(
        np.clip(lift_raw, -c["LIFT_CORRECTION_LIMIT"], c["LIFT_CORRECTION_LIMIT"]),
        a1, b0, b1,
    )
    positive = ((time >= c["BFF_RAP_START"]) & (time < c["BFF_RAP_START"] + 0.5 * c["BFF_RAP_DURATION"])).astype(float)
    negative = ((time >= c["BFF_RAP_START"] + 0.5 * c["BFF_RAP_DURATION"]) & (time < c["BFF_RAP_END"])).astype(float)
    rap_raw = c["BFF_RAP_AMPLITUDE"] * (positive - negative)
    rap_act = filter_first_order(rap_raw, a1, b0, b1)
    live = {
        "BFL": c["BODY_TRIM_SURFACE"] + long_act + dir_act,
        "BFR": c["BODY_TRIM_SURFACE"] + long_act - dir_act,
        "WF1L": c["TRIM_SURFACE"] + lat_act + lift_act,
        "WF1R": c["TRIM_SURFACE"] - lat_act + lift_act,
        "WF2L": c["TRIM_SURFACE"] + lat_act + lift_act,
        "WF2R": c["TRIM_SURFACE"] - lat_act + lift_act,
        "WF3L": c["TRIM_SURFACE"] + lift_act, "WF3R": c["TRIM_SURFACE"] + lift_act,
        "WF4L": np.full_like(time, c["TRIM_SURFACE"]),
        "WF4R": np.full_like(time, c["TRIM_SURFACE"]),
    }
    trigger = (time < c["SAS_OFF_START"]) | (time >= c["SAS_ON_START"])
    held = {
        name: sample_and_hold(value, trigger, value[0])
        for name, value in live.items()
    }
    # WF4 is frozen like every other surface; the prescribed doublet is the
    # only plant input during SAS-off and is added to its held release value.
    surfaces = dict(held)
    surfaces["WF4L"] = held["WF4L"] + rap_act
    surfaces["WF4R"] = held["WF4R"] + rap_act
    applied = {
        "longitudinal_feedback": np.where(trigger, long_act, 0.0),
        "lateral_feedback": np.where(trigger, lat_act, 0.0),
        "directional_feedback": np.where(trigger, dir_act, 0.0),
        "lift_feedback": np.where(trigger, lift_act, 0.0),
        "rap": rap_act,
        "safety": np.zeros_like(time),
    }
    actuator_states = {
        "longitudinal": long_act, "lateral": lat_act, "directional": dir_act,
        "lift": lift_act, "rap": rap_act,
    }
    return surfaces, applied, actuator_states, live


def write_window_csv(results: Path, stem: str, window: float, relative_time: np.ndarray, relative: dict[str, np.ndarray], use: np.ndarray) -> None:
    folder = results / "windows"
    folder.mkdir(exist_ok=True)
    columns = {
        "time_from_freeze_s": relative_time,
        "delta_theta_deg": relative["delta_theta"][use], "q_deg_s": relative["q"][use],
        "alpha_deg": relative["alpha_raw"][use], "delta_alpha_deg": relative["alpha"][use],
        "delta_Vz_mps": relative["delta_vz"][use], "delta_Z_m": relative["delta_z"][use],
        "swb1_modal_coordinate_relative": relative["swb1"][use],
        "swb1_modal_velocity": relative["swb1_velocity"][use],
        "left_tip_relative_m": relative["left_tip"][use], "right_tip_relative_m": relative["right_tip"][use],
        "symmetric_tip_relative_m": relative["symmetric_tip"][use],
    }
    name = f"{stem}_window_{int(round(1000 * window)):04d}ms.csv"
    with (folder / name).open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(columns); writer.writerows(zip(*columns.values()))


def write_pole_tables(results: Path, stem: str, rows: list[dict], clusters: dict[str, list[dict]]) -> None:
    if rows:
        with (results / f"{stem}_all_poles.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    flat = []
    for values in clusters.values():
        for value in values:
            row = dict(value)
            for key in ("windows_s", "orders", "detrends"):
                row[key] = "|".join(map(str, row[key]))
            flat.append(row)
    if flat:
        with (results / f"{stem}_pole_clusters.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(flat[0])); writer.writeheader(); writer.writerows(flat)


def build_trim_check(
    time: np.ndarray, indices: tuple[int, int], labels: tuple[str, str], pitch: np.ndarray,
    q: np.ndarray, vz: np.ndarray, position: np.ndarray, surface_commands: dict[str, np.ndarray],
    surface_kinematics: dict[str, np.ndarray],
    pid: dict[int, np.ndarray], actuator_states: dict[str, np.ndarray], residuals: dict,
) -> dict:
    result = {}
    for label, index in zip(labels, indices):
        force, moment = residuals[label]
        result[label] = {
            "time_s": float(time[index]), "theta_deg": float(pitch[index]), "q_deg_s": float(q[index]),
            "Vz_mps": float(vz[index]), "Z_m": float(position[index, 2] * IN_TO_M),
            "Fz_residual_N": float(force[2]), "My_residual_Nm": float(moment[1]),
            "surface_commands_deg": {name: float(value[index] * RAD_TO_DEG) for name, value in surface_commands.items()},
            "actual_surfaces_deg": {name: float(value[index]) for name, value in surface_kinematics.items()},
            "pid_outputs_deg": {name: float(pid[number][index] * RAD_TO_DEG) for number, name in zip(PID_LABELS, PID_NAMES)},
            "actuator_states_deg": {name: float(value[index] * RAD_TO_DEG) for name, value in actuator_states.items()},
        }
    return result


def make_identification_plots(
    stem: str, results: Path, time: np.ndarray, open_time: float, relative: dict[str, np.ndarray],
    local_valid: dict[float, bool], pole_rows: list[dict], candidate: dict, short_period: dict,
) -> None:
    plots = results / "plots"; plots.mkdir(exist_ok=True)
    use = (time >= open_time) & (time <= open_time + max(WINDOWS_S) + 1e-9)
    tau = time[use] - open_time
    fig, axes = plt.subplots(5, 2, figsize=(12, 14), sharex=True)
    raw_series = (
        ("theta_raw", "theta [deg]"), ("q", "q [deg/s]"), ("alpha_raw", "alpha [deg]"),
        ("vz_raw", "Vz [m/s]"),
        ("z_raw", "Z [m]"), ("swb1_raw", "SWB1"), ("swb1_velocity", "SWB1dot"),
        ("left_tip_raw", "left tip Z [m]"), ("right_tip_raw", "right tip Z [m]"),
    )
    for axis, (name, ylabel) in zip(axes.flat, raw_series):
        axis.plot(tau, relative[name][use]); axis.set_ylabel(ylabel); axis.grid(alpha=.2)
    axes[-1, 0].set_xlabel("time from freeze [s]"); axes[-1, 1].set_xlabel("time from freeze [s]")
    fig.tight_layout(); fig.savefig(plots / f"{stem}_raw_histories.png", dpi=170); plt.close(fig)

    fig, axes = plt.subplots(5, 2, figsize=(12, 14), sharex=True)
    series = (
        ("delta_theta", "Δtheta [deg]"), ("q", "q [deg/s]"), ("alpha", "Δalpha [deg]"),
        ("delta_vz", "ΔVz [m/s]"),
        ("delta_z", "ΔZ [m]"), ("swb1", "ΔSWB1"), ("swb1_velocity", "SWB1dot"),
        ("left_tip", "left tip ΔZ [m]"), ("right_tip", "right tip ΔZ [m]"),
    )
    for axis, (name, ylabel) in zip(axes.flat, series):
        axis.plot(tau, relative[name][use]); axis.set_ylabel(ylabel); axis.grid(alpha=.2)
    axes[-1, 0].set_xlabel("time from freeze [s]"); axes[-1, 1].set_xlabel("time from freeze [s]")
    fig.tight_layout(); fig.savefig(plots / f"{stem}_relative_histories.png", dpi=170); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    tip = detrend(relative["symmetric_tip"][use]); mv = detrend(relative["swb1"][use])
    tip /= max(np.max(np.abs(tip)), 1e-12); mv /= max(np.max(np.abs(mv)), 1e-12)
    ax.plot(tau, tip, label="symmetric tip deformation, detrended/normalized"); ax.plot(tau, mv, label="SWB1 detrended/normalized")
    ax.set_xlabel("time from freeze [s]"); ax.set_ylabel("normalized"); ax.legend(); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(plots / f"{stem}_tip_vs_swb1.png", dpi=170); plt.close(fig)

    for signal_name, key in (("symmetric_tip", "symmetric_tip"), ("swb1", "swb1")):
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
        for axis, window in zip(axes.flat, WINDOWS_S):
            window_use = (time >= open_time) & (time <= open_time + window + 1e-9)
            tw = time[window_use] - open_time; y = relative[key][window_use]
            if len(tw) < 20:
                continue
            dt = float(np.median(np.diff(tw))); nfft = max(4096, 2 ** int(np.ceil(np.log2(len(tw) * 16))))
            frequency = np.fft.rfftfreq(nfft, dt)
            show = frequency <= min(15.0, 0.45 / dt)
            for kind in DETRENDS:
                value = detrend_variant(tw, y, kind) * np.hanning(len(y))
                amplitude = np.abs(np.fft.rfft(value, nfft)); amplitude /= max(np.max(amplitude[show]), 1e-30)
                axis.plot(frequency[show], amplitude[show], label=kind)
            axis.set_title(f"{window:.2f} s — {'valid' if local_valid.get(window) else 'departure'}")
            axis.grid(alpha=.2)
        axes[0, 0].legend(); axes[-1, 0].set_xlabel("frequency [Hz]"); axes[-1, 1].set_xlabel("frequency [Hz]"); axes[-1, 2].set_xlabel("frequency [Hz]")
        fig.supylabel("normalized spectrum"); fig.tight_layout(); fig.savefig(plots / f"{stem}_spectrum_{signal_name}.png", dpi=170); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = {"swb1": "tab:orange", "swb1_velocity": "tab:green", "symmetric_tip": "tab:red"}
    for name, color in colors.items():
        rows = [row for row in pole_rows if row["signal"] == name and row["is_complex"] and row["local_window_valid"]]
        if rows:
            ax.scatter([row["frequency_hz"] for row in rows], [row["sigma_per_s"] for row in rows], s=8, alpha=.22, color=color, label=name)
    if candidate.get("matches"):
        for name, value in candidate["matches"].items():
            ax.scatter(value["frequency_hz"], value["sigma_per_s"], marker="*", s=180, edgecolor="k", color=colors[name])
    ax.axhline(0.0, color="k", lw=1); ax.set_xlim(left=0.0); ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("sigma [1/s]"); ax.legend(); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(plots / f"{stem}_matrix_pencil_poles.png", dpi=170); plt.close(fig)

    rigid_duration = max((window for window, valid in local_valid.items() if valid), default=0.0)
    rigid_use = (time >= open_time) & (time <= open_time + rigid_duration + 1e-9)
    rigid_tau = time[rigid_use] - open_time
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    for key, label in (("q", "q"), ("alpha", "effective alpha"), ("delta_theta", "pitch")):
        value = detrend(relative[key][rigid_use]); value /= max(np.max(np.abs(value)), 1e-12)
        axes[0].plot(rigid_tau, value, label=label)
    axes[0].set_xlabel("time from identification start [s]"); axes[0].set_ylabel("normalized")
    axes[0].legend(); axes[0].grid(alpha=.2)
    rigid_colors = {"q": "tab:blue", "alpha": "tab:purple", "pitch": "tab:orange"}
    for name, color in rigid_colors.items():
        rows = [row for row in pole_rows if row["signal"] == name and row["is_complex"] and row["local_window_valid"]]
        if rows:
            axes[1].scatter([row["frequency_hz"] for row in rows], [row["sigma_per_s"] for row in rows], s=8, alpha=.22, color=color, label=name)
    if short_period.get("matches"):
        for name, value in short_period["matches"].items():
            axes[1].scatter(value["frequency_hz"], value["sigma_per_s"], marker="*", s=180, edgecolor="k", color=rigid_colors[name])
    axes[1].axhline(0., color="k", lw=1); axes[1].set_xlim(left=0.)
    axes[1].set_xlabel("frequency [Hz]"); axes[1].set_ylabel("sigma [1/s]")
    axes[1].legend(); axes[1].grid(alpha=.2)
    fig.tight_layout(); fig.savefig(plots / f"{stem}_short_period_q_alpha.png", dpi=170); plt.close(fig)


def analyze_case(
    nc_path: Path,
    results_dir: Path | None = None,
    tracking_reference_frequency_hz: float | None = None,
) -> dict:
    nc_path = Path(nc_path)
    results = Path(results_dir or nc_path.parent)
    results.mkdir(parents=True, exist_ok=True)
    input_path = nc_path.with_suffix(".mbd")
    text = input_path.read_text()
    # run_case.py inlines the scheduled constants so each result remains
    # auditable after later tuning.  Legacy cases still need the shared file.
    if "set: const real NASTRAN_DLM_ROM_K7" not in text:
        text += "\n" + (ROOT / "INCLUDE/setconst.mbd").read_text()
    c = model_constants(text)
    swb = identify_swb1_from_fem()
    with Dataset(nc_path) as data:
        time = arr(data, "time")
        if len(time) == 0 or float(time[-1]) < c["SAS_ON_START"] - 1.5 * c["TIME_STEP"]:
            last_time = None if len(time) == 0 else float(time[-1])
            raise RuntimeError(
                f"{nc_path.name}: run incompleta (ultimo tempo={last_time}, "
                f"riaccensione SAS={c['SAS_ON_START']}); nessuna identificazione BFF valida"
            )
        position = arr(data, "node.struct.990000.X")
        velocity = arr(data, "node.struct.990000.XP")
        phi = arr(data, "node.struct.990000.Phi")
        omega = arr(data, "node.struct.990000.Omega")
        base_rotation = Rotation.from_rotvec(phi)
        euler = np.unwrap(base_rotation.as_euler("xyz"), axis=0)
        air_velocity = np.zeros_like(velocity)
        air_velocity[:, 0] = c["VINF"]
        relative_air_body = base_rotation.inv().apply(air_velocity - velocity)
        alpha_deg = np.arctan2(relative_air_body[:, 2], relative_air_body[:, 0]) * RAD_TO_DEG
        modal = arr(data, "elem.joint.5.a")
        modal_velocity = arr(data, "elem.joint.5.aPrime")
        pid = {label: arr(data, f"elem.loadable.{label}.output") for label in PID_LABELS}
        kine = surface_kinematics(data)
        # Physical bending observable: tip vertical motion in the floating
        # body frame.  Removing base translation and rotation prevents rigid
        # heave/pitch from being mistaken for wing bending.
        tip_l_position = arr(data, "node.struct.990020.X")
        tip_r_position = arr(data, "node.struct.991020.X")
        tip_l = base_rotation.inv().apply(tip_l_position - position)[:, 2] * IN_TO_M
        tip_r = base_rotation.inv().apply(tip_r_position - position)[:, 2] * IN_TO_M
        pre_rap_index = max(0, int(np.searchsorted(time, c["BFF_RAP_START"])) - 1)
        pre_freeze_index = max(0, int(np.searchsorted(time, c["SAS_OFF_START"])) - 1)
        residuals = {}
        for label, index in (("pre_rap", pre_rap_index), ("pre_freeze", pre_freeze_index)):
            residuals[label] = point_force_moment_residual(data, index, position[index])

    surfaces, applied, actuator_states, live_surfaces = reconstruct_commands(time, pid, c)
    column = int(swb["modal_column"])
    trim_use = (time >= c["SAS_OFF_START"] - 1.0) & (time < c["SAS_OFF_START"])
    release_index = int(np.searchsorted(time, c["SAS_OFF_START"]))
    recovery_index = int(np.searchsorted(time, c["SAS_ON_START"]))
    id_index = int(np.searchsorted(time, c["IDENTIFICATION_START"]))
    id_time = float(time[id_index])
    id_end = min(c["IDENTIFICATION_END"], float(time[-1]))
    id_use = (time >= id_time) & (time <= id_end + 1e-9)
    surface_command_ranges = {
        name: float(np.ptp(value[id_use]) * RAD_TO_DEG) for name, value in surfaces.items()
    }
    surface_kinematic_ranges = {name: float(np.ptp(value[id_use])) for name, value in kine.items()}
    last_closed_loop_step_changes = {
        name: float(abs(value[release_index] - value[max(release_index - 1, 0)]) * RAD_TO_DEG)
        for name, value in surfaces.items()
    }
    # MBDyn's sample-and-hold updates its stored value after convergence of a
    # triggered step.  The first open sample must therefore equal the live
    # command from the last closed step, not the extrapolated live value at
    # the new time.  This distinguishes a real hold discontinuity from the
    # normal actuator motion over the final closed-loop time step.
    command_hold_capture_errors = {
        name: float(abs(surfaces[name][release_index] - live_surfaces[name][max(release_index - 1, 0)]) * RAD_TO_DEG)
        for name in surfaces
    }
    command_release_changes = {
        name: float(abs(value[min(release_index + 1, len(value) - 1)] - value[release_index]) * RAD_TO_DEG)
        for name, value in surfaces.items()
    }
    kinematic_release_changes = {
        name: float(abs(value[min(release_index + 1, len(value) - 1)] - value[release_index]))
        for name, value in kine.items()
    }
    applied_max = {name: float(np.max(np.abs(value[id_use]))) for name, value in applied.items()}
    no_notch = "NOTCH" not in text and "BFF_BS" not in text
    no_safety_path = "SAFETY_ACT_DRIVE" not in text
    excitation_ended = c["BFF_RAP_END"] < c["IDENTIFICATION_START"]
    # The open-loop audit acts on all ten commanded hinge deflections.  A hinge
    # angle reconstructed from flexible node rotations also contains elastic
    # deformation and must remain a diagnostic, not a false feedback detector.
    command_hold_tolerance_deg = c.get("SURFACE_HOLD_TOLERANCE_DEG", 0.001)
    held_surface_names = tuple(surfaces)
    all_surfaces_constant = max(surface_command_ranges[name] for name in held_surface_names) <= command_hold_tolerance_deg
    release_continuous = max(
        max(command_hold_capture_errors[name], command_release_changes[name])
        for name in held_surface_names
    ) <= 1e-6
    all_feedback_off = max(
        applied_max[name] for name in (
            "longitudinal_feedback", "lateral_feedback",
            "directional_feedback", "lift_feedback", "safety",
        )
    ) <= 1e-12
    sas_modal_force = np.where(
        (time < c["SAS_OFF_START"]) | (time >= c["SAS_ON_START"]),
        -c.get("SAS_MODAL_DAMPING", 0.0) * modal_velocity[:, column],
        0.0,
    )
    sas_modal_force_max_open = float(np.max(np.abs(sas_modal_force[id_use])))
    sas_modal_damper_off = sas_modal_force_max_open <= 1e-12
    completed = bool(time[-1] >= c["FINAL_TIME"] - 1.5 * c["TIME_STEP"])
    bounded_window = c["SAS_OFF_START"] < c["BFF_RAP_START"] < c["BFF_RAP_END"] < c["IDENTIFICATION_START"] < c["IDENTIFICATION_END"] < c["SAS_ON_START"]
    open_loop_valid = bool(
        completed and bounded_window and no_notch and no_safety_path
        and excitation_ended and all_surfaces_constant
        and release_continuous and all_feedback_off and sas_modal_damper_off
    )

    roll_deg = euler[:, 0] * RAD_TO_DEG
    pitch_deg = euler[:, 1] * RAD_TO_DEG
    p_deg_s = omega[:, 0] * RAD_TO_DEG
    q_deg_s = omega[:, 1] * RAD_TO_DEG
    vz_mps = velocity[:, 2] * IN_TO_M
    trim_stationary = bool(
        np.count_nonzero(trim_use) >= 50
        and np.std(pitch_deg[trim_use]) < 0.05
        and np.std(q_deg_s[trim_use]) < 0.20
        and np.std(roll_deg[trim_use]) < 0.10
        and np.std(p_deg_s[trim_use]) < 0.30
        and np.std(vz_mps[trim_use]) < 0.05
        and abs(q_deg_s[release_index]) < 0.10
        and abs(p_deg_s[release_index]) < 0.20
        and abs(vz_mps[release_index]) < 0.05
    )
    freeze_values = {
        "theta_deg": pitch_deg[id_index], "q_deg_s": q_deg_s[id_index],
        "alpha_deg": alpha_deg[id_index],
        "Vz_mps": vz_mps[id_index], "Z_m": position[id_index, 2] * IN_TO_M,
        "swb1": modal[id_index, column], "swb1_velocity": modal_velocity[id_index, column],
        "left_tip_m": tip_l[id_index], "right_tip_m": tip_r[id_index],
    }
    relative = {
        "theta_raw": pitch_deg,
        "delta_theta": pitch_deg - freeze_values["theta_deg"],
        "q": q_deg_s,
        "alpha_raw": alpha_deg,
        "alpha": alpha_deg - freeze_values["alpha_deg"],
        "vz_raw": vz_mps,
        "delta_vz": vz_mps - freeze_values["Vz_mps"],
        "z_raw": position[:, 2] * IN_TO_M,
        "delta_z": position[:, 2] * IN_TO_M - freeze_values["Z_m"],
        "swb1_raw": modal[:, column],
        "swb1": modal[:, column] - freeze_values["swb1"],
        "swb1_velocity": modal_velocity[:, column],
        "left_tip_raw": tip_l,
        "right_tip_raw": tip_r,
        "left_tip": tip_l - freeze_values["left_tip_m"],
        "right_tip": tip_r - freeze_values["right_tip_m"],
        "symmetric_tip": 0.5 * ((tip_l - freeze_values["left_tip_m"]) + (tip_r - freeze_values["right_tip_m"])),
    }
    local_valid: dict[float, bool] = {}
    window_metrics: dict[str, dict] = {}
    for window in WINDOWS_S:
        use = (time >= id_time) & (time <= min(id_time + window, id_end) + 1e-9)
        if id_time + window > id_end + 1.5 * c["TIME_STEP"] or np.count_nonzero(use) < 20:
            local_valid[window] = False
            continue
        metrics = {
            "max_abs_delta_theta_deg": float(np.max(np.abs(relative["delta_theta"][use]))),
            "max_abs_q_deg_s": float(np.max(np.abs(relative["q"][use]))),
            "max_abs_delta_vz_mps": float(np.max(np.abs(relative["delta_vz"][use]))),
            "max_abs_delta_z_m": float(np.max(np.abs(relative["delta_z"][use]))),
        }
        valid = bool(
            metrics["max_abs_delta_theta_deg"] <= 10.0
            and metrics["max_abs_q_deg_s"] <= 40.0
            and metrics["max_abs_delta_vz_mps"] <= 8.0
            and metrics["max_abs_delta_z_m"] <= 3.0
        )
        local_valid[window] = valid
        metrics["local_valid"] = valid
        window_metrics[f"{window:.2f}"] = metrics
        write_window_csv(results, nc_path.stem, window, time[use] - id_time, relative, use)
    valid_windows = [window for window, valid in local_valid.items() if valid]
    selected_window = max(valid_windows) if valid_windows else None
    open_use = id_use
    pencil_signals = {
        "symmetric_tip": relative["symmetric_tip"][open_use],
        "swb1": relative["swb1"][open_use],
        "swb1_velocity": relative["swb1_velocity"][open_use],
        "q": relative["q"][open_use],
        "alpha": relative["alpha"][open_use],
        "pitch": relative["delta_theta"][open_use],
    }
    pole_rows, clusters, candidate, short_period = identify_multimodal(
        time[open_use], pencil_signals, local_valid,
        tracking_reference_frequency_hz=tracking_reference_frequency_hz,
    )
    symmetric_fem = [item for item in swb["candidates"] if item[2] < 0.20]
    if short_period.get("accepted") and symmetric_fem:
        fem_frequency, fem_mode, fem_symmetry = min(
            symmetric_fem, key=lambda item: abs(item[0] - short_period["frequency_hz"])
        )
        short_period["nearest_symmetric_fem_mode"] = fem_mode
        short_period["nearest_symmetric_fem_dry_frequency_hz"] = fem_frequency
        short_period["nearest_symmetric_fem_frequency_difference_hz"] = abs(fem_frequency - short_period["frequency_hz"])
        short_period["nearest_symmetric_fem_symmetry_residual"] = fem_symmetry
    write_pole_tables(results, nc_path.stem, pole_rows, clusters)
    secondary_fit = {"symmetric_tip": {}, "swb1": {}}
    correlation = math.nan
    if selected_window is not None:
        selected = (time >= id_time) & (time <= id_time + selected_window + 1e-9)
        secondary_fit = {
            "symmetric_tip": fit_damped_sine(time[selected], relative["symmetric_tip"][selected]),
            "swb1": fit_damped_sine(time[selected], relative["swb1"][selected]),
        }
        correlation = float(np.corrcoef(detrend(relative["symmetric_tip"][selected]), detrend(relative["swb1"][selected]))[0, 1])
    sigma = candidate.get("sigma_per_s") if candidate.get("accepted", False) else None
    frequency = candidate.get("frequency_hz") if candidate.get("accepted", False) else None
    damping = candidate.get("damping_ratio") if candidate.get("accepted", False) else None
    leakage_limit_deg = 0.025
    feedback_names = (
        "longitudinal_feedback", "lateral_feedback",
        "directional_feedback", "lift_feedback", "safety",
    )
    control_band_leakage_deg: dict[str, float | None] = {}
    if frequency is not None:
        for name in feedback_names:
            control_band_leakage_deg[name] = detrended_harmonic_amplitude(
                time[id_use], applied[name][id_use], frequency,
            ) * RAD_TO_DEG
    else:
        control_band_leakage_deg = {
            name: None for name in feedback_names
        }
    all_feedback_band_clean = bool(
        frequency is not None
        and max(value for value in control_band_leakage_deg.values() if value is not None) <= leakage_limit_deg
    )
    mode9_frequency = next((item[0] for item in swb["candidates"] if item[1] == 9), None)
    q_bff_amplitude = None if frequency is None else detrended_harmonic_amplitude(time[id_use], q_deg_s[id_use], frequency)
    q_mode9_amplitude = None if mode9_frequency is None else detrended_harmonic_amplitude(time[id_use], q_deg_s[id_use], mode9_frequency)
    high_frequency_q_audit = {
        "bff_frequency_hz": frequency,
        "q_amplitude_at_bff_deg_s": q_bff_amplitude,
        "first_symmetric_torsion_fem_frequency_hz": mode9_frequency,
        "q_amplitude_at_first_symmetric_torsion_deg_s": q_mode9_amplitude,
        "torsion_to_bff_q_amplitude_ratio": None if not q_bff_amplitude else q_mode9_amplitude / q_bff_amplitude,
        "rap_doublet_spectral_peak_hz": 0.742 / c["BFF_RAP_DURATION"],
        "rap_target_frequency_hz": c.get("BFF_RAP_TARGET_FREQUENCY"),
    }
    identification_valid = bool(
        open_loop_valid and trim_stationary and candidate.get("accepted", False)
        and all_feedback_band_clean
    )
    altitude_m = position[:, 2] * IN_TO_M
    release_altitude = float(altitude_m[release_index])
    sas_off_use = (time >= c["SAS_OFF_START"]) & (time <= c["SAS_ON_START"] + 1e-9)
    recovery_use = time >= c["SAS_ON_START"]
    final_use = time >= max(c["SAS_ON_START"], time[-1] - 1.0)
    recovery = {
        "release_theta_deg": float(pitch_deg[release_index]),
        "release_q_deg_s": float(q_deg_s[release_index]),
        "release_roll_deg": float(roll_deg[release_index]),
        "release_p_deg_s": float(p_deg_s[release_index]),
        "release_vz_mps": float(vz_mps[release_index]),
        "altitude_change_at_rap_start_m": float(altitude_m[pre_rap_index] - release_altitude),
        "altitude_change_at_sas_reengage_m": float(altitude_m[recovery_index] - release_altitude),
        "max_altitude_loss_m": float(max(0.0, release_altitude - np.min(altitude_m[release_index:]))),
        "max_altitude_loss_during_sas_off_m": float(max(0.0, release_altitude - np.min(altitude_m[sas_off_use]))),
        "max_abs_pitch_deg": float(np.max(np.abs(pitch_deg[release_index:]))),
        "max_abs_q_deg_s": float(np.max(np.abs(q_deg_s[release_index:]))),
        "max_abs_roll_deg": float(np.max(np.abs(roll_deg[release_index:]))),
        "max_abs_p_deg_s": float(np.max(np.abs(p_deg_s[release_index:]))),
        "final_altitude_error_m": float(np.mean(altitude_m[final_use]) - release_altitude),
        "final_pitch_deg": float(np.mean(pitch_deg[final_use])),
        "final_q_deg_s": float(np.mean(q_deg_s[final_use])),
        "final_roll_deg": float(np.mean(roll_deg[final_use])),
        "final_p_deg_s": float(np.mean(p_deg_s[final_use])),
        "final_vz_mps": float(np.mean(vz_mps[final_use])),
        "sas_reengaged": bool(np.max(np.abs(applied["longitudinal_feedback"][recovery_use])) > 1e-8),
    }
    final_altitude_tolerance_m = 0.25
    recovery["satisfactory"] = bool(
        abs(recovery["release_q_deg_s"]) < 0.10
        and abs(recovery["release_vz_mps"]) < 0.05
        and recovery["max_altitude_loss_m"] < 0.50
        and recovery["max_altitude_loss_during_sas_off_m"] < 0.05
        and recovery["max_abs_pitch_deg"] < 3.0
        and recovery["max_abs_q_deg_s"] < 10.0
        and recovery["max_abs_roll_deg"] < 5.0
        and recovery["max_abs_p_deg_s"] < 10.0
        and abs(recovery["final_altitude_error_m"]) < final_altitude_tolerance_m
        and abs(recovery["final_q_deg_s"]) < 0.05
        and abs(recovery["final_p_deg_s"]) < 0.05
        and abs(recovery["final_vz_mps"]) < 0.05
        and recovery["sas_reengaged"]
    )
    sas_limit = c["SURFACE_CORRECTION_LIMIT"]
    sas_saturation_fraction = float(np.mean(np.abs(actuator_states["longitudinal"]) >= sas_limit - 1e-9))
    q_dyn_lbf_in2 = 0.5 * c["RHO_AIR"] * c["VINF"]**2
    required_cz = WEIGHT_LBF / (q_dyn_lbf_in2 * S_REF_IN2)
    pre_freeze_fz_residual_n = float(residuals["pre_freeze"][0][2])
    achieved_cz = (WEIGHT_LBF + pre_freeze_fz_residual_n / LBF_TO_N) / (q_dyn_lbf_in2 * S_REF_IN2)
    trim_lift_audit = {
        "release_alpha_deg": float(alpha_deg[release_index]),
        "release_pitch_deg": float(pitch_deg[release_index]),
        "required_CZ_for_weight": float(required_cz),
        "achieved_CZ_pre_freeze": float(achieved_cz),
        "vertical_force_residual_N_pre_freeze": pre_freeze_fz_residual_n,
        "vertical_force_residual_fraction_of_weight": float(pre_freeze_fz_residual_n / (WEIGHT_LBF * LBF_TO_N)),
        "nastran_rigid_zero_surface_alpha_for_weight_deg": float(required_cz / NASTRAN_CZ_ALPHA_PER_DEG),
        "note": "The NASTRAN alpha is a zero-surface rigid-polar reference; the simulated trim includes surface bias and static aeroelastic deformation.",
    }
    aerobody_text = (ROOT / "INCLUDE/aerobody.mbd").read_text()
    rendered_model_text = nc_path.with_suffix(".mbd").read_text(errors="ignore")
    unsteady_count = aerobody_text.count("theodorsen, c81")
    three_quarter_terms = aerobody_text.count("-0.5*(")
    dlm_rom_enabled = "force: NASTRAN_DLM_ROM_FORCE, modal, MODAL_JOINT" in rendered_model_text
    try:
        nastran_reference = interpolate_at_velocity(parse_flutter_point(DEFAULT_F06, point=7), c["V_INF"])
    except (OSError, ValueError):
        nastran_reference = None
    if nastran_reference is not None and frequency is not None and sigma is not None:
        nastran_reference["mbdyn_frequency_error_hz"] = float(frequency - nastran_reference["frequency_hz"])
        nastran_reference["mbdyn_sigma_error_per_s"] = float(sigma - nastran_reference["sigma_per_s"])
    trim_check = build_trim_check(
        time, (pre_rap_index, pre_freeze_index), ("pre_rap", "pre_freeze"),
        pitch_deg, q_deg_s, vz_mps, position, surfaces, kine, pid, actuator_states, residuals,
    )
    (results / f"{nc_path.stem}_trim_check.json").write_text(json.dumps(trim_check, indent=2))
    summary = {
        "case": nc_path.stem, "velocity_mps": c["V_INF"], "test_sequence": "bounded_nasa_style_sas_off",
        "sas_off_start_s": c["SAS_OFF_START"], "sas_on_start_s": c["SAS_ON_START"],
        "open_loop_valid": open_loop_valid,
        "identification_valid": identification_valid,
        "identification_scope": "bounded local candidate; run_sweep.py adds adjacent-speed continuity",
        "trim_stationary": trim_stationary,
        "no_notch": no_notch, "no_safety_path": no_safety_path,
        "all_surfaces_held": all_surfaces_constant,
        "all_surface_feedback_off": all_feedback_off,
        "sas_modal_damper_off": sas_modal_damper_off,
        "sas_modal_generalized_force_max_open": sas_modal_force_max_open,
        "bff_pitch_q_feedback_off": all_feedback_off,
        "retained_controls": [],
        "feedback_bff_band_leakage_deg": control_band_leakage_deg,
        "retained_control_bff_band_leakage_deg": control_band_leakage_deg,
        "retained_control_bff_band_limit_deg": leakage_limit_deg,
        "all_feedback_band_clean": all_feedback_band_clean,
        "high_frequency_q_audit": high_frequency_q_audit,
        "retained_controls_band_clean": all_feedback_band_clean,
        "feedback_applied_max_rad": applied_max,
        "surface_range_open_loop_deg": surface_command_ranges,
        "surface_command_audit_tolerance_deg": command_hold_tolerance_deg,
        "surface_kinematic_range_open_loop_deg": surface_kinematic_ranges,
        "last_closed_loop_step_change_deg": last_closed_loop_step_changes,
        "hold_capture_error_deg": command_hold_capture_errors,
        "release_change_open_to_next_deg": command_release_changes,
        "kinematic_release_change_open_to_next_deg": kinematic_release_changes,
        "swb1_fem_mode": swb["fem_mode"], "swb1_dry_frequency_hz": swb["dry_frequency_hz"],
        "swb1_symmetry_residual": swb["symmetry_residual"],
        "adaptive_windows": window_metrics, "selected_window_s": selected_window,
        "matrix_pencil_candidate": candidate,
        "short_period_candidate": short_period,
        "identification_channels": {
            "bff": ["symmetric_tip_body_frame", "swb1", "swb1_velocity"],
            "short_period": ["q", "effective_alpha", "pitch"],
            "effective_alpha_definition": "atan2(Vair_body_z, Vair_body_x), Vair_global=[VINF,0,0]-base_velocity",
        },
        "matrix_pencil_clusters": clusters,
        "secondary_single_mode_fit": secondary_fit,
        "symmetric_tip_swb1_linear_detrend_correlation": correlation,
        "recovery": recovery,
        "recovery_final_altitude_tolerance_m": final_altitude_tolerance_m,
        "sas_saturation_fraction": sas_saturation_fraction,
        "trim_lift_audit": trim_lift_audit,
        "aerodynamic_model": {
            "type": "hybrid sectional Wagner/Theodorsen+C81 with NASTRAN SOL145 reduced-order mode-7 correction",
            "theodorsen_element_count": unsteady_count,
            "total_aerodynamic_element_count": len(aerodynamic_spans()),
            "all_elements_unsteady": unsteady_count == len(aerodynamic_spans()),
            "three_quarter_chord_collocation_element_count": three_quarter_terms // 2,
            "all_collocation_points_at_three_quarter_chord": three_quarter_terms == 2 * len(aerodynamic_spans()),
            "nastran_dlm_rom_correction": {
                "enabled": dlm_rom_enabled,
                "fem_mode": 7,
                "stiffness_correction_per_s2": c.get("NASTRAN_DLM_ROM_K7"),
                "damping_correction_per_s": c.get("NASTRAN_DLM_ROM_C7"),
                "trim_modal_coordinate": c.get("NASTRAN_DLM_ROM_Q7_EQ"),
                "initialization_ramp_s": c.get("NASTRAN_DLM_ROM_RAMP_TIME"),
                "reference": "released X-56A SOL 145 point 7; correction is aerodynamic, not SAS feedback",
            },
        },
        "sigma_swb1_per_s": sigma, "frequency_swb1_hz": frequency,
        "damping_ratio_swb1": damping,
        "classification": "no_local_candidate" if not identification_valid else ("local_unstable_candidate" if sigma > 0.0 else "local_stable_candidate" if sigma < 0.0 else "local_neutral_candidate"),
        "trim_check_file": f"{nc_path.stem}_trim_check.json",
        "nastran_bff_point7_interpolated": nastran_reference,
    }
    (results / f"{nc_path.stem}_summary.json").write_text(json.dumps(summary, indent=2))

    columns: dict[str, np.ndarray] = {
        "time_s": time, "X_m": position[:, 0] * IN_TO_M, "Y_m": position[:, 1] * IN_TO_M, "Z_m": position[:, 2] * IN_TO_M,
        "roll_deg": roll_deg, "p_deg_s": p_deg_s,
        "pitch_deg": pitch_deg, "q_deg_s": q_deg_s, "vertical_velocity_mps": vz_mps,
        "alpha_deg": alpha_deg, "delta_alpha_deg": relative["alpha"],
        "time_from_identification_s": time - id_time,
        "delta_theta_deg": relative["delta_theta"], "delta_Vz_mps": relative["delta_vz"], "delta_Z_m": relative["delta_z"],
        "swb1_modal_coordinate": modal[:, column], "swb1_modal_velocity": modal_velocity[:, column],
        "delta_swb1_modal_coordinate": relative["swb1"],
        "left_wingtip_Z_m": tip_l, "right_wingtip_Z_m": tip_r, "symmetric_wingtip_Z_m": 0.5 * (tip_l + tip_r),
        "left_wingtip_relative_m": relative["left_tip"], "right_wingtip_relative_m": relative["right_tip"],
        "symmetric_wingtip_relative_m": relative["symmetric_tip"],
        "sas_active_gate": ((time < c["SAS_OFF_START"]) | (time >= c["SAS_ON_START"])).astype(float),
    }
    for name, value in surfaces.items():
        columns[f"surface_command_{name}_deg"] = value * RAD_TO_DEG
        columns[f"surface_kinematic_{name}_deg"] = kine[name]
    for name, value in applied.items():
        columns[f"applied_{name}_deg"] = value * RAD_TO_DEG
    for name, value in actuator_states.items():
        columns[f"actuator_state_{name}_deg"] = value * RAD_TO_DEG
    for label, name in zip(PID_LABELS, PID_NAMES):
        columns[f"pid_raw_{name}_deg"] = pid[label] * RAD_TO_DEG
    with (results / f"{nc_path.stem}_timeseries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
    make_case_plots(nc_path.stem, results, time, id_use, pitch_deg, q_deg_s, modal[:, column], relative["symmetric_tip"], surfaces, pid, applied, c, summary)
    make_identification_plots(nc_path.stem, results, time, id_time, relative, local_valid, pole_rows, candidate, short_period)
    make_test_dashboard(
        nc_path.stem, results, time, altitude_m - release_altitude, roll_deg, p_deg_s, pitch_deg, q_deg_s,
        vz_mps, alpha_deg, modal[:, column], relative["symmetric_tip"], surfaces, applied, c, summary,
    )
    return summary


def make_case_plots(stem: str, results: Path, time: np.ndarray, id_use: np.ndarray, pitch: np.ndarray, q: np.ndarray, mode: np.ndarray, symmetric_tip: np.ndarray, surfaces: dict[str, np.ndarray], pid: dict[int, np.ndarray], applied: dict[str, np.ndarray], c: dict[str, float], summary: dict) -> None:
    plots = results / "plots"
    plots.mkdir(exist_ok=True)
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(time, pitch, label="pitch θ"); ax[0].plot(time, q, label="pitch rate q"); ax[0].legend(); ax[0].set_ylabel("deg, deg/s")
    tipn = detrend(symmetric_tip[id_use]); mn = detrend(mode[id_use]); tipn /= max(np.max(np.abs(tipn)), 1e-12); mn /= max(np.max(np.abs(mn)), 1e-12)
    ax[1].plot(time[id_use], tipn, label="symmetric tip deformation normalized"); ax[1].plot(time[id_use], mn, label=f"SWB1 FEM {summary['swb1_fem_mode']} normalized"); ax[1].legend(); ax[1].set_ylabel("normalized"); ax[1].set_xlabel("time [s]")
    for axis in ax:
        axis.axvspan(c["SAS_OFF_START"], c["SAS_ON_START"], color="tab:red", alpha=0.10)
        axis.axvline(c["SAS_ON_START"], color="tab:green", ls="--", lw=1)
    fig.tight_layout(); fig.savefig(plots / f"{stem}_free_response.png", dpi=170); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for name, value in surfaces.items(): ax.plot(time, value * RAD_TO_DEG, label=name)
    ax.axvspan(c["SAS_OFF_START"], c["SAS_ON_START"], color="tab:red", alpha=0.10, label="SAS off")
    ax.axvspan(c["IDENTIFICATION_START"], c["IDENTIFICATION_END"], color="0.75", alpha=0.35, label="identification")
    ax.set_xlabel("time [s]"); ax.set_ylabel("surface command [deg]"); ax.legend(ncol=5, fontsize=8); fig.tight_layout(); fig.savefig(plots / f"{stem}_surfaces.png", dpi=170); plt.close(fig)

    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for label, name in zip(PID_LABELS, PID_NAMES): ax[0].plot(time, pid[label] * RAD_TO_DEG, label=name)
    for name, value in applied.items(): ax[1].plot(time, value * RAD_TO_DEG, label=name)
    for axis in ax:
        axis.axvspan(c["SAS_OFF_START"], c["SAS_ON_START"], color="tab:red", alpha=0.10)
    ax[0].set_ylabel("raw PID [deg]"); ax[1].set_ylabel("applied [deg]"); ax[1].set_xlabel("time [s]"); ax[0].legend(ncol=3, fontsize=7); ax[1].legend(ncol=3, fontsize=8); fig.tight_layout(); fig.savefig(plots / f"{stem}_controllers.png", dpi=170); plt.close(fig)


def make_test_dashboard(
    stem: str, results: Path, time: np.ndarray, altitude_delta: np.ndarray,
    roll: np.ndarray, p: np.ndarray, pitch: np.ndarray, q: np.ndarray,
    vz: np.ndarray, alpha: np.ndarray, swb1: np.ndarray, symmetric_tip: np.ndarray,
    surfaces: dict[str, np.ndarray], applied: dict[str, np.ndarray],
    c: dict[str, float], summary: dict,
) -> None:
    plots = results / "plots"
    plots.mkdir(exist_ok=True)
    fig, axes = plt.subplots(5, 1, figsize=(12, 14), sharex=True, constrained_layout=True)
    axes[0].plot(time, altitude_delta, label="altitude delta", color="tab:blue")
    axes[0].plot(time, vz, label="vertical speed", color="tab:orange", alpha=0.85)
    axes[0].axhline(0.0, color="0.25", lw=0.7)
    axes[0].set_ylabel("m, m/s"); axes[0].legend(loc="best")

    axes[1].plot(time, pitch, label="pitch", color="tab:blue")
    axes[1].plot(time, q, label="pitch rate", color="tab:orange")
    axes[1].plot(time, alpha, label="effective alpha", color="tab:green", alpha=.8)
    axes[1].axhline(0.0, color="0.25", lw=0.7)
    axes[1].set_ylabel("deg, deg/s"); axes[1].legend(loc="best")

    axes[2].plot(time, roll, label="roll angle", color="tab:purple")
    axes[2].plot(time, p, label="roll rate", color="tab:brown")
    axes[2].axhline(0.0, color="0.25", lw=0.7)
    axes[2].set_ylabel("deg, deg/s"); axes[2].legend(loc="best")

    identify = (time >= c["IDENTIFICATION_START"]) & (time <= c["IDENTIFICATION_END"])
    tip_id = detrend(symmetric_tip[identify]); swb_id = detrend(swb1[identify])
    tip_id /= max(np.max(np.abs(tip_id)), 1e-12); swb_id /= max(np.max(np.abs(swb_id)), 1e-12)
    axes[3].plot(time[identify], tip_id, label="symmetric tip deformation normalized")
    axes[3].plot(time[identify], swb_id, label="SWB1 normalized")
    axes[3].set_ylabel("normalized"); axes[3].legend(loc="best")

    axes[4].plot(time, applied["longitudinal_feedback"] * RAD_TO_DEG, label="pitch/q SAS")
    axes[4].plot(time, applied["lift_feedback"] * RAD_TO_DEG, label="altitude/Vz hold")
    axes[4].plot(time, applied["lateral_feedback"] * RAD_TO_DEG, label="roll/p SAS")
    axes[4].plot(time, applied["rap"] * RAD_TO_DEG, label="WF4 rap", alpha=0.8)
    axes[4].set_ylabel("deg"); axes[4].set_xlabel("time [s]"); axes[4].legend(loc="best", ncol=2)

    for axis in axes:
        axis.axvspan(c["SAS_OFF_START"], c["SAS_ON_START"], color="tab:red", alpha=0.10)
        axis.axvspan(c["IDENTIFICATION_START"], c["IDENTIFICATION_END"], color="0.65", alpha=0.18)
        axis.axvline(c["SAS_ON_START"], color="tab:green", ls="--", lw=1)
        axis.grid(alpha=0.22); axis.margins(x=0)
    candidate = summary["matrix_pencil_candidate"]
    short_period = summary["short_period_candidate"]
    fig.suptitle(
        f"NASA-style bounded SAS-off — V={summary['velocity_mps']:.1f} m/s — "
        f"BFF f={candidate.get('frequency_hz', math.nan):.3f} Hz, "
        f"sigma={candidate.get('sigma_per_s', math.nan):+.3f} 1/s — "
        f"rigid-observable f={short_period.get('frequency_hz', math.nan):.3f} Hz"
    )
    fig.savefig(plots / f"{stem}_test_dashboard.png", dpi=190)
    plt.close(fig)

    # Focused view of the actual open/reengage event. This makes a millimetric
    # release transient visible without letting the full 15 s scale hide it.
    detail = (time >= c["SAS_OFF_START"] - 0.20) & (time <= c["SAS_ON_START"] + 0.35)
    fig, axes = plt.subplots(4, 1, figsize=(12, 11), sharex=True, constrained_layout=True)
    axes[0].plot(time[detail], altitude_delta[detail] * 100.0, label="altitude delta")
    axes[0].plot(time[detail], vz[detail] * 100.0, label="vertical speed")
    axes[0].set_ylabel("cm, cm/s"); axes[0].legend(loc="best")
    axes[1].plot(time[detail], pitch[detail], label="pitch")
    axes[1].plot(time[detail], q[detail], label="q")
    axes[1].set_ylabel("deg, deg/s"); axes[1].legend(loc="best")
    axes[2].plot(time[detail], roll[detail], label="roll")
    axes[2].plot(time[detail], p[detail], label="p")
    axes[2].set_ylabel("deg, deg/s"); axes[2].legend(loc="best")
    axes[3].plot(time[detail], applied["longitudinal_feedback"][detail] * RAD_TO_DEG, label="pitch/q SAS")
    axes[3].plot(time[detail], applied["lift_feedback"][detail] * RAD_TO_DEG, label="altitude/Vz hold")
    axes[3].plot(time[detail], applied["lateral_feedback"][detail] * RAD_TO_DEG, label="roll/p SAS")
    axes[3].plot(time[detail], applied["rap"][detail] * RAD_TO_DEG, label="WF4 rap")
    axes[3].set_ylabel("deg"); axes[3].set_xlabel("time [s]"); axes[3].legend(loc="best", ncol=2)
    for axis in axes:
        axis.axvspan(c["SAS_OFF_START"], c["SAS_ON_START"], color="tab:red", alpha=0.10)
        axis.axvspan(c["IDENTIFICATION_START"], c["IDENTIFICATION_END"], color="0.65", alpha=0.18)
        axis.axvline(c["BFF_RAP_START"], color="tab:orange", ls=":", lw=1)
        axis.axvline(c["SAS_ON_START"], color="tab:green", ls="--", lw=1)
        axis.axhline(0.0, color="0.25", lw=0.6); axis.grid(alpha=0.22); axis.margins(x=0)
    fig.suptitle(f"Full-surface hold SAS-off detail — V={summary['velocity_mps']:.1f} m/s")
    fig.savefig(plots / f"{stem}_sas_off_detail.png", dpi=190)
    plt.close(fig)


def write_global_plots(summaries: list[dict], results: Path) -> None:
    valid = [s for s in summaries if s.get("identification_valid")]
    if not valid:
        return
    velocity = np.array([s["velocity_mps"] for s in valid])
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    axes[0].plot(velocity, [s["sigma_swb1_per_s"] for s in valid], "o-", label="MBDyn identified"); axes[0].axhline(0.0, color="k", lw=1); axes[0].set_ylabel("sigma [1/s]")
    axes[1].plot(velocity, [s["frequency_swb1_hz"] for s in valid], "o-", label="BFF tip/SWB1")
    sp_frequency = [s.get("short_period_candidate", {}).get("frequency_hz") if s.get("short_period_candidate", {}).get("accepted") else np.nan for s in valid]
    axes[1].plot(velocity, sp_frequency, "s--", label="rigid q/alpha/pitch candidate")
    axes[1].set_ylabel("frequency [Hz]"); axes[1].legend()
    axes[2].plot(velocity, [s["damping_ratio_swb1"] for s in valid], "o-", label="MBDyn identified"); axes[2].axhline(0.0, color="k", lw=1); axes[2].set_ylabel("damping ratio"); axes[2].set_xlabel("V_INF [m/s]")
    try:
        reference = parse_flutter_point(DEFAULT_F06, point=7)
        reference = [row for row in reference if velocity.min() - 3.0 <= row["velocity_mps"] <= velocity.max() + 3.0]
        if reference:
            ref_v = np.array([row["velocity_mps"] for row in reference])
            ref_sigma = np.array([row["sigma_per_s"] for row in reference])
            ref_f = np.array([row["frequency_hz"] for row in reference])
            ref_zeta = -ref_sigma / np.hypot(ref_sigma, 2.0 * math.pi * ref_f)
            axes[0].plot(ref_v, ref_sigma, "ks--", ms=4, label="NASTRAN SOL 145 point 7")
            axes[1].plot(ref_v, ref_f, "ks--", ms=4, label="NASTRAN SOL 145 point 7")
            axes[2].plot(ref_v, ref_zeta, "ks--", ms=4, label="NASTRAN SOL 145 point 7")
        crossing = zero_crossing(parse_flutter_point(DEFAULT_F06, point=7))
        if crossing:
            for axis in axes:
                axis.axvline(crossing["velocity_mps"], color="tab:red", ls=":", lw=1.2, label="NASTRAN BFF crossing")
    except (OSError, ValueError):
        pass
    for axis in axes:
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            unique = dict(zip(labels, handles))
            axis.legend(unique.values(), unique.keys(), fontsize=8)
    fig.tight_layout(); fig.savefig(results / "local_candidate_overview.png", dpi=170); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--output", type=Path, default=Path(os.environ.get("BFF_OPEN_LOOP_OUTPUT_DIR", DEFAULT_OUTPUT)))
    args = parser.parse_args()
    files = args.files or sorted(args.output.glob("NASA_OL_V_*.nc"))
    if not files:
        raise SystemExit(f"nessun NetCDF in {args.output}")
    summaries = [analyze_case(path, args.output) for path in files]
    write_global_plots(summaries, args.output)
    for value in summaries:
        print(f"{value['case']}: OPEN_LOOP_VALID={value['open_loop_valid']}, ID_VALID={value['identification_valid']}, f={value['frequency_swb1_hz']}, sigma={value['sigma_swb1_per_s']}")


if __name__ == "__main__":
    main()
