#!/usr/bin/env python3
"""Analisi delle run X-56 BFF MBDyn con output NetCDF piatti."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_bff_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from cycler import cycler

# Increase typography sizes and apply thesis color cycle
TITLE_SIZE = 24
LABEL_SIZE = 20
LEGEND_SIZE = 18
TICK_SIZE = 16
THESIS_COLORS = [
    "#6A0DAD", "#008B8B", "#808000", "#C71585", "#4B4B4B",
    "#8B4513", "#2E8B57", "#4682B4", "#B8860B", "#7B68EE",
]
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.titlesize": TITLE_SIZE,
    "figure.titlesize": TITLE_SIZE,
    "axes.labelsize": LABEL_SIZE,
    "legend.fontsize": LEGEND_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "axes.prop_cycle": cycler(color=THESIS_COLORS),
})
from scipy.signal import coherence, detrend, freqz
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path("/mnt/c/Users/Utente/Desktop/bff")
RAW_OUTPUT = Path(os.environ.get("BFF_OUTPUT_DIR", str(DEFAULT_OUTPUT)))
RESULTS = Path(os.environ.get("BFF_RESULTS_DIR", str(RAW_OUTPUT)))
CLEAN_PLOTS = RESULTS / "plots_clean"
RHO_IPS = 9.7284e-8
RHO_SI = 1.039663910516137
NASTRAN_VF = 60.8421
NASTRAN_FF = 2.0597
IN_TO_M = 0.0254
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167
RAD_TO_DEG = 180.0 / math.pi
MASS_IPS = 1.0860117276031622
WEIGHT_LBF = 419.63493154586183
CG_FEM = np.array([163.187383385809, 0.110529571088, 101.239797358848])
GAUSS_W = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
PID_LABELS = [9101, 9102, 9103, 9104, 9201, 9202, 9301, 9302, 9303]


@dataclass
class Pole:
    frequency_hz: float = math.nan
    sigma_per_s: float = math.nan
    damping_ratio: float = math.nan
    cycles: float = 0.0
    frequency_ci_low: float = math.nan
    frequency_ci_high: float = math.nan
    sigma_ci_low: float = math.nan
    sigma_ci_high: float = math.nan
    reliable: bool = False


@dataclass
class Summary:
    velocity_mps: float
    density_kg_m3: float
    dynamic_pressure_pa: float
    mean_vrel_mps: float
    completed: bool
    trim_valid: bool
    identification_valid: bool
    rigid_identification_valid: bool
    safety_triggered: bool
    max_altitude_error_m: float
    max_lateral_drift_m: float
    max_pitch_deg: float
    max_roll_deg: float
    max_yaw_deg: float
    max_surface_deg: float
    saturation_fraction: float
    trim_mean_fx_n: float
    trim_mean_fy_n: float
    trim_mean_fz_n: float
    trim_mean_mx_nm: float
    trim_mean_my_nm: float
    trim_mean_mz_nm: float
    trim_mean_rx_n: float
    short_period_frequency_hz: float
    short_period_sigma_per_s: float
    short_period_damping_ratio: float
    short_period_frequency_ci_low: float
    short_period_frequency_ci_high: float
    short_period_sigma_ci_low: float
    short_period_sigma_ci_high: float
    short_period_reliable: bool
    bff_frequency_hz: float
    bff_sigma_per_s: float
    bff_damping_ratio: float
    bff_frequency_ci_low: float
    bff_frequency_ci_high: float
    bff_sigma_ci_low: float
    bff_sigma_ci_high: float
    bff_reliable: bool
    sas_amplitude_at_nastran_deg: float
    sas_phase_to_q1_deg: float
    sas_to_q1_deg_per_modal_unit: float
    sas_q1_coherence_at_nastran: float
    notch_attenuation_db: float
    bff_classification: str
    rigid_classification: str
    classification: str


def constant(text: str, name: str, default: float = math.nan) -> float:
    pattern = rf"(?m)^\s*set:\s*const\s+(?:real|integer)\s+{re.escape(name)}\s*=\s*([^;]+);"
    match = re.search(pattern, text)
    if not match:
        return default
    expression = match.group(1).strip()
    try:
        return float(expression)
    except ValueError:
        # Le costanti MBDyn usate dal post-processore possono essere scritte
        # anche come, per esempio, ``0.08*deg2rad``.  Il namespace ristretto
        # evita di trasformare questa piccola comodita' in un eval generico.
        try:
            return float(eval(expression, {"__builtins__": {}}, {"deg2rad": math.pi / 180.0, "pi": math.pi}))
        except (NameError, SyntaxError, TypeError, ValueError, ZeroDivisionError):
            return default


def token(value: float) -> str:
    return f"{value:08.4f}".replace(".", "p")


def paths_for(velocity: float) -> tuple[Path, Path]:
    stem = f"V_{token(velocity)}"
    return RAW_OUTPUT / f"{stem}.mbd", RAW_OUTPUT / f"{stem}.nc"


def arr(data: Dataset, name: str) -> np.ndarray:
    return np.asarray(data[name][:]).squeeze()


def filter_first_order(signal: np.ndarray, a1: float, b0: float, b1: float) -> np.ndarray:
    output = np.zeros_like(signal, dtype=float)
    previous_input = 0.0
    for i, value in enumerate(signal):
        previous_output = output[i - 1] if i else 0.0
        output[i] = a1 * previous_output + b0 * value + b1 * previous_input
        previous_input = value
    return output


def smooth_burst(time: np.ndarray, start: float, frequency: float, cycles: float, amplitude: float) -> np.ndarray:
    duration = cycles / frequency
    tau = time - start
    use = (tau >= 0.0) & (tau <= duration)
    result = np.zeros_like(time)
    result[use] = (
        amplitude
        * np.sin(2.0 * np.pi * frequency * tau[use])
        * 0.5
        * (1.0 - np.cos(2.0 * np.pi * tau[use] / duration))
    )
    return result


def aerodynamic_spans() -> dict[int, float]:
    text = (ROOT / "INCLUDE/aerobody.mbd").read_text()
    spans: dict[int, float] = {}
    for block in re.split(r"(?=\s*aerodynamic body:)", text):
        match = re.search(r"aerodynamic body:\s*(\d+)", block)
        if not match:
            continue
        lines = [line.strip() for line in block.splitlines()]
        refs = [i for i, line in enumerate(lines) if line.startswith("reference, node")]
        expression = lines[refs[1] + 1].split(",", 1)[0]
        spans[int(match.group(1))] = float(eval(expression, {"__builtins__": {}}, {}))
    if len(spans) != 58:
        raise RuntimeError(f"attesi 58 elementi aerodinamici (50 ala + 8 winglet), trovati {len(spans)}")
    return spans


def force_residuals(data: Dataset, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = len(position)
    force = np.zeros((n, 3))
    moment = np.zeros((n, 3))
    for label, span in aerodynamic_spans().items():
        for point, weight in enumerate(GAUSS_W):
            scale = 0.5 * span * weight
            f = arr(data, f"elem.aerodynamic.{label}.F_{point}") * scale
            m = arr(data, f"elem.aerodynamic.{label}.M_{point}") * scale
            x = arr(data, f"elem.aerodynamic.{label}.X_{point}")
            force += f
            moment += m + np.cross(x - position, f)
    force[:, 2] -= WEIGHT_LBF
    return force, moment


def matrix_pencil(time: np.ndarray, signals: np.ndarray, rank_max: int = 16) -> tuple[np.ndarray, np.ndarray]:
    n = len(time)
    if n < 100:
        return np.array([], dtype=complex), np.empty((signals.shape[0], 0), complex)
    dt = float(np.median(np.diff(time)))
    clean = []
    scales = []
    for signal in signals:
        value = detrend(np.asarray(signal, float))
        scale = np.sqrt(np.mean(value**2))
        if not np.isfinite(scale) or scale < 1e-12:
            scale = 1.0
        clean.append(value / scale)
        scales.append(scale)
    clean = np.asarray(clean)
    length = min(max(30, n // 3), 180)
    blocks0, blocks1 = [], []
    for signal in clean:
        windows = np.lib.stride_tricks.sliding_window_view(signal, length + 1)
        blocks0.append(windows[:, :-1].T)
        blocks1.append(windows[:, 1:].T)
    h0 = np.vstack(blocks0)
    h1 = np.vstack(blocks1)
    u, singular, vh = np.linalg.svd(h0, full_matrices=False)
    rank = int(min(rank_max, max(4, np.count_nonzero(singular > singular[0] * 2e-3))))
    operator = u[:, :rank].conj().T @ h1 @ vh[:rank].conj().T @ np.diag(1.0 / singular[:rank])
    eigenvalues = np.linalg.eigvals(operator)
    poles = np.log(eigenvalues) / dt
    vandermonde = eigenvalues[None, :] ** np.arange(n)[:, None]
    amplitudes = np.vstack([np.linalg.lstsq(vandermonde, signal, rcond=None)[0] for signal in clean])
    return poles, amplitudes


def select_pole(
    time: np.ndarray,
    signals: np.ndarray,
    low: float,
    high: float,
    elastic: bool,
    target_frequency: float | None = None,
) -> Pole:
    poles, amplitudes = matrix_pencil(time, signals)
    if not len(poles):
        return Pole()
    frequency = np.abs(np.imag(poles)) / (2.0 * np.pi)
    candidates = np.flatnonzero((np.imag(poles) > 0.0) & (frequency >= low) & (frequency <= high))
    if not len(candidates):
        return Pole()
    if elastic:
        score = np.abs(amplitudes[0, candidates]) + 0.5 * np.abs(amplitudes[2, candidates])
    else:
        score = np.abs(amplitudes[1, candidates])
    # Il BFF viene seguito per continuita' a partire dal caso valido a
    # velocita' maggiore. Questo impedisce salti verso armoniche 2.4--2.8 Hz
    # quando il ramo elastico e' ancora poco osservabile.
    if elastic and target_frequency is not None and np.isfinite(target_frequency):
        distance = np.abs(frequency[candidates] - target_frequency)
        keep = distance <= 0.45
        if not np.any(keep):
            return Pole()
        candidates = candidates[keep]
        distance = distance[keep]
        score = score[keep] * np.exp(-0.5 * (distance / 0.20) ** 2)
    chosen = candidates[int(np.argmax(score))]
    pole = poles[chosen]
    f = float(frequency[chosen])
    sigma = float(np.real(pole))
    duration = float(time[-1] - time[0])
    estimates = []
    for cut_start in (0.0, 0.10, 0.20):
        for cut_end in (0.0, 0.10, 0.20):
            use = (time >= time[0] + cut_start) & (time <= time[-1] - cut_end)
            p2, a2 = matrix_pencil(time[use], signals[:, use], rank_max=14)
            if not len(p2):
                continue
            f2 = np.abs(np.imag(p2)) / (2.0 * np.pi)
            can = np.flatnonzero((np.imag(p2) > 0.0) & (f2 >= low) & (f2 <= high) & (np.abs(f2 - f) < 0.35))
            if len(can):
                j = can[np.argmin(np.abs(f2[can] - f))]
                estimates.append((float(f2[j]), float(np.real(p2[j]))))
    if estimates:
        values = np.asarray(estimates)
        f_lo, f_hi = np.percentile(values[:, 0], [5, 95])
        s_lo, s_hi = np.percentile(values[:, 1], [5, 95])
    else:
        f_lo = f_hi = s_lo = s_hi = math.nan
    cycles = f * duration
    zeta = -sigma / math.sqrt(sigma**2 + (2.0 * math.pi * f) ** 2)
    # Un modo monotono non passa: servono almeno tre cicli e stabilita' del fit.
    reliable = bool(
        cycles >= 3.0
        and len(estimates) >= 5
        and np.isfinite(f_lo)
        and f_hi - f_lo < 0.35
        and abs(s_hi - s_lo) < 4.0
        and f_lo - 0.05 <= f <= f_hi + 0.05
        and s_lo - 0.05 <= sigma <= s_hi + 0.05
    )
    return Pole(f, sigma, zeta, cycles, f_lo, f_hi, s_lo, s_hi, reliable)


def harmonic(signal: np.ndarray, time: np.ndarray, frequency: float) -> complex:
    clean = detrend(np.asarray(signal, float))
    window = np.hanning(len(clean))
    return 2.0 * np.sum(clean * window * np.exp(-2j * np.pi * frequency * time)) / np.sum(window)


def scalar_at(frequency: np.ndarray, value: np.ndarray, target: float) -> float:
    return float(value[np.argmin(np.abs(frequency - target))])


def analyze_case(nc_path: Path, bff_target: float | None = None) -> Summary:
    stem = nc_path.stem
    input_path = nc_path.with_suffix(".mbd")
    if not input_path.exists():
        raise FileNotFoundError(f"manca l'input associato {input_path}")
    text = input_path.read_text()
    velocity_inf = constant(text, "V_INF")
    final_time = constant(text, "FINAL_TIME")
    excite = constant(text, "EXCITATION_START")
    excite_f = constant(text, "EXCITATION_FREQUENCY")
    excite_cycles = constant(text, "EXCITATION_CYCLES")
    identify = excite + excite_cycles / excite_f + 0.30
    safety_arm = constant(text, "SAFETY_ARM_TIME")
    if not np.isfinite(safety_arm):
        safety_arm = excite + excite_cycles / excite_f + 4.0
    trim_scale = (1.146e-7 / RHO_IPS) * (63.0 / velocity_inf) ** 2
    trim_surface = (0.11983843 * trim_scale + 0.66144968) / RAD_TO_DEG
    # Compatibilita' con le run precedenti: la correzione body-flap e'
    # applicata soltanto se la relativa costante compare nell'input salvato.
    body_trim_correction = 0.8400 * (trim_scale - 1.2) if "BODY_TRIM_CORRECTION" in text else 0.0
    body_trim_surface = trim_surface + (-0.37830 * trim_scale + body_trim_correction) / RAD_TO_DEG
    with Dataset(nc_path) as data:
        time = arr(data, "time")
        if len(time) < 3:
            raise RuntimeError(f"{nc_path}: nessun passo temporale utilizzabile")
        position = arr(data, "node.struct.990000.X")
        phi = arr(data, "node.struct.990000.Phi")
        velocity = arr(data, "node.struct.990000.XP")
        omega = arr(data, "node.struct.990000.Omega")
        rotation = Rotation.from_rotvec(phi).as_matrix()
        euler = np.unwrap(Rotation.from_rotvec(phi).as_euler("xyz"), axis=0)
        pid = {label: arr(data, f"elem.loadable.{label}.output") for label in PID_LABELS}
        modal_q = arr(data, "elem.joint.5.a")
        modal_qd = arr(data, "elem.joint.5.aPrime")
        tip_l_velocity = arr(data, "node.struct.990020.XP") * IN_TO_M
        tip_r_velocity = arr(data, "node.struct.991020.XP") * IN_TO_M
        aero_labels = sorted(aerodynamic_spans())
        aero_alpha = np.mean(
            np.column_stack([arr(data, f"elem.aerodynamic.{label}.alpha_1") for label in aero_labels]), axis=1
        )
        aero_beta = np.mean(
            np.column_stack([arr(data, f"elem.aerodynamic.{label}.gamma_1") for label in aero_labels]), axis=1
        )
        # Controlli prima e dopo la dinamica attuatore.
        a1 = constant(text, "ACT_A1")
        b0 = constant(text, "ACT_B0")
        b1 = constant(text, "ACT_B1")
        limit = constant(text, "SURFACE_CORRECTION_LIMIT")
        # Stessa deadband del safety MBDyn. I contributi filtrati separatamente
        # sono salvati per audit; la superficie reale usa la somma saturata e
        # filtrata, quindi resta esattamente quella simulata.
        q_deadband = np.maximum(0.0, omega[:, 1] - 18.0 / RAD_TO_DEG) + np.minimum(
            0.0, omega[:, 1] + 18.0 / RAD_TO_DEG
        )
        pitch_deadband = np.maximum(0.0, euler[:, 1] - 12.0 / RAD_TO_DEG) + np.minimum(
            0.0, euler[:, 1] + 12.0 / RAD_TO_DEG
        )
        safety_raw = np.where(time >= safety_arm, np.clip(-0.35 * q_deadband - 0.25 * pitch_deadband,
                                                          -4.0 / RAD_TO_DEG, 4.0 / RAD_TO_DEG), 0.0)
        long_raw = pid[9103] + pid[9104] + safety_raw
        lat_raw = pid[9201] + pid[9202]
        direction_raw = pid[9301] + pid[9302] + pid[9303]
        long_act = filter_first_order(np.clip(long_raw, -limit, limit), a1, b0, b1)
        lat_act = filter_first_order(np.clip(lat_raw, -limit, limit), a1, b0, b1)
        direction_act = filter_first_order(np.clip(direction_raw, -limit, limit), a1, b0, b1)
        long_baseline = filter_first_order(pid[9104], a1, b0, b1)
        long_outer = filter_first_order(pid[9103], a1, b0, b1)
        lat_baseline = filter_first_order(pid[9202], a1, b0, b1)
        lat_outer = filter_first_order(pid[9201], a1, b0, b1)
        dir_baseline = filter_first_order(pid[9302], a1, b0, b1)
        dir_outer = filter_first_order(pid[9301], a1, b0, b1)
        dir_vy_loop = filter_first_order(pid[9303], a1, b0, b1)
        safety_act = filter_first_order(safety_raw, a1, b0, b1)
        burst = smooth_burst(time, excite, excite_f, excite_cycles, constant(text, "EXCITATION_AMPLITUDE"))
        surfaces = {
            "BFL": body_trim_surface + long_act + direction_act,
            "BFR": body_trim_surface + long_act - direction_act,
            "WF1L": trim_surface + lat_act,
            "WF1R": trim_surface - lat_act,
            "WF2L": trim_surface + lat_act,
            "WF2R": trim_surface - lat_act,
            "WF3L": np.full_like(time, trim_surface),
            "WF3R": np.full_like(time, trim_surface),
            "WF4L": trim_surface + burst,
            "WF4R": trim_surface + burst,
        }
        force_lbf, moment_lbfin = force_residuals(data, position)
        constraint_reaction_lbf = arr(data, "elem.joint.1.F")

    velocity_si = velocity * IN_TO_M
    flow_global = np.column_stack((np.full(len(time), velocity_inf), np.zeros((len(time), 2)))) - velocity_si
    flow_body = np.einsum("nji,nj->ni", rotation, flow_global)
    vrel = np.linalg.norm(flow_body, axis=1)
    alpha = np.arctan2(flow_body[:, 2], flow_body[:, 0]) * RAD_TO_DEG
    beta = np.arctan2(flow_body[:, 1], np.hypot(flow_body[:, 0], flow_body[:, 2])) * RAD_TO_DEG
    dynamic_pressure = 0.5 * RHO_SI * vrel**2
    euler_deg = euler * RAD_TO_DEG
    omega_deg = omega * RAD_TO_DEG
    altitude = (position[:, 2] - position[0, 2]) * IN_TO_M
    lateral = (position[:, 1] - position[0, 1]) * IN_TO_M
    vertical_speed = velocity_si[:, 2]
    tip_l_acc = np.gradient(tip_l_velocity[:, 2], time)
    tip_r_acc = np.gradient(tip_r_velocity[:, 2], time)
    tip_symmetric = 0.5 * (tip_l_acc + tip_r_acc)

    saturation = (
        (np.abs(long_raw) >= limit)
        | (np.abs(lat_raw) >= limit)
        | (np.abs(direction_raw) >= limit)
        | np.logical_or.reduce([np.abs(pid[label]) >= (5.99 / RAD_TO_DEG) for label in (9103, 9104, 9201, 9202)])
    )
    safety = (time >= safety_arm) & (
        (np.abs(euler_deg[:, 0]) > 20.0)
        | (np.abs(euler_deg[:, 1]) > 20.0)
        | (np.abs(euler_deg[:, 2]) > 30.0)
        | (np.abs(altitude) > 20.0)
        | (np.abs(omega_deg[:, 1]) > 18.0)
        | (np.abs(modal_q[:, 0]) > 1.5)
        | (np.maximum(np.abs(tip_l_acc), np.abs(tip_r_acc)) > 20.0)
        | (vrel < 0.5 * velocity_inf)
        | (vrel > 1.5 * velocity_inf)
    )
    safety_indices = np.flatnonzero(safety)
    safe_end = float(time[safety_indices[0]]) if len(safety_indices) else float(time[-1])
    id_use = (time >= identify) & (time <= safe_end)
    # Due secondi interamente pre-burst: la media copre almeno un periodo del
    # modo rigido lento e non dipende piu' dalla sua fase istantanea.
    trim_use = (time >= max(0.0, excite - 2.02)) & (time <= excite - 0.02)
    completed = bool(time[-1] >= final_time - 0.03)
    trim_valid = bool(
        np.count_nonzero(trim_use) > 20
        and np.max(np.abs(euler_deg[trim_use, 0])) < 3.0
        # Il volo orizzontale richiede flight-path angle circa nullo, non
        # pitch nullo: a 30 m/s l'X-56 necessita fisiologicamente 4--5 deg di
        # incidenza. L'inviluppo viene controllato anche tramite alpha.
        and np.max(np.abs(euler_deg[trim_use, 1])) < 8.0
        and np.max(np.abs(alpha[trim_use])) < 12.0
        and np.max(np.abs(euler_deg[trim_use, 2])) < 5.0
        and np.max(np.abs(vertical_speed[trim_use])) < 2.0
        and np.max(np.abs(vrel[trim_use] - velocity_inf)) < 0.10 * velocity_inf
        and np.max(np.abs(altitude[trim_use])) < 5.0
        and np.max(np.abs(lateral[trim_use])) < 1.0e-6
        and not np.any(saturation[trim_use])
        # Fx e Fy sono bilanciate dalle reazioni ammesse e non entrano nel trim.
        and abs(np.mean(force_lbf[trim_use, 2]) * LBF_TO_N) < 0.05 * WEIGHT_LBF * LBF_TO_N
        and np.all(np.abs(np.mean(moment_lbfin[trim_use], axis=0) * LBFIN_TO_NM) < 75.0)
        and np.max(np.abs(position[trim_use, 0] - position[0, 0])) < 1.0e-6
        and np.max(np.abs(position[trim_use, 1] - position[0, 1])) < 1.0e-6
        and np.max(np.abs(constraint_reaction_lbf[trim_use, 2])) < 1.0e-8
    )
    if np.count_nonzero(id_use) >= 100:
        identification_signals = np.vstack((modal_q[id_use, 0], euler[id_use, 1], tip_symmetric[id_use]))
        bff = select_pole(time[id_use], identification_signals, 1.0, 4.0, True, bff_target)
        short = select_pole(time[id_use], identification_signals, 0.10, 1.0, False)
        h_q1 = harmonic(modal_q[id_use, 0], time[id_use], NASTRAN_FF)
        h_sas = harmonic(long_act[id_use], time[id_use], NASTRAN_FF)
        sas_amp = abs(h_sas) * RAD_TO_DEG
        sas_phase = float(np.angle(h_sas / h_q1, deg=True)) if abs(h_q1) > 0 else math.nan
        sas_ratio = float(abs(h_sas / h_q1) * RAD_TO_DEG) if abs(h_q1) > 0 else math.nan
        fs_id = 1.0 / float(np.median(np.diff(time[id_use])))
        coh_f, coh = coherence(detrend(long_act[id_use]), detrend(modal_q[id_use, 0]), fs=fs_id,
                               nperseg=min(128, np.count_nonzero(id_use)))
        sas_coherence = scalar_at(coh_f, coh, NASTRAN_FF)
    else:
        bff = Pole()
        short = Pole()
        sas_amp = sas_phase = sas_ratio = sas_coherence = math.nan
    notch_b = [constant(text, "NOTCH_B0"), constant(text, "NOTCH_B1"), constant(text, "NOTCH_B2")]
    notch_a = [1.0, -constant(text, "NOTCH_A1"), -constant(text, "NOTCH_A2")]
    w, response = freqz(notch_b, notch_a, worN=8192, fs=50.0)
    notch_db = 20.0 * math.log10(max(abs(response[np.argmin(abs(w - NASTRAN_FF))]), 1e-12))
    # Una soglia di sicurezza invalida solo i campioni successivi. Se prima
    # dell'evento rimangono almeno tre cicli e il polo e' robusto alle finestre,
    # l'identificazione resta utilizzabile (ma il caso e' segnalato come
    # troncato). Una run MBDyn incompleta, invece, non e' mai valida.
    usable_cycles = bff.frequency_hz * max(0.0, safe_end - identify)
    rigid_usable_cycles = short.frequency_hz * max(0.0, safe_end - identify)
    id_saturation_fraction = float(np.mean(saturation[id_use])) if np.any(id_use) else 1.0
    identification_valid = bool(
        completed and trim_valid and bff.reliable and usable_cycles >= 3.0 and id_saturation_fraction < 0.01
    )
    rigid_identification_valid = bool(
        completed and trim_valid and short.reliable and rigid_usable_cycles >= 3.0
        and id_saturation_fraction < 0.01
    )

    def pole_state(pole: Pole) -> str:
        if not pole.reliable:
            return "insufficient evidence"
        if pole.sigma_ci_low > 0.0:
            return "unstable"
        if pole.sigma_ci_high < 0.0:
            return "stable"
        return "neutral/CI crosses zero"

    invalid_reason = ""
    if not completed:
        invalid_reason = "invalid: incomplete run"
    elif not trim_valid:
        invalid_reason = "invalid: trim/flight envelope"
    elif id_saturation_fraction >= 0.01:
        invalid_reason = "invalid: actuator saturation in identification window"

    if invalid_reason:
        bff_classification = rigid_classification = classification = invalid_reason
    else:
        bff_state = pole_state(bff)
        rigid_state = pole_state(short)
        bff_classification = f"BFF {bff_state}"
        rigid_classification = f"rigid longitudinal {rigid_state}"
        unstable_modes = []
        if rigid_state == "unstable":
            unstable_modes.append("rigid longitudinal")
        if bff_state == "unstable":
            unstable_modes.append("BFF")
        if unstable_modes:
            classification = "unstable: " + " + ".join(unstable_modes)
        elif rigid_state == "stable" and bff_state == "stable":
            classification = "stable: rigid longitudinal + BFF"
        elif "neutral/CI crosses zero" in (rigid_state, bff_state):
            classification = "neutral: at least one CI crosses zero"
        else:
            classification = "inconclusive: " + rigid_classification + "; " + bff_classification

    # Dataset derivato completo, in SI salvo gli angoli esplicitamente in deg.
    columns: dict[str, np.ndarray] = {
        "time_s": time,
        "V_inf_mps": np.full_like(time, velocity_inf),
        "q_Pa": dynamic_pressure,
        "Vrel_mps": vrel,
        "X_m": position[:, 0] * IN_TO_M,
        "Y_m": position[:, 1] * IN_TO_M,
        "Z_m": position[:, 2] * IN_TO_M,
        "altitude_delta_m": altitude,
        "lateral_delta_m": lateral,
        "roll_deg": euler_deg[:, 0], "pitch_deg": euler_deg[:, 1], "yaw_deg": euler_deg[:, 2],
        "p_deg_s": omega_deg[:, 0], "q_deg_s": omega_deg[:, 1], "r_deg_s": omega_deg[:, 2],
        "alpha_deg": alpha, "beta_deg": beta, "alpha_aero_mean_deg": aero_alpha, "beta_aero_mean_deg": aero_beta,
        "vertical_speed_mps": vertical_speed,
        "trim_surface_deg": np.full_like(time, trim_surface * RAD_TO_DEG),
        "long_sas_deg": long_act * RAD_TO_DEG, "lat_sas_deg": lat_act * RAD_TO_DEG,
        "dir_sas_deg": direction_act * RAD_TO_DEG, "burst_deg": burst * RAD_TO_DEG,
        "long_baseline_sas_deg": long_baseline * RAD_TO_DEG,
        "long_outer_loop_deg": long_outer * RAD_TO_DEG,
        "lat_baseline_sas_deg": lat_baseline * RAD_TO_DEG,
        "lat_outer_loop_deg": lat_outer * RAD_TO_DEG,
        "dir_baseline_sas_deg": dir_baseline * RAD_TO_DEG,
        "dir_outer_loop_deg": dir_outer * RAD_TO_DEG,
        "dir_vy_loop_deg": dir_vy_loop * RAD_TO_DEG,
        "flutter_safety_deg": safety_act * RAD_TO_DEG,
        "constraint_Rx_N": constraint_reaction_lbf[:, 0] * LBF_TO_N,
        "constraint_Ry_N": constraint_reaction_lbf[:, 1] * LBF_TO_N,
        "constraint_Rz_N": constraint_reaction_lbf[:, 2] * LBF_TO_N,
        "tip_left_accel_z_mps2": tip_l_acc, "tip_right_accel_z_mps2": tip_r_acc,
        "Fx_residual_N": force_lbf[:, 0] * LBF_TO_N, "Fy_residual_N": force_lbf[:, 1] * LBF_TO_N,
        "Fz_residual_N": force_lbf[:, 2] * LBF_TO_N, "Mx_residual_Nm": moment_lbfin[:, 0] * LBFIN_TO_NM,
        "My_residual_Nm": moment_lbfin[:, 1] * LBFIN_TO_NM, "Mz_residual_Nm": moment_lbfin[:, 2] * LBFIN_TO_NM,
        "saturation": saturation.astype(float), "safety": safety.astype(float),
    }
    for name, signal in surfaces.items():
        columns[f"surface_{name}_deg"] = signal * RAD_TO_DEG
    for index in range(modal_q.shape[1]):
        columns[f"modal_q_FEM_{index + 7}"] = modal_q[:, index]
        columns[f"modal_qdot_FEM_{index + 7}"] = modal_qd[:, index]
    for label, name in zip(PID_LABELS, ("alt", "vz", "pitch", "q", "roll", "p", "yaw", "r", "vy")):
        columns[f"pid_{name}_output"] = pid[label]
    with (RESULTS / f"{stem}_timeseries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))

    make_clean_plots(
        stem, velocity_inf, 0.5 * RHO_SI * velocity_inf**2,
        time, euler_deg, omega_deg, altitude, surfaces, modal_q,
        excite, excite + excite_cycles / excite_f, identify,
        bff, completed, trim_valid, id_saturation_fraction,
    )
    result = Summary(
        velocity_inf, RHO_SI, 0.5 * RHO_SI * velocity_inf**2, float(np.mean(vrel[trim_use])), completed,
        trim_valid, identification_valid, rigid_identification_valid,
        bool(len(safety_indices)), float(np.max(np.abs(altitude))),
        float(np.max(np.abs(lateral))),
        float(np.max(np.abs(euler_deg[:, 1]))), float(np.max(np.abs(euler_deg[:, 0]))),
        float(np.max(np.abs(euler_deg[:, 2]))), float(max(np.max(np.abs(v)) for v in surfaces.values()) * RAD_TO_DEG),
        id_saturation_fraction,
        *[float(value) for value in (np.mean(force_lbf[trim_use], axis=0) * LBF_TO_N)],
        *[float(value) for value in (np.mean(moment_lbfin[trim_use], axis=0) * LBFIN_TO_NM)],
        float(np.mean(constraint_reaction_lbf[trim_use, 0]) * LBF_TO_N),
        short.frequency_hz, short.sigma_per_s, short.damping_ratio,
        short.frequency_ci_low, short.frequency_ci_high, short.sigma_ci_low, short.sigma_ci_high,
        short.reliable,
        bff.frequency_hz, bff.sigma_per_s, bff.damping_ratio, bff.frequency_ci_low, bff.frequency_ci_high,
        bff.sigma_ci_low, bff.sigma_ci_high, bff.reliable, sas_amp, sas_phase, sas_ratio, sas_coherence,
        notch_db, bff_classification, rigid_classification, classification,
    )
    return result


def make_clean_plots(
    stem: str,
    velocity: float,
    q_nominal: float,
    time: np.ndarray,
    euler: np.ndarray,
    omega: np.ndarray,
    altitude: np.ndarray,
    surfaces: dict[str, np.ndarray],
    modal: np.ndarray,
    excite_start: float,
    excite_end: float,
    identify_start: float,
    bff: Pole,
    completed: bool,
    trim_valid: bool,
    saturation_fraction: float,
) -> None:
    """Due sole figure leggibili per caso, raccolte in ``plots_clean``."""
    CLEAN_PLOTS.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True, constrained_layout=True)
    axes[0].plot(time, altitude, lw=1.6)
    axes[0].axhline(0.0, color="0.25", lw=0.7)
    axes[0].set_ylabel("altitude delta [m]")

    axes[1].plot(time, euler[:, 1], lw=1.5)
    axes[1].axhline(0.0, color="0.25", lw=0.7)
    axes[1].set_ylabel("pitch [deg]")

    for index, label in ((0, "p"), (1, "q"), (2, "r")):
        axes[2].plot(time, omega[:, index], label=label, lw=1.25)
    axes[2].set_ylabel("rates [deg/s]")
    axes[2].legend(ncol=3, loc="upper right")

    surface_series = (
        ("BFL", "body flap L", "-"),
        ("BFR", "body flap R", "--"),
        ("WF1L", "WF1/WF2 L", "-"),
        ("WF1R", "WF1/WF2 R", "--"),
        ("WF4L", "WF4 L=R (burst)", "-"),
    )
    for name, label, style in surface_series:
        series = surfaces[name] * RAD_TO_DEG
        series = series - series[0]
        axes[3].plot(time, series, label=label, ls=style, lw=1.15)
    axes[3].set_ylabel(r"\Delta\delta_control_surfaces")
    axes[3].set_xlabel("time [s]")
    axes[3].legend(ncol=3, loc="upper right", fontsize=8)

    for axis in axes:
        axis.axvspan(excite_start, min(excite_end, time[-1]), color="0.85", alpha=0.45, lw=0)
        axis.grid(alpha=0.22)
        axis.margins(x=0)
    fig.suptitle(f"V={velocity:.2f} m/s   q={q_nominal:.1f} Pa")
    fig.savefig(CLEAN_PLOTS / f"{stem}_flight.png", dpi=160)
    plt.close(fig)

    post = time >= identify_start
    if np.count_nonzero(post) < 20:
        return
    ti = time[post]
    mode7 = detrend(modal[post, 0])
    pitch = detrend(euler[post, 1])

    identification_valid = (
        completed
        and trim_valid
        and saturation_fraction < 0.01
        and bff.reliable
        and np.isfinite(bff.frequency_hz)
        and np.isfinite(bff.damping_ratio)
    )
    if identification_valid:
        details = f"f = {bff.frequency_hz:.3f} Hz\ndamping ζ = {100.0 * bff.damping_ratio:+.2f} %"
    else:
        details = "f = n.d.\ndamping ζ = n.d."

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.0), sharex=True, constrained_layout=True)
    axes[0].plot(ti, mode7, lw=1.45)
    axes[0].axhline(0.0, color="0.25", lw=0.7)
    axes[0].text(
        0.02, 0.96, details, transform=axes[0].transAxes,
        va="top", ha="left", color="black",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.45", "alpha": 0.90},
    )
    axes[0].set_ylabel("Mode 7")

    axes[1].plot(ti, pitch, lw=1.25)
    axes[1].axhline(0.0, color="0.25", lw=0.7)
    axes[1].set_ylabel("Δpitch [deg]")
    axes[1].set_xlabel("time post-burst [s]")

    for axis in axes:
        axis.grid(alpha=0.22)
        axis.margins(x=0)
    fig.suptitle(f"Mode 7 and pitch — V={velocity:.2f} m/s, q={q_nominal:.1f} Pa")
    fig.savefig(CLEAN_PLOTS / f"{stem}_mode7.png", dpi=160)
    plt.close(fig)


def write_study(summaries: list[Summary]) -> None:
    summaries.sort(key=lambda value: value.velocity_mps)
    RESULTS.mkdir(exist_ok=True)
    CLEAN_PLOTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(summaries[0])))
        writer.writeheader(); writer.writerows(asdict(value) for value in summaries)
    q = np.asarray([value.dynamic_pressure_pa for value in summaries])
    f = np.asarray([value.bff_frequency_hz for value in summaries])
    sigma = np.asarray([value.bff_sigma_per_s for value in summaries])
    bff_valid = np.asarray([value.identification_valid for value in summaries])
    rigid_valid = np.asarray([value.rigid_identification_valid for value in summaries])
    rigid_sigma = np.asarray([value.short_period_sigma_per_s for value in summaries])
    fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    ax[0].plot(q[bff_valid], f[bff_valid], "o-", label="BFF branch valid")
    ax[0].plot(q[~bff_valid], f[~bff_valid], "x", color="0.55", label="case excluded")
    ax[0].axhline(NASTRAN_FF, color="k", ls="--", label="NASTRAN"); ax[0].set_ylabel("BFF frequency [Hz]"); ax[0].legend()
    ax[1].plot(q[bff_valid], sigma[bff_valid], "o-", label="BFF")
    ax[1].plot(q[~bff_valid], sigma[~bff_valid], "x", color="0.55")
    ax[1].axhline(0, color="k", lw=.8); ax[1].set_ylabel("BFF sigma [1/s]")
    ax[2].plot(q[rigid_valid], rigid_sigma[rigid_valid], "o-", label="rigid longitudinal")
    ax[2].plot(q[~rigid_valid], rigid_sigma[~rigid_valid], "x", color="0.55")
    ax[2].axhline(0, color="k", lw=.8); ax[2].set_ylabel("rigid sigma [1/s]")
    ax[2].set_xlabel("nominal q [Pa]")
    for axis in ax: axis.grid(alpha=.25)
    fig.savefig(CLEAN_PLOTS / "BFF_summary.png", dpi=160); plt.close(fig)
    report = [
        "# X-56 body-freedom flutter — studio MBDyn 4-DOF rigidi", "",
        f"Densità fissa: `{RHO_SI:.9f} kg/m³` (`{RHO_IPS:.7g}` IPS). Variabile esterna unica: `V_INF`.",
        "Il riferimento NASTRAN separato è `Vf=60.8421 m/s, f=2.0597 Hz` (DLM/PK open-loop); MBDyn è strip-C81 quasi-stazionario, non lineare e closed-loop.",
        "I modi FEM 1–6 sono rigidi e sono esclusi perché il moto rigido è rappresentato dal floating frame MBDyn; il joint esterno rimuove poi soltanto X e Y. I modi flessibili 7–12 coprono 3.217–12.759 Hz.",
        "Il solo joint esterno vincola le traslazioni globali X e Y al CG. Z, roll, pitch e yaw restano liberi; non esistono spinta, throttle o airspeed hold. Rx e Ry sono reazioni ideali e non segnali di controllo.", "",
        "| V [m/s] | q [Pa] | trim | max abs(Z-Z0) [m] | rigido f [Hz] | sigma rigido [1/s] | BFF f [Hz] | sigma BFF [1/s] | Rx trim [N] | sat. | esito globale |",
        "|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for s in summaries:
        report.append(f"| {s.velocity_mps:.4f} | {s.dynamic_pressure_pa:.1f} | {'OK' if s.trim_valid else 'NO'} | {s.max_altitude_error_m:.3f} | {s.short_period_frequency_hz:.3f} | {s.short_period_sigma_per_s:+.3f} | {s.bff_frequency_hz:.3f} | {s.bff_sigma_per_s:+.3f} | {s.trim_mean_rx_n:.1f} | {100*s.saturation_fraction:.1f}% | {s.classification} |")
    bracket = None
    valid = [s for s in summaries if s.identification_valid]
    for left, right in zip(valid, valid[1:]):
        if left.bff_sigma_per_s * right.bff_sigma_per_s <= 0.0:
            bracket = (left, right)
            break
    if bracket:
        report += ["", f"Frontiera del ramo BFF closed-loop racchiusa tra **{bracket[0].velocity_mps:.4f} e {bracket[1].velocity_mps:.4f} m/s** (q={bracket[0].dynamic_pressure_pa:.1f}–{bracket[1].dynamic_pressure_pa:.1f} Pa). La stabilita' globale richiede anche sigma rigido < 0."]
    report += [
        "", "## Controllo e filtri", "",
        "Le superfici ricevono `delta = delta_trim + delta_SAS + delta_outer + delta_burst + delta_safety`; il burst e' l'unico termine rimosso. Tutti i comandi passano in un attuatore Tustin del primo ordine con tau=0.01 s e limite di correzione ±8 deg.", "",
        "- scheduling longitudinale: S=1.5 q60/(q_inf+0.5 q60), quindi S=2.0 a 30 m/s, 1.0 a 60 m/s e 0.81 a 70 m/s;",
        "- quota: Kp=4.00e-4 S rad/in, Ki=1.00e-5 S rad/(in s), LP1 a 0.15 Hz; Vz: Kp=1.20e-3 S rad/(in/s), LP1 a 0.35 Hz; entrambi limitati a ±2 deg;",
        "- pitch: Kp=-0.50 S, Ki=-0.008 S; pitch rate: Kp=-0.60 S s; limiti individuali ±6 deg;",
        "- roll: Kp=-0.80, Ki=-0.040; roll rate: Kp=-0.80 s; limiti ±6 deg;",
        "- yaw: Kp=+1.80, Ki=+0.150; yaw rate: Kp=+3.00 s; Vy: Kp=6.0e-4 rad/(in/s), Ki=2.0e-5;",
        "- notch digitale (fs=50 Hz): f0=2.0597 Hz, Q=3, B=[0.95862109,-1.85337926,0.95862109], A=[1,-1.85337926,0.91724218];",
        "- i low-pass longitudinali sono del primo ordine per ridurre il ritardo di fase; il precedente Butterworth a 0.8 Hz resta soltanto sul canale laterale Vy;",
        "- attuatore Tustin: y[k]=0.5 u[k]+0.5 u[k-1];",
        "- safety longitudinale: armato 4 s dopo il burst, deadband |q|=18 deg/s e |pitch|=12 deg, guadagni -0.35 s e -0.25, limite ±4 deg.",
        "", "## Criteri automatici", "",
        "Trim valido: X e Y costanti entro 1e-6 in; Rz=0 entro 1e-8 lbf; |roll|<3 deg, |pitch|<8 deg, |alpha|<12 deg, |yaw|<5 deg, |Vz|<2 m/s, errore Vrel<10%, quota<5 m, nessuna saturazione, |Fz| medio <5% del peso e ogni momento medio <75 N m. Il limite pitch ammette l'incidenza fisica necessaria alle basse velocita'; Fx e Fy sono equilibrate esclusivamente da Rx e Ry.",
        "Identificazione valida: run completa, trim valido, almeno 3 cicli, matrix-pencil multi-segnale robusto a nove perturbazioni della finestra, polo nominale interno alla propria CI e saturazione <1%.",
        "Il BFF e' seguito per continuita' in frequenza procedendo dalle velocita' maggiori verso le minori. La classificazione globale e' `stable` soltanto quando sia il modo longitudinale rigido sia il BFF sono identificati e hanno tutta la CI di sigma sotto zero.",
    ]
    (RESULTS / "REPORT.md").write_text("\n".join(report) + "\n")
    audit = """# Audit del modello

- Nodo modale: `990000`; giunto modale: `5`. Il total-pin `1` ha due componenti attive: traslazioni globali X e Y.
- Gradi rigidi liberi: Z, roll, pitch e yaw. X e Y sono le sole componenti traslazionali `active`; Z e tutte le rotazioni del joint `1` sono `inactive`.
- Modi FEM attivi: 7–12. Modi 1–6: rigidi esclusi. Frequenze secche: 3.2171, 5.3027, 8.7051, 11.1640, 12.2571, 12.7589 Hz.
- Superfici: BFL/BFR `1004/2004`; WF1 `1008/2008`; WF2 `1011/2011`; WF3 `1014/2014`; WF4 `1017/2017`.
- Combinazione simmetrica: stesso segno numerico L/R. Combinazione differenziale: segno opposto L/R.
- Rappresentazione del moto: velocità strutturale iniziale nulla e vento uniforme `+VINF`; `Vrel` viene sempre ricalcolata come vento meno velocità inerziale.
- Controllo: quota/Vz → riferimento pitch; pitch/q → body flap simmetrici; roll/p → WF1/WF2 differenziali; yaw/r/Vy → body flap differenziali.
- Non esistono propulsione, throttle, airspeed hold o feedback delle reazioni X/Y.
- Il trim e' un termine separato; le correzioni degli attuatori partono da zero e non esiste alcun rilascio programmato del SAS.
- Burst: WF4 simmetrico, finestra Hann, 4 cicli a 2.0597 Hz. Il SAS resta attivo e notchato durante e dopo il burst.
- Aerodinamica: 58 elementi C81 quasi-stazionari. Otto patch passive delle winglet riproducono area, rastremazione e quarto di corda dei CAERO1 NASTRAN 146001–149001 e 246001–249001.
"""
    (ROOT / "MODEL_AUDIT.md").write_text(audit)


def write_variable_inventory(paths: list[Path]) -> None:
    """Inventaria ogni variabile NetCDF, anche quelle escluse dal CSV ridotto."""
    lines = [
        "# Inventario completo delle variabili NetCDF MBDyn",
        "# Il CSV conserva le grandezze fisiche utili; il NetCDF raw resta integro.",
        "",
    ]
    for path in paths:
        lines.append(f"[{path.name}]")
        with Dataset(path) as data:
            for name in sorted(data.variables):
                variable = data.variables[name]
                units = getattr(variable, "units", "-")
                description = getattr(variable, "description", "-")
                shape = "x".join(str(value) for value in variable.shape) or "scalar"
                lines.append(
                    f"{name}\tdtype={variable.dtype}\tshape={shape}\tunits={units}\tdescription={description}"
                )
        lines.append("")
    (RESULTS / "netcdf_variables.txt").write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--velocities", nargs="+", type=float, help="analizza solo queste velocita'; default: tutti i V_*.nc")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    if args.velocities:
        paths = [paths_for(velocity)[1] for velocity in args.velocities]
    else:
        paths = sorted(RAW_OUTPUT.glob("V_*.nc"))
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("mancano: " + ", ".join(str(path) for path in missing))
    if not paths:
        raise FileNotFoundError(f"nessuna run V_*.nc in {RAW_OUTPUT}; eseguire prima run_sweep.py")
    write_variable_inventory(paths)
    # Il primo caso utilizzabile ad alta velocita' costituisce l'ancora del
    # ramo elastico; i casi successivi lo seguono verso velocita' decrescenti.
    summaries = []
    bff_target: float | None = None
    ordered_paths = sorted(
        paths,
        key=lambda value: constant(value.with_suffix(".mbd").read_text(), "V_INF"),
        reverse=True,
    )
    # NetCDF su /mnt/c e' molto lento quando vengono lette centinaia di
    # variabili piccole. Ogni caso viene copiato temporaneamente sul filesystem
    # Linux; CSV, figure e report continuano invece a finire in RESULTS.
    with tempfile.TemporaryDirectory(prefix="mbdyn_bff_analysis_", dir="/tmp") as scratch_name:
        scratch = Path(scratch_name)
        for index, path in enumerate(ordered_paths, 1):
            local_nc = scratch / path.name
            local_mbd = local_nc.with_suffix(".mbd")
            print(f"[analyse {index:02d}/{len(ordered_paths):02d}] {path.stem}: staging NetCDF...", flush=True)
            shutil.copyfile(path, local_nc)
            shutil.copyfile(path.with_suffix(".mbd"), local_mbd)
            summary = analyze_case(local_nc, bff_target)
            summaries.append(summary)
            if summary.identification_valid and np.isfinite(summary.bff_frequency_hz):
                bff_target = summary.bff_frequency_hz
            print(
                f"[analysed] V={summary.velocity_mps:.4f} m/s, "
                f"rigid={summary.short_period_frequency_hz:.3f}/{summary.short_period_sigma_per_s:+.3f}, "
                f"BFF={summary.bff_frequency_hz:.3f}/{summary.bff_sigma_per_s:+.3f}",
                flush=True,
            )
            local_nc.unlink()
            local_mbd.unlink()
    write_study(summaries)
    for value in sorted(summaries, key=lambda item: item.velocity_mps):
        print(f"V={value.velocity_mps:8.4f} m/s q={value.dynamic_pressure_pa:8.2f} Pa "
              f"rigid={value.short_period_frequency_hz:6.3f}/{value.short_period_sigma_per_s:+7.3f} "
              f"BFF={value.bff_frequency_hz:6.3f}/{value.bff_sigma_per_s:+7.3f} {value.classification}")
    print(f"Risultati: {RESULTS}")


if __name__ == "__main__":
    main()
