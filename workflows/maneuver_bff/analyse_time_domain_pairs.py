#!/usr/bin/env python3
"""Compare steady and maneuver BFF responses from paired MBDyn trajectories."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.signal import butter, sosfiltfilt
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT / "study_config.json").read_text())
BASELINE_VALUE = Path(CONFIG["baseline_directory"])
BASELINE = BASELINE_VALUE if BASELINE_VALUE.is_absolute() else (ROOT / BASELINE_VALUE).resolve()
sys.path.insert(0, str(BASELINE))

from run_case import dlm_rom_values  # noqa: E402
from compare_paired_response import growth_metrics, paired_delta  # noqa: E402

G = 9.81
IN_TO_M = 0.0254


def detrended_harmonic_amplitude(
    time: np.ndarray, signal: np.ndarray, frequency: float,
) -> float:
    shifted_time = np.asarray(time, float) - float(time[0])
    values = np.asarray(signal, float)
    design = np.column_stack((
        np.ones_like(shifted_time), shifted_time,
        np.sin(2.0 * math.pi * frequency * shifted_time),
        np.cos(2.0 * math.pi * frequency * shifted_time),
    ))
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return float(math.hypot(coefficients[2], coefficients[3]))


def read_symmetric_tip(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return bending-like symmetric tip displacement in the body frame."""
    with Dataset(path) as data:
        time = np.asarray(data["time"][:]).squeeze()
        base = np.asarray(data["node.struct.990000.X"][:])
        phi = np.asarray(data["node.struct.990000.Phi"][:])
        left = np.asarray(data["node.struct.990020.X"][:])
        right = np.asarray(data["node.struct.991020.X"][:])
    rotation = Rotation.from_rotvec(phi)
    left_body = rotation.inv().apply(left - base)[:, 2] * IN_TO_M
    right_body = rotation.inv().apply(right - base)[:, 2] * IN_TO_M
    return time, 0.5 * (left_body + right_body)


def paired_tip_delta(shadow: Path, excited: Path) -> tuple[np.ndarray, np.ndarray]:
    ts, ys = read_symmetric_tip(shadow)
    te, ye = read_symmetric_tip(excited)
    if len(ts) != len(te) or not np.allclose(ts, te, atol=1e-9, rtol=0.0):
        raise ValueError("assi temporali tip shadow/excited differenti")
    return ts, ye - ys


def read_surface_angle(path: Path, parent: int, control: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the relative surface rotation magnitude in degrees."""
    with Dataset(path) as data:
        time = np.asarray(data["time"][:]).squeeze()
        parent_rotation = Rotation.from_rotvec(np.asarray(data[f"node.struct.{parent}.Phi"][:]))
        control_rotation = Rotation.from_rotvec(np.asarray(data[f"node.struct.{control}.Phi"][:]))
    relative = (parent_rotation.inv() * control_rotation).as_rotvec()
    return time, np.degrees(np.linalg.norm(relative, axis=1))


def paired_body_flap_leakage(
    shadow: Path, excited: Path, start: float, end: float, frequency: float,
) -> float:
    amplitudes = []
    for parent, control in ((990004, 880004), (991004, 881004)):
        ts, ys = read_surface_angle(shadow, parent, control)
        te, ye = read_surface_angle(excited, parent, control)
        if len(ts) != len(te) or not np.allclose(ts, te, atol=1e-9, rtol=0.0):
            raise ValueError("assi temporali body-flap shadow/excited differenti")
        use = (ts >= start) & (ts <= end)
        amplitudes.append(detrended_harmonic_amplitude(ts[use], (ye - ys)[use], frequency))
    return float(max(amplitudes))


def case_metadata(mbd: Path) -> dict:
    match = re.search(r"(?m)^# MANEUVER_METADATA (\{.*\})$", mbd.read_text())
    if not match:
        raise ValueError(f"metadata mancanti in {mbd}")
    return json.loads(match.group(1))


def constant(mbd: Path, name: str) -> float:
    match = re.search(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*([-+0-9.eE]+)\s*;",
        mbd.read_text(),
    )
    if not match:
        raise ValueError(f"{name} non numerico o assente in {mbd}")
    return float(match.group(1))


def trajectory_metrics(nc: Path, mbd: Path) -> dict:
    with Dataset(nc) as data:
        time = np.asarray(data["time"][:]).squeeze()
        position = np.asarray(data["node.struct.990000.X"][:]) * IN_TO_M
        velocity = np.asarray(data["node.struct.990000.XP"][:]) * IN_TO_M
        rotation_vector = np.asarray(data["node.struct.990000.Phi"][:])
        omega = np.asarray(data["node.struct.990000.Omega"][:])
    release = constant(mbd, "SAS_OFF_START")
    duration = constant(mbd, "SAS_OFF_DURATION")
    if len(time) < 20 or time[-1] < release + duration:
        raise ValueError(f"run incompleta: {nc}")
    dt = float(np.median(np.diff(time)))
    acceleration = np.gradient(velocity, dt, axis=0)
    specific_global = acceleration + np.array([0.0, 0.0, G])
    body_specific = Rotation.from_rotvec(rotation_vector).inv().apply(specific_global)
    nz_raw = body_specific[:, 2] / G
    # Validate the slow maneuver trajectory without counting the ~2 Hz BFF
    # oscillation itself as load-factor tracking error.
    sos = butter(3, 0.75, btype="lowpass", fs=1.0 / dt, output="sos")
    nz = sosfiltfilt(sos, nz_raw)
    euler_deg = Rotation.from_rotvec(rotation_vector).as_euler("xyz", degrees=True)
    pitch_deg = sosfiltfilt(sos, euler_deg[:, 1])
    pre = (time >= release - 0.5) & (time < release)
    sas_off = (time >= release) & (time <= release + duration)
    off_time = time[sas_off]
    off_nz = nz[sas_off]
    off_slope = np.polyfit(off_time - off_time[0], off_nz, 1)[0]
    off_pitch = pitch_deg[sas_off]
    pitch_slope = np.polyfit(off_time - off_time[0], off_pitch, 1)[0]
    at_release = min(int(np.searchsorted(time, release)), len(time) - 1)
    return {
        "achieved_n_mean_pre_release": float(np.mean(nz[pre])),
        "achieved_n_std_pre_release": float(np.std(nz[pre])),
        "achieved_n_mean_sas_off": float(np.mean(off_nz)),
        "achieved_n_std_sas_off": float(np.std(off_nz)),
        "achieved_n_min_sas_off": float(np.min(off_nz)),
        "achieved_n_max_sas_off": float(np.max(off_nz)),
        "achieved_n_slope_sas_off_per_s": float(off_slope),
        "q_mean_sas_off_deg_s": float(math.degrees(np.mean(omega[sas_off, 1]))),
        "q_std_sas_off_deg_s": float(math.degrees(np.std(omega[sas_off, 1]))),
        "q_release_deg_s": float(math.degrees(omega[at_release, 1])),
        "p_release_deg_s": float(math.degrees(omega[at_release, 0])),
        "pitch_mean_pre_release_deg": float(np.mean(pitch_deg[pre])),
        "pitch_mean_sas_off_deg": float(np.mean(off_pitch)),
        "pitch_std_sas_off_deg": float(np.std(off_pitch)),
        "pitch_slope_sas_off_deg_s": float(pitch_slope),
        "altitude_change_sas_off_m": float(position[sas_off][-1, 2] - position[sas_off][0, 2]),
        "vertical_speed_mean_sas_off_m_s": float(np.mean(velocity[sas_off, 2])),
    }


def discover(
    cases: Path, allowed_stems: set[str] | None = None,
) -> dict[tuple, dict[bool, tuple[Path, Path, dict]]]:
    grouped: dict[tuple, dict[bool, tuple[Path, Path, dict]]] = {}
    for nc in sorted(cases.glob("*.nc")):
        if allowed_stems is not None and nc.stem not in allowed_stems:
            continue
        mbd = nc.with_suffix(".mbd")
        if not mbd.is_file():
            continue
        try:
            with Dataset(nc) as data:
                time = np.asarray(data["time"][:]).squeeze()
            required_end = constant(mbd, "SAS_OFF_START") + constant(mbd, "SAS_OFF_DURATION")
            if time.ndim != 1 or len(time) < 20 or float(time[-1]) < required_end:
                print(f"[skip] run incompleta: {nc.name}")
                continue
        except Exception as error:
            print(f"[skip] NetCDF non leggibile {nc.name}: {error}")
            continue
        meta = case_metadata(mbd)
        key = (
            meta["family"], float(meta["velocity_mps"]), float(meta["load_factor"]),
            float(meta["bank_angle_deg"]),
        )
        grouped.setdefault(key, {})[bool(meta["excited"])] = (nc, mbd, meta)
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--manifest", type=Path, action="append",
        help="analizza soltanto gli stem presenti nei manifest (ripetibile)",
    )
    args = parser.parse_args()
    output = args.output or args.cases.parent / "time_domain_analysis"
    output.mkdir(parents=True, exist_ok=True)
    allowed_stems = None
    if args.manifest:
        allowed_stems = set()
        for manifest in args.manifest:
            with manifest.open(newline="") as stream:
                allowed_stems.update(row["stem"] for row in csv.DictReader(stream))
    grouped = discover(args.cases, allowed_stems)
    rows: list[dict] = []
    for key, pair in sorted(grouped.items()):
        family, velocity, load_factor, bank = key
        if set(pair) != {False, True}:
            print(f"[skip] coppia incompleta: {key}")
            continue
        steady_key = ("pullup", velocity, 1.0, 0.0)
        steady = grouped.get(steady_key)
        if steady is None or set(steady) != {False, True}:
            print(f"[skip] riferimento steady mancante: {key}")
            continue
        ts, ds = paired_delta(steady[False][0], steady[True][0])
        tm, dm = paired_delta(pair[False][0], pair[True][0])
        if len(ts) != len(tm) or not np.allclose(ts, tm, atol=1e-9, rtol=0.0):
            print(f"[skip] griglie temporali incompatibili: {key}")
            continue
        mbd = pair[True][1]
        frequency = dlm_rom_values(velocity)[3]
        release = constant(mbd, "SAS_OFF_START")
        duration = constant(mbd, "SAS_OFF_DURATION")
        rap_duration = 0.742 / frequency
        start = release + 0.05 + rap_duration + 0.05
        end = release + duration - 0.05
        steady_growth, _, _ = growth_metrics(ts, ds, start, end, frequency)
        maneuver_growth, _, _ = growth_metrics(tm, dm, start, end, frequency)
        delta_sigma = maneuver_growth["sigma_per_s"] - steady_growth["sigma_per_s"]
        _, steady_tip = paired_tip_delta(steady[False][0], steady[True][0])
        _, maneuver_tip = paired_tip_delta(pair[False][0], pair[True][0])
        steady_tip_growth, _, _ = growth_metrics(ts, steady_tip, start, end, frequency)
        maneuver_tip_growth, _, _ = growth_metrics(tm, maneuver_tip, start, end, frequency)
        tip_delta_sigma = maneuver_tip_growth["sigma_per_s"] - steady_tip_growth["sigma_per_s"]
        audit = trajectory_metrics(pair[False][0], pair[False][1])
        keeper_leakage = paired_body_flap_leakage(
            pair[False][0], pair[True][0], start, end, frequency,
        )
        trajectory_valid = bool(
            abs(audit["achieved_n_mean_sas_off"] - load_factor) <= 0.075
            and audit["achieved_n_std_sas_off"] <= 0.10
            and abs(audit["achieved_n_slope_sas_off_per_s"]) <= 0.10
        )
        mode_tip_consistent = abs(
            maneuver_growth["sigma_per_s"] - maneuver_tip_growth["sigma_per_s"]
        ) <= 0.25
        identification_valid = bool(
            trajectory_valid and keeper_leakage <= 0.025 and mode_tip_consistent
        )
        physical_verdict = (
            "amplified" if delta_sigma > 0.05 else
            "suppressed" if delta_sigma < -0.05 else
            "indistinguishable"
        )
        rows.append({
            "family": family,
            "velocity_mps": velocity,
            "load_factor_command": load_factor,
            "bank_angle_deg": bank,
            **audit,
            "load_factor_error": audit["achieved_n_mean_sas_off"] - load_factor,
            "steady_sigma_per_s": steady_growth["sigma_per_s"],
            "maneuver_sigma_per_s": maneuver_growth["sigma_per_s"],
            "delta_sigma_per_s": delta_sigma,
            "steady_tip_sigma_per_s": steady_tip_growth["sigma_per_s"],
            "maneuver_tip_sigma_per_s": maneuver_tip_growth["sigma_per_s"],
            "tip_delta_sigma_per_s": tip_delta_sigma,
            "maneuver_keeper_bff_leakage_deg": keeper_leakage,
            "trajectory_valid": trajectory_valid,
            "mode7_tip_sigma_consistent": mode_tip_consistent,
            "identification_valid": identification_valid,
            "verdict": physical_verdict if identification_valid else "invalid",
        })
    if not rows:
        raise SystemExit("nessuna coppia completa analizzabile")
    csv_path = output / "time_domain_maneuver_comparison.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "time_domain_maneuver_comparison.json").write_text(json.dumps(rows, indent=2))
    print(csv_path)


if __name__ == "__main__":
    main()
