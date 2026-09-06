#!/usr/bin/env python3
"""Analyze balanced dive--pull-up runs and their paired BFF response."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bbf_dive_pullup_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.spatial.transform import Rotation

from analyse_time_domain_pairs import (
    case_metadata,
    constant,
    paired_body_flap_leakage,
    paired_tip_delta,
    trajectory_metrics,
)
from compare_paired_response import growth_metrics, paired_delta
from maneuver_case import CONFIG, baseline_module

IN_TO_M = 0.0254


def discover(cases: Path) -> dict[tuple[float, float], dict[bool, tuple[Path, Path, dict]]]:
    groups: dict[tuple[float, float], dict[bool, tuple[Path, Path, dict]]] = {}
    for nc in sorted(cases.glob("*.nc")):
        mbd = nc.with_suffix(".mbd")
        if not mbd.is_file():
            continue
        try:
            with Dataset(nc) as data:
                time = np.asarray(data["time"][:]).squeeze()
            if len(time) < 20 or time[-1] < constant(mbd, "FINAL_TIME") - 0.05:
                continue
            meta = case_metadata(mbd)
        except Exception:
            continue
        if meta["family"] != "dive_pullup":
            continue
        key = (float(meta["velocity_mps"]), float(meta["nominal_load_factor"]))
        groups.setdefault(key, {})[bool(meta["excited"])] = (nc, mbd, meta)
    return groups


def flight_metrics(nc: Path, mbd: Path) -> dict[str, float]:
    with Dataset(nc) as data:
        time = np.asarray(data["time"][:]).squeeze()
        position = np.asarray(data["node.struct.990000.X"][:]) * IN_TO_M
        velocity = np.asarray(data["node.struct.990000.XP"][:]) * IN_TO_M
        rotvec = np.asarray(data["node.struct.990000.Phi"][:])
        omega = np.degrees(np.asarray(data["node.struct.990000.Omega"][:])[:, 1])
    pitch = Rotation.from_rotvec(rotvec).as_euler("xyz", degrees=True)[:, 1]
    start = constant(mbd, "MANEUVER_START")
    end = constant(mbd, "DIVE_MANEUVER_END")
    final_time = constant(mbd, "FINAL_TIME")
    initial = (time >= start - 1.0) & (time < start)
    active = (time >= start) & (time <= end)
    final = (time >= final_time - 1.0) & (time <= final_time)
    h0 = float(np.mean(position[initial, 2]))
    theta0 = float(np.mean(pitch[initial]))
    return {
        "initial_pitch_deg": theta0,
        "minimum_pitch_increment_deg": float(np.min(pitch[active]) - theta0),
        "maximum_pitch_increment_deg": float(np.max(pitch[active]) - theta0),
        "minimum_altitude_change_m": float(np.min(position[active, 2]) - h0),
        "final_altitude_error_m": float(np.mean(position[final, 2]) - h0),
        "final_pitch_error_deg": float(np.mean(pitch[final]) - theta0),
        "final_q_mean_deg_s": float(np.mean(omega[final])),
        "final_vertical_speed_mean_m_s": float(np.mean(velocity[final, 2])),
    }


def onset_rows(rows: list[dict]) -> list[dict]:
    grouped: dict[float, list[dict]] = defaultdict(list)
    for row in rows:
        if row["identification_valid"]:
            grouped[row["nominal_load_class"]].append(row)
    result = []
    for load_factor in sorted({row["nominal_load_class"] for row in rows}):
        points = sorted(grouped[load_factor], key=lambda row: row["velocity_mps"])
        bracket = None
        estimate = None
        for left, right in zip(points, points[1:]):
            sl, sr = left["maneuver_sigma_per_s"], right["maneuver_sigma_per_s"]
            if sl == 0.0 or sl * sr <= 0.0:
                vl, vr = left["velocity_mps"], right["velocity_mps"]
                bracket = (vl, vr)
                estimate = vl if sl == 0.0 else vl - sl * (vr - vl) / (sr - sl)
                break
        result.append({
            "nominal_load_class": load_factor,
            "achieved_n_mean_across_speeds": (
                None if not points else float(np.mean([row["achieved_n_mean_sas_off"] for row in points]))
            ),
            "status": "crossing_found" if bracket else "no_bracket",
            "lower_velocity_mps": None if bracket is None else bracket[0],
            "upper_velocity_mps": None if bracket is None else bracket[1],
            "onset_velocity_linear_mps": estimate,
            "valid_speed_count": len(points),
        })
    return result


def analyze(campaign: Path) -> tuple[list[dict], list[dict], Path]:
    groups = discover(campaign / "cases")
    analysis = campaign / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for (velocity, nominal_load), pair in sorted(groups.items()):
        if set(pair) != {False, True}:
            continue
        tm, maneuver_mode = paired_delta(pair[False][0], pair[True][0])
        mbd = pair[True][1]
        frequency = baseline_module().dlm_rom_values(velocity)[3]
        release = constant(mbd, "SAS_OFF_START")
        duration = constant(mbd, "SAS_OFF_DURATION")
        fit_start = release + 0.05 + 0.742 / frequency + 0.05
        fit_end = release + duration - 0.05
        mg, _, _ = growth_metrics(tm, maneuver_mode, fit_start, fit_end, frequency)
        _, maneuver_tip = paired_tip_delta(pair[False][0], pair[True][0])
        mtg, _, _ = growth_metrics(tm, maneuver_tip, fit_start, fit_end, frequency)
        audit = trajectory_metrics(pair[False][0], pair[False][1])
        flight = flight_metrics(pair[False][0], pair[False][1])
        amplitude = float(pair[False][2].get("pitch_angle_deg", 0.0))
        q_command = float(pair[False][2].get("pitch_rate_deg_s", 0.0))
        leakage = paired_body_flap_leakage(pair[False][0], pair[True][0], fit_start, fit_end, frequency)
        trajectory_valid = bool(
            flight["minimum_pitch_increment_deg"] <= -0.35 * amplitude
            and flight["maximum_pitch_increment_deg"] >= 0.20 * amplitude
            and audit["achieved_n_mean_sas_off"] > 1.02
            and audit["achieved_n_std_sas_off"] <= 0.12
            and abs(audit["achieved_n_slope_sas_off_per_s"]) <= 0.15
            and abs(audit["q_mean_sas_off_deg_s"] - q_command) <= max(0.75, 0.30 * q_command)
            and abs(flight["final_pitch_error_deg"]) <= 1.0
            and abs(flight["final_q_mean_deg_s"]) <= 0.8
            and abs(flight["final_vertical_speed_mean_m_s"]) <= 1.0
            and abs(flight["final_altitude_error_m"]) <= 8.0
        )
        mode_tip_consistent = abs(mg["sigma_per_s"] - mtg["sigma_per_s"]) <= 0.25
        valid = bool(trajectory_valid and mode_tip_consistent and leakage <= 0.025)
        rows.append({
            "velocity_mps": velocity,
            "nominal_load_class": nominal_load,
            "pitch_rate_command_deg_s": q_command,
            "pitch_amplitude_deg": amplitude,
            **flight,
            **audit,
            "maneuver_sigma_per_s": mg["sigma_per_s"],
            "maneuver_tip_sigma_per_s": mtg["sigma_per_s"],
            "surface_bff_leakage_deg": leakage,
            "trajectory_valid": trajectory_valid,
            "mode7_tip_sigma_consistent": mode_tip_consistent,
            "identification_valid": valid,
            "verdict": "invalid" if not valid else (
                "unstable" if mg["sigma_per_s"] > 0.05 else
                "stable" if mg["sigma_per_s"] < -0.05 else "near_onset"
            ),
        })
    if not rows:
        raise SystemExit("nessuna coppia dive--pull-up completa")
    csv_path = analysis / "dive_pullup_comparison.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (analysis / "dive_pullup_comparison.json").write_text(json.dumps(rows, indent=2))
    onsets = onset_rows(rows)
    with (analysis / "dive_pullup_onset.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(onsets[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(onsets)
    return rows, onsets, analysis


def plot_summary(rows: list[dict], onsets: list[dict], analysis: Path) -> None:
    load_factors = sorted({row["nominal_load_class"] for row in rows})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for load_factor in load_factors:
        group = sorted((r for r in rows if r["nominal_load_class"] == load_factor), key=lambda r: r["velocity_mps"])
        velocity = [r["velocity_mps"] for r in group]
        label = f"classe nominale {load_factor:.2f}"
        axes[0, 0].plot(velocity, [r["maneuver_sigma_per_s"] for r in group], "o-", label=label)
        axes[0, 1].plot(velocity, [r["achieved_n_mean_sas_off"] for r in group], "o-", label=label)
        axes[1, 0].plot(velocity, [r["final_altitude_error_m"] for r in group], "o-", label=label)
    found = [row for row in onsets if row["onset_velocity_linear_mps"] is not None]
    if found:
        axes[1, 1].plot([r["achieved_n_mean_across_speeds"] for r in found], [r["onset_velocity_linear_mps"] for r in found], "o-")
    axes[0, 0].set(xlabel="TAS [m/s]", ylabel="sigma BFF [1/s]", title="Crescita BFF nella richiamata")
    axes[0, 1].set(xlabel="TAS [m/s]", ylabel="n medio misurato [-]", title="Carico realmente ottenuto dalla shadow")
    axes[1, 0].set(xlabel="TAS [m/s]", ylabel="errore quota finale [m]", title="Chiusura della traiettoria")
    axes[1, 1].set(xlabel="n medio misurato [-]", ylabel="onset [m/s]", title="Onset BBF durante la manovra")
    axes[0, 0].axvline(float(CONFIG["dive_pullup"]["reference_onset_mps"]), color="0.4", ls="--", label="onset longitudinale")
    for ax in axes.flat:
        ax.axhline(0.0, color="black", lw=0.8); ax.grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(analysis / "dive_pullup_envelope.png", dpi=190); plt.close(fig)


def plot_cases(campaign: Path, rows: list[dict], analysis: Path) -> None:
    groups = discover(campaign / "cases")
    for row in rows:
        key = (row["velocity_mps"], row["nominal_load_class"])
        pair = groups[key]
        with Dataset(pair[False][0]) as data:
            time = np.asarray(data["time"][:]).squeeze()
            position = np.asarray(data["node.struct.990000.X"][:]) * IN_TO_M
            pitch = Rotation.from_rotvec(np.asarray(data["node.struct.990000.Phi"][:])).as_euler("xyz", degrees=True)[:, 1]
            q = np.degrees(np.asarray(data["node.struct.990000.Omega"][:])[:, 1])
        _, delta = paired_delta(pair[False][0], pair[True][0])
        fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
        axes[0].plot(time, pitch); axes[0].set_ylabel("pitch [deg]")
        axes[1].plot(time, q); axes[1].set_ylabel("q [deg/s]")
        axes[2].plot(time, position[:, 2] - position[0, 2]); axes[2].set_ylabel("Delta h [m]")
        axes[3].plot(time, delta); axes[3].set(ylabel="Delta q7", xlabel="tempo [s]")
        release = constant(pair[False][1], "SAS_OFF_START")
        end = release + constant(pair[False][1], "SAS_OFF_DURATION")
        for ax in axes:
            ax.axvspan(release, end, color="tab:red", alpha=0.10); ax.grid(True, alpha=0.3)
        fig.suptitle(
            f"Dive--pull-up: V={key[0]:g} m/s, classe={key[1]:.2f}, "
            f"n misurato={row['achieved_n_mean_sas_off']:.3f}"
        )
        fig.tight_layout(); fig.savefig(analysis / f"case_V_{key[0]:06.2f}_class_{key[1]:04.2f}.png", dpi=180); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--no-case-plots", action="store_true")
    args = parser.parse_args()
    campaign = args.campaign.expanduser().resolve()
    rows, onsets, analysis = analyze(campaign)
    plot_summary(rows, onsets, analysis)
    if not args.no_case_plots:
        plot_cases(campaign, rows, analysis)
    (analysis / "dive_pullup_summary.json").write_text(json.dumps({
        "valid_points": sum(row["identification_valid"] for row in rows),
        "total_points": len(rows),
        "onsets": onsets,
        "reference_longitudinal_onset_mps": CONFIG["dive_pullup"]["reference_onset_mps"],
        "load_note": "n is measured from each shadow run; nominal classes are not feedback targets",
        "interpretation": CONFIG["rom"]["warning"],
    }, indent=2))
    print(analysis / "dive_pullup_envelope.png")


if __name__ == "__main__":
    main()
