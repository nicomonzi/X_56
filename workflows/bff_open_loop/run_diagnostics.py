#!/usr/bin/env python3
"""Run and compare only the eight prescribed open-loop diagnostic cases."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bff_open_loop_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyse_open_loop import analyze_case
from run_case import DEFAULT_OUTPUT, run_case

# Le varianti immediate/delayed appartengono al protocollo precedente. I dati
# storici restano leggibili, ma nuove diagnostiche devono usare run_sweep.py.
FREEZE_VARIANTS: tuple[str, ...] = ()

VELOCITIES = (40.0, 45.0, 50.0, 55.0)


def compatible(a: dict, b: dict) -> tuple[bool, float, float]:
    df = abs(float(a["frequency_hz"]) - float(b["frequency_hz"]))
    ds = abs(float(a["sigma_per_s"]) - float(b["sigma_per_s"]))
    fmean = 0.5 * (float(a["frequency_hz"]) + float(b["frequency_hz"]))
    smean = 0.5 * (abs(float(a["sigma_per_s"])) + abs(float(b["sigma_per_s"])))
    return df <= max(0.30, 0.12 * fmean) and ds <= max(1.8, 0.50 * smean), df, ds


def cross_freeze_points(summaries: list[dict]) -> list[dict]:
    indexed = {(float(item["velocity_mps"]), item["freeze_variant"]): item for item in summaries}
    points = []
    for velocity in VELOCITIES:
        immediate = indexed[(velocity, "immediate")]
        delayed = indexed[(velocity, "delayed")]
        ca, cb = immediate["matrix_pencil_candidate"], delayed["matrix_pencil_candidate"]
        accepted = bool(immediate["identification_valid"] and delayed["identification_valid"])
        agree, df, ds = compatible(ca, cb) if accepted else (False, math.nan, math.nan)
        qualities = np.array([float(ca.get("quality_score", 0.0)), float(cb.get("quality_score", 0.0))])
        weights = np.maximum(qualities, 1e-6)
        point = {
            "velocity_mps": velocity,
            "cross_freeze_accepted": bool(accepted and agree),
            "frequency_hz": float(np.average([ca.get("frequency_hz", math.nan), cb.get("frequency_hz", math.nan)], weights=weights)) if accepted else None,
            "sigma_per_s": float(np.average([ca.get("sigma_per_s", math.nan), cb.get("sigma_per_s", math.nan)], weights=weights)) if accepted else None,
            "frequency_difference_hz": None if not accepted else df,
            "sigma_difference_per_s": None if not accepted else ds,
            "mean_local_quality": float(np.mean(qualities)),
            "immediate_case": immediate["case"], "delayed_case": delayed["case"],
        }
        if point["frequency_hz"] is not None:
            omega_n = math.hypot(point["sigma_per_s"], 2.0 * math.pi * point["frequency_hz"])
            point["damping_ratio"] = -point["sigma_per_s"] / omega_n if omega_n else None
        else:
            point["damping_ratio"] = None
        points.append(point)
    return points


def mark_continuous_track(points: list[dict]) -> None:
    for index, point in enumerate(points):
        neighbors = []
        for other_index in (index - 1, index + 1):
            if 0 <= other_index < len(points) and point["cross_freeze_accepted"] and points[other_index]["cross_freeze_accepted"]:
                agree, df, ds = compatible(point, points[other_index])
                neighbors.append({"velocity_mps": points[other_index]["velocity_mps"], "compatible": agree, "df_hz": df, "ds_per_s": ds})
        tracked = bool(any(item["compatible"] for item in neighbors))
        point["adjacent_comparisons"] = neighbors
        point["bff_tracked"] = tracked
        point["classification"] = (
            "no_robust_bff_candidate" if not tracked else
            "tracked_unstable_bff_candidate" if point["sigma_per_s"] > 0.0 else
            "tracked_stable_bff_candidate" if point["sigma_per_s"] < 0.0 else
            "tracked_neutral_bff_candidate"
        )


def write_outputs(summaries: list[dict], points: list[dict], output: Path) -> None:
    rows = []
    point_by_speed = {point["velocity_mps"]: point for point in points}
    for summary in sorted(summaries, key=lambda item: (item["velocity_mps"], item["freeze_variant"])):
        candidate = summary["matrix_pencil_candidate"]
        tracked = point_by_speed[summary["velocity_mps"]]
        rows.append({
            "velocity_mps": summary["velocity_mps"], "freeze_variant": summary["freeze_variant"],
            "open_loop_valid": summary["open_loop_valid"], "trim_stationary": summary["trim_stationary"],
            "local_candidate_accepted": candidate.get("accepted", False),
            "local_confidence": candidate.get("confidence", "none"),
            "local_quality": candidate.get("quality_score"),
            "local_frequency_hz": candidate.get("frequency_hz"), "local_sigma_per_s": candidate.get("sigma_per_s"),
            "selected_window_s": summary.get("selected_window_s"),
            "cross_freeze_accepted": tracked["cross_freeze_accepted"], "bff_tracked": tracked["bff_tracked"],
            "tracked_frequency_hz": tracked["frequency_hz"] if tracked["bff_tracked"] else None,
            "tracked_sigma_per_s": tracked["sigma_per_s"] if tracked["bff_tracked"] else None,
            "classification": tracked["classification"],
        })
    with (output / "diagnostic_cases.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    trim_rows = []
    for summary in sorted(summaries, key=lambda item: (item["velocity_mps"], item["freeze_variant"])):
        snapshots = json.loads((output / summary["trim_check_file"]).read_text())
        for snapshot_name, snapshot in snapshots.items():
            row = {
                "velocity_mps": summary["velocity_mps"], "freeze_variant": summary["freeze_variant"],
                "snapshot": snapshot_name, "time_s": snapshot["time_s"], "theta_deg": snapshot["theta_deg"],
                "q_deg_s": snapshot["q_deg_s"], "Vz_mps": snapshot["Vz_mps"], "Z_m": snapshot["Z_m"],
                "Fz_residual_N": snapshot["Fz_residual_N"], "My_residual_Nm": snapshot["My_residual_Nm"],
            }
            for group in ("actual_surfaces_deg", "surface_commands_deg", "pid_outputs_deg", "actuator_states_deg"):
                for name, value in snapshot[group].items():
                    row[f"{group}.{name}"] = value
            trim_rows.append(row)
    trim_fields = list(dict.fromkeys(key for row in trim_rows for key in row))
    with (output / "trim_snapshots.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=trim_fields); writer.writeheader(); writer.writerows(trim_rows)
    report = {
        "scope": "diagnostic only: 40/45/50/55 m/s, immediate and delayed freeze; not a flutter-boundary sweep",
        "selection_rule": "common persistent q/SWB1/SWB1dot/symmetric-tip pole, cross-freeze agreement, and adjacent-speed continuity",
        "trim_table": "trim_snapshots.csv (complete actual surfaces, commanded surfaces, PID outputs, and actuator states)",
        "points": points,
        "cases": rows,
    }
    (output / "diagnostic_summary.json").write_text(json.dumps(report, indent=2))

    lines = [
        "# Open-loop longitudinal diagnostic comparison", "",
        "These are diagnostic cases only; no 50--70 m/s boundary sweep was run.", "",
        "| V [m/s] | Freeze | OL valid | Local candidate | Confidence | f [Hz] | sigma [1/s] | Tracked BFF |",
        "|---:|:---|:---:|:---:|:---:|---:|---:|:---:|",
    ]
    for row in rows:
        f = "—" if row["local_frequency_hz"] is None else f"{row['local_frequency_hz']:.4f}"
        s = "—" if row["local_sigma_per_s"] is None else f"{row['local_sigma_per_s']:.4f}"
        lines.append(f"| {row['velocity_mps']:.0f} | {row['freeze_variant']} | {row['open_loop_valid']} | {row['local_candidate_accepted']} | {row['local_confidence']} | {f} | {s} | {row['bff_tracked']} |")
    lines.extend(["", "## Cross-freeze consensus", "", "| V [m/s] | Cross-freeze | f [Hz] | sigma [1/s] | Classification |", "|---:|:---:|---:|---:|:---|"])
    for point in points:
        f = "—" if point["frequency_hz"] is None else f"{point['frequency_hz']:.4f}"
        s = "—" if point["sigma_per_s"] is None else f"{point['sigma_per_s']:.4f}"
        lines.append(f"| {point['velocity_mps']:.0f} | {point['cross_freeze_accepted']} | {f} | {s} | {point['classification']} |")
    lines.extend([
        "", "## Trim snapshots", "",
        "The compact table below shows flight state and force/moment residuals. `trim_snapshots.csv` and each case `*_trim_check.json` also contain all actual and commanded surfaces, PID outputs, and actuator states.", "",
        "| V | Freeze | Snapshot | theta [deg] | q [deg/s] | Vz [m/s] | Z [m] | Fz res. [N] | My res. [Nm] |",
        "|---:|:---|:---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in trim_rows:
        lines.append(f"| {row['velocity_mps']:.0f} | {row['freeze_variant']} | {row['snapshot']} | {row['theta_deg']:.4f} | {row['q_deg_s']:.4f} | {row['Vz_mps']:.4f} | {row['Z_m']:.4f} | {row['Fz_residual_N']:.2f} | {row['My_residual_Nm']:.2f} |")
    (output / "DIAGNOSTIC_REPORT.md").write_text("\n".join(lines) + "\n")

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    for variant, marker in (("immediate", "o"), ("delayed", "s")):
        selected = sorted([item for item in summaries if item["freeze_variant"] == variant and item["identification_valid"]], key=lambda item: item["velocity_mps"])
        if selected:
            axes[0].plot([x["velocity_mps"] for x in selected], [x["frequency_swb1_hz"] for x in selected], marker + "--", alpha=.55, label=f"{variant} local")
            axes[1].plot([x["velocity_mps"] for x in selected], [x["sigma_swb1_per_s"] for x in selected], marker + "--", alpha=.55, label=f"{variant} local")
            axes[2].plot([x["velocity_mps"] for x in selected], [x["damping_ratio_swb1"] for x in selected], marker + "--", alpha=.55, label=f"{variant} local")
    tracked = [point for point in points if point["bff_tracked"]]
    if tracked:
        axes[0].plot([x["velocity_mps"] for x in tracked], [x["frequency_hz"] for x in tracked], "k*-", ms=10, label="cross-freeze tracked")
        axes[1].plot([x["velocity_mps"] for x in tracked], [x["sigma_per_s"] for x in tracked], "k*-", ms=10, label="cross-freeze tracked")
        axes[2].plot([x["velocity_mps"] for x in tracked], [x["damping_ratio"] for x in tracked], "k*-", ms=10, label="cross-freeze tracked")
    axes[0].set_ylabel("frequency [Hz]"); axes[1].set_ylabel("sigma [1/s]"); axes[2].set_ylabel("damping ratio"); axes[2].set_xlabel("V_INF [m/s]")
    axes[1].axhline(0.0, color="k", lw=1); axes[2].axhline(0.0, color="k", lw=1)
    for axis in axes: axis.grid(alpha=.2); axis.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(output / "diagnostic_tracking.png", dpi=180); plt.close(fig)


def main() -> None:
    raise SystemExit(
        "run_diagnostics.py e' un driver legacy immediate/delayed; "
        "usare run_case.py o run_sweep.py per il protocollo bounded SAS-off"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--jobs", type=int, default=2, choices=(1, 2, 3, 4))
    parser.add_argument("--analyse-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = [(velocity, variant) for velocity in VELOCITIES for variant in FREEZE_VARIANTS]
    if not args.analyse_only:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run_case, velocity, args.output, args.overwrite, variant): (velocity, variant) for velocity, variant in cases}
            for future in as_completed(futures):
                future.result()
    summaries = []
    for velocity, variant in cases:
        nc_path = args.output / (f"OL_{variant}_V_{velocity:08.4f}".replace(".", "p") + ".nc")
        if not nc_path.exists():
            raise FileNotFoundError(nc_path)
        summaries.append(analyze_case(nc_path, args.output))
    points = cross_freeze_points(summaries)
    mark_continuous_track(points)
    write_outputs(summaries, points, args.output)
    for summary in sorted(summaries, key=lambda item: (item["velocity_mps"], item["freeze_variant"])):
        snapshots = json.loads((args.output / summary["trim_check_file"]).read_text())
        for name, value in snapshots.items():
            print(
                f"[trim] V={summary['velocity_mps']:.0f} {summary['freeze_variant']} {name}: "
                f"theta={value['theta_deg']:.4f} deg, q={value['q_deg_s']:.4f} deg/s, "
                f"Vz={value['Vz_mps']:.4f} m/s, Z={value['Z_m']:.4f} m, "
                f"Fz_res={value['Fz_residual_N']:.2f} N, My_res={value['My_residual_Nm']:.2f} Nm; "
                f"surfaces/PID/actuators={summary['trim_check_file']}"
            )
    print(f"Diagnostic report: {args.output / 'DIAGNOSTIC_REPORT.md'}")


if __name__ == "__main__":
    main()
