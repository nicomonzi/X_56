#!/usr/bin/env python3
"""Plan or run the focused dive--pull-up BFF campaign (no run by default)."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

from maneuver_case import CONFIG, ManeuverPoint, model_revision, run_point

ROOT = Path(__file__).resolve().parent


def write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "case_index", "family", "velocity_mps", "load_factor",
        "bank_angle_deg", "excited", "pitch_angle_deg", "pitch_rate_deg_s",
        "nominal_load_factor", "stem", "status", "message",
    )
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    path.with_suffix(".json").write_text(json.dumps(rows, indent=2))


def execute_point(point: ManeuverPoint, cases: Path, overwrite: bool) -> tuple[str, str]:
    target = cases / f"{point.stem}.nc"
    input_file = target.with_suffix(".mbd")
    input_text = input_file.read_text() if input_file.is_file() else ""
    reusable = bool(
        target.is_file()
        and input_file.is_file()
        and f"# MANEUVER_MODEL_REVISION {model_revision(point)}" in input_text
    )
    if reusable and not overwrite:
        return "complete", "risultato definitivo esistente riutilizzato"
    try:
        run_point(point, cases, overwrite=overwrite or target.exists() or input_file.exists())
    except Exception as error:
        return "failed", str(error)
    return "complete", "run MBDyn completata"


def command_for_class(velocity: float, nominal_load: float) -> tuple[float, float]:
    """Return q command and dive amplitude; actual n is never fed back."""
    q_deg_s = (nominal_load - 1.0) * 9.81 / velocity * 180.0 / math.pi
    timing = CONFIG["dive_pullup_timing"]
    centered_duration = (
        float(timing["pullup_release_lead_s"])
        + float(timing["sas_off_duration_s"])
    )
    amplitude_deg = 0.5 * q_deg_s * centered_duration
    return q_deg_s, amplitude_deg


def build_points(velocities: list[float], nominal_loads: list[float], shadow_only: bool = False) -> list[ManeuverPoint]:
    points: list[ManeuverPoint] = []
    for velocity in velocities:
        for nominal_load in nominal_loads:
            q_command, amplitude = command_for_class(velocity, nominal_load)
            points.extend(
                ManeuverPoint(
                    "dive_pullup", velocity, 1.0, 0.0, excited, amplitude,
                    q_command, nominal_load,
                )
                for excited in ((False,) if shadow_only else (False, True))
            )
    return points


def parse_args() -> argparse.Namespace:
    cfg = CONFIG["dive_pullup"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(cfg["default_output"]))
    parser.add_argument(
        "--velocities", type=float, nargs="+", default=cfg["velocities_mps"],
        help="velocita' mirate attorno all'onset longitudinale",
    )
    parser.add_argument(
        "--nominal-loads", type=float, nargs="+", default=cfg["nominal_load_classes"],
        help="classi nominali usate solo per ricavare q; n viene misurato dalla shadow",
    )
    parser.add_argument("--shadow-only", action="store_true",
                        help="prima fase: esegue soltanto le shadow, riutilizzabili dopo")
    parser.add_argument("--execute", action="store_true", help="senza questa opzione crea solo il manifest")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--no-analysis", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.jobs < 1 or not args.velocities or not args.nominal_loads:
        raise SystemExit("lista dei casi o numero di job non valido")
    velocities = sorted(set(float(value) for value in args.velocities))
    nominal_loads = sorted(set(float(value) for value in args.nominal_loads))
    if any(not 50.0 <= value <= 70.0 for value in velocities):
        raise SystemExit("le velocita' devono essere comprese tra 50 e 70 m/s")
    if any(not 1.0 < value <= 1.8 for value in nominal_loads):
        raise SystemExit("le classi nominali devono soddisfare 1.0 < n_nominale <= 1.8")
    planned = build_points(velocities, nominal_loads, args.shadow_only)
    output = args.output.expanduser().resolve()
    manifest = output / "dive_pullup_manifest.csv"
    rows = []
    for index, point in enumerate(planned, start=1):
        row = asdict(point)
        row.update(case_index=index, stem=point.stem, status="planned", message="")
        rows.append(row)
    write_manifest(rows, manifest)
    print(
        f"[plan] {len(velocities)} velocita' x {len(nominal_loads)} aggressivita' "
        f"x {'solo shadow' if args.shadow_only else 'shadow/excited'} = {len(planned)} traiettorie"
    )
    for velocity in velocities:
        mapping = ", ".join(
            (
                f"classe={load:.2f}: qcmd={command_for_class(velocity, load)[0]:.3f} deg/s, "
                f"A={command_for_class(velocity, load)[1]:.3f} deg"
            )
            for load in nominal_loads
        )
        print(f"[mapping] V={velocity:.1f} m/s: {mapping}")
    print(f"[output] {output}")
    if not args.execute:
        print("[plan-only] nessuna run avviata; aggiungere --execute")
        return

    cases = output / "cases"
    cases.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        pending = {
            executor.submit(execute_point, point, cases, args.overwrite): index
            for index, point in enumerate(planned)
        }
        for completed, future in enumerate(as_completed(pending), start=1):
            index = pending[future]
            status, message = future.result()
            rows[index]["status"] = status
            rows[index]["message"] = message
            write_manifest(rows, manifest)
            print(f"[{completed}/{len(planned)}] {rows[index]['stem']}: {status}", flush=True)
    failures = sum(row["status"] == "failed" for row in rows)
    print(f"[done] complete={len(rows) - failures}, failed={failures}")
    if not args.no_analysis and not args.shadow_only:
        subprocess.run(
            [sys.executable, str(ROOT / "analyse_dive_pullup.py"), str(output)],
            check=True,
        )


if __name__ == "__main__":
    main()
