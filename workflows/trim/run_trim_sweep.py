#!/usr/bin/env python3
"""Completa lo sweep MBDyn del trim accoppiato 2x2 dell'X-56.

Default: 30--70 m/s inclusi, incremento 2.5 m/s, densita' NASTRAN.
I casi gia' terminati vengono conservati e saltati; ogni nuovo caso salva input,
stdout/stderr e output MBDyn in una cartella propria. Al termine vengono creati
automaticamente CSV e grafici riepilogativi.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "main_trim.mbd"
DEFAULT_RESULTS = Path("/mnt/c/Users/Utente/Desktop/TRIM")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--vmin", type=float, default=30.0)
    parser.add_argument("--vmax", type=float, default=70.0)
    parser.add_argument("--dv", type=float, default=2.5)
    parser.add_argument(
        "--rho",
        type=float,
        default=9.7284e-8,
        help="densita' coerente IPS [lbf s^2/in^4] (default: NASTRAN)",
    )
    parser.add_argument("--mbdyn", type=Path, default=None)
    resume = parser.add_mutually_exclusive_group()
    resume.add_argument(
        "--resume",
        dest="resume",
        action="store_true",
        help="salta i casi completi (comportamento predefinito)",
    )
    resume.add_argument(
        "--force",
        dest="resume",
        action="store_false",
        help="riesegue anche i casi gia' completi, sovrascrivendone i risultati",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="continua lo sweep se un singolo caso fallisce",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="genera e controlla gli input senza lanciare MBDyn",
    )
    parser.add_argument(
        "--analyse",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="esegue analyze_trim_sweep.py al termine dello sweep",
    )
    parser.set_defaults(resume=True)
    return parser.parse_args()


def locate_mbdyn(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.expanduser())
    if os.environ.get("MBDYN_BIN"):
        candidates.append(Path(os.environ["MBDYN_BIN"]).expanduser())
    found = shutil.which("mbdyn")
    if found:
        candidates.append(Path(found))
    candidates.extend(
        [Path("/usr/local/mbdyn/bin/mbdyn"), Path("/usr/local/bin/mbdyn")]
    )
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise FileNotFoundError("MBDyn non trovato; usare --mbdyn /percorso/mbdyn")


def velocity_grid(vmin: float, vmax: float, step: float) -> list[float]:
    if vmin <= 0.0 or vmax < vmin or step <= 0.0:
        raise ValueError("Richiesti vmin>0, vmax>=vmin e dv>0")
    count = int(math.floor((vmax - vmin) / step + 1.0e-10)) + 1
    values = [vmin + index * step for index in range(count)]
    if values[-1] < vmax - 1.0e-9:
        values.append(vmax)
    else:
        values[-1] = vmax
    return values


def replace_real_constant(source: str, name: str, value: float) -> str:
    pattern = re.compile(
        rf"(set\s*:\s*(?:ifndef\s+)?const\s+real\s+{re.escape(name)}\s*=\s*)([^;]+)(;)",
        flags=re.IGNORECASE,
    )
    replaced, count = pattern.subn(
        lambda match: f"{match.group(1)}{value:.16g}{match.group(3)}",
        source,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"Dichiarazione unica di {name} non trovata nell'input")
    return replaced


def generated_source(template: str, velocity: float, density: float) -> str:
    source = replace_real_constant(template, "V_INF", velocity)
    return replace_real_constant(source, "RHO_AIR", density)


def completed_case(case_directory: Path) -> bool:
    output = case_directory / "result.out"
    netcdf = case_directory / "result.nc"
    if not output.is_file() or not netcdf.is_file():
        return False
    tail = output.read_text(encoding="utf-8", errors="replace")[-8000:]
    if "MBDyn terminated normally" in tail:
        return True
    steps = re.findall(
        r"^Step\s+\d+\s+([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)",
        tail,
        flags=re.MULTILINE,
    )
    return bool(steps) and float(steps[-1]) >= 24.0 - 1.0e-9


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["velocity_mps", "density_ips", "status", "returncode", "elapsed_s"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run() -> int:
    args = arguments()
    input_path = args.input.expanduser().resolve()
    output_directory = args.output_dir.expanduser()
    if not output_directory.is_absolute():
        output_directory = (ROOT / output_directory).resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if args.rho <= 0.0:
        raise ValueError("rho deve essere positiva")

    mbdyn = None if args.dry_run else locate_mbdyn(args.mbdyn)
    velocities = velocity_grid(args.vmin, args.vmax, args.dv)
    template = input_path.read_text(encoding="utf-8")
    output_directory.mkdir(parents=True, exist_ok=True)

    print(f"Input: {input_path}")
    print(f"Risultati: {output_directory}")
    print(f"Densita': {args.rho:.8e} IPS")
    print("Velocita': " + ", ".join(f"{value:g}" for value in velocities))

    manifest: list[dict[str, object]] = []
    for index, velocity in enumerate(velocities, start=1):
        case_directory = output_directory / f"V_{velocity:05.1f}_mps"
        case_directory.mkdir(parents=True, exist_ok=True)
        archive_input = case_directory / "case_input.mbd"

        if args.resume and completed_case(case_directory):
            print(f"[{index:02d}/{len(velocities):02d}] {velocity:5.1f} m/s: gia' completo")
            manifest.append(
                dict(
                    velocity_mps=velocity,
                    density_ips=f"{args.rho:.12e}",
                    status="skipped_complete",
                    returncode=0,
                    elapsed_s=0.0,
                )
            )
            write_manifest(output_directory / "sweep_manifest.csv", manifest)
            continue

        archive_input.write_text(
            generated_source(template, velocity, args.rho), encoding="utf-8"
        )

        if args.dry_run:
            print(f"[{index:02d}/{len(velocities):02d}] {velocity:5.1f} m/s: input generato")
            manifest.append(
                dict(
                    velocity_mps=velocity,
                    density_ips=f"{args.rho:.12e}",
                    status="dry_run",
                    returncode=0,
                    elapsed_s=0.0,
                )
            )
            continue

        temporary_input = ROOT / f".__trim_V_{velocity:05.1f}.mbd"
        temporary_input.write_text(archive_input.read_text(encoding="utf-8"), encoding="utf-8")
        output_prefix = case_directory / "result"
        command = [str(mbdyn), "-f", temporary_input.name, "-o", str(output_prefix)]
        print(f"[{index:02d}/{len(velocities):02d}] {velocity:5.1f} m/s: avvio")
        started = time.monotonic()
        try:
            process = subprocess.run(
                command,
                cwd=ROOT,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
        finally:
            temporary_input.unlink(missing_ok=True)
        elapsed = time.monotonic() - started
        (case_directory / "console.log").write_text(
            process.stdout, encoding="utf-8", errors="replace"
        )
        status = "complete" if process.returncode == 0 else "failed"
        manifest.append(
            dict(
                velocity_mps=velocity,
                density_ips=f"{args.rho:.12e}",
                status=status,
                returncode=process.returncode,
                elapsed_s=f"{elapsed:.3f}",
            )
        )
        write_manifest(output_directory / "sweep_manifest.csv", manifest)
        print(f"    {status}, {elapsed:.1f} s")
        if process.returncode != 0 and not args.continue_on_error:
            print(f"Errore: vedere {case_directory / 'console.log'}", file=sys.stderr)
            return process.returncode or 1

    write_manifest(output_directory / "sweep_manifest.csv", manifest)
    if args.analyse and not args.dry_run:
        return subprocess.call(
            [sys.executable, str(ROOT / "analyze_trim_sweep.py"), "--results", str(output_directory)],
            cwd=ROOT,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
