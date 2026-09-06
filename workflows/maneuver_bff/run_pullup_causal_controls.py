#!/usr/bin/env python3
"""Three no-rap controls completing a maneuver x SAS-release comparison.

Reuse primary V67.25 n1.6 shadow as the fourth cell. No Nastran calls.
New results only; refuse to overwrite existing inputs/results.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import re
import signal
import subprocess
import threading
import time

from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / "bff_open_loop"
SOURCE = Path("/mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/primary/cases")
GATE = "((Time<SAS_OFF_START)||(Time>=SAS_ON_START))"
STOP = threading.Event()
PROCESSES = set()
PROCESS_LOCK = threading.Lock()


def stop_children(signum, frame):
    STOP.set()
    with PROCESS_LOCK:
        for process in tuple(PROCESSES):
            if process.poll() is None:
                process.terminate()
    print(f"Stop requested ({signum}); no further cases will start.", flush=True)


def replace_constant(text, name, value):
    text, count = re.subn(rf"(?m)^set: const real {name} = [^;]+;",
                          f"set: const real {name} = {value};", text)
    if count != 1:
        raise ValueError((name, count))
    return text


def prepare(output):
    candidates = list(SOURCE.glob("*V_067p250*nnom_01p60_shadow.mbd"))
    if len(candidates) != 1:
        raise RuntimeError(candidates)
    source = candidates[0]
    original = source.read_text()
    if original.count(GATE) != 2:
        raise RuntimeError("Unexpected SAS gates; inspect rather than silently modify")
    metadata = json.loads(re.search(r"(?m)^# STIFFNESS_STUDY_METADATA (.+)$", original)[1])
    release = float(re.search(r"set: const real SAS_OFF_START = ([^;]+);", original)[1])
    final = float(re.search(r"set: const real FINAL_TIME = ([^;]+);", original)[1])
    rows = []
    for name, maneuver, release_sas in (("sham_release", False, True),
                                       ("pullup_sas_continuous", True, False),
                                       ("sham_sas_continuous", False, False)):
        text = original
        if not maneuver:
            for label in ("PULLUP_PITCH_ANGLE", "DIVE_PITCH_RATE_COMMAND"):
                text = replace_constant(text, label, "0.")
            text = replace_constant(text, "DIVE_NOMINAL_LOAD_CLASS", "1.")
        if not release_sas:
            text = text.replace(GATE, "1.")
        # Same integration and maneuver clock; save from before the dive.
        text = text.replace("closest next, SAS_OFF_START - STUDY_OUTPUT_PREHISTORY,",
                            "closest next, 7.,")
        meta = dict(metadata, campaign="causal_controls", excited=False,
                    nominal_load_factor=1.6 if maneuver else 1.0,
                    pitch_amplitude_deg=metadata["pitch_amplitude_deg"] if maneuver else 0.,
                    pitch_rate_command_deg_s=metadata["pitch_rate_command_deg_s"] if maneuver else 0.,
                    modification=name)
        text = re.sub(r"(?m)^# STIFFNESS_STUDY_METADATA .+$",
                      "# STIFFNESS_STUDY_METADATA " + json.dumps(meta, sort_keys=True), text)
        # Remove stale generator metadata that describes the unmodified command.
        text = re.sub(r"(?m)^# MANEUVER_METADATA .+\n", "", text)
        info = dict(name=name, maneuver=maneuver, sas_release=release_sas,
                    reference_release_s=release, required_final_s=final,
                    source=str(source), source_sha256=hashlib.sha256(original.encode()).hexdigest(),
                    note="No rap; DLM and structural stiffness unchanged; sham preserves all event timings.")
        text = "# CAUSAL_CONTROL_METADATA " + json.dumps(info, sort_keys=True) + "\n" + text
        case = output / "cases" / name
        case.parent.mkdir(parents=True, exist_ok=True)
        path = case.with_suffix(".mbd")
        if path.exists() and path.read_text() != text:
            raise RuntimeError(f"Refusing changed input overwrite: {path}")
        if not path.exists():
            path.write_text(text)
        rows.append(dict(info, prefix=str(case), input_sha256=hashlib.sha256(text.encode()).hexdigest()))
    manifest = output / "design.json"
    payload = dict(reused_pullup_release=str(source.with_suffix('.nc')), cases=rows,
                   interpretation="Separates maneuver history and SAS-release protocol. Continuous-SAS cells are closed-loop responses, not open-loop flutter tests.")
    serialized = json.dumps(payload, indent=2) + "\n"
    if manifest.exists() and manifest.read_text() != serialized:
        raise RuntimeError(f"Refusing manifest overwrite: {manifest}")
    if not manifest.exists():
        manifest.write_text(serialized)
    readme = output / "README.md"
    if not readme.exists():
        readme.write_text((ROOT / 'SWEEP_PULLUP_BFF.md').read_text())
    return rows


def run(row):
    if STOP.is_set():
        return dict(name=row['name'], complete=False, cancelled_before_start=True)
    prefix = Path(row["prefix"])
    status = prefix.with_suffix(".status.json")
    if status.exists():
        saved = json.loads(status.read_text())
        if saved.get("complete"):
            return saved
        raise RuntimeError(f"Existing unsuccessful case requires inspection: {prefix}")
    if any(prefix.with_suffix(s).exists() for s in (".nc", ".out", ".log", ".stdout")):
        raise RuntimeError(f"Refusing overwrite of existing solver output: {prefix}")
    start = time.monotonic()
    with prefix.with_suffix(".stdout").open("x") as stream:
        with PROCESS_LOCK:
            if STOP.is_set():
                return dict(name=row['name'], complete=False, cancelled_before_start=True)
            process = subprocess.Popen(["/usr/local/mbdyn/bin/mbdyn", "-s", "-f", str(prefix.with_suffix('.mbd')),
                                        "-o", str(prefix)], cwd=BASE, stdout=stream, stderr=subprocess.STDOUT)
            PROCESSES.add(process)
        returncode = process.wait()
        with PROCESS_LOCK:
            PROCESSES.discard(process)
    last = None
    if prefix.with_suffix('.nc').exists():
        try:
            with Dataset(prefix.with_suffix('.nc')) as nc:
                if len(nc['time']):
                    last = float(nc['time'][-1])
        except (OSError, RuntimeError):
            pass  # A terminated solver may leave an unreadable NetCDF.
    outcome = dict(name=row['name'], returncode=returncode,
                   wall_seconds=time.monotonic()-start, last_saved_time_s=last,
                   complete=returncode == 0 and last is not None and last >= row['required_final_s']-.021)
    status.write_text(json.dumps(outcome, indent=2)+'\n')
    return outcome


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=Path('/home/nicomonzi/ZENO/BFF_PULLUP_CAUSAL_READY'))
    p.add_argument('--execute', action='store_true')
    p.add_argument('--jobs', type=int, default=2)
    args = p.parse_args()
    if args.jobs < 1:
        p.error('--jobs must be positive')
    rows = prepare(args.output.resolve())
    print(f"Prepared {len(rows)} controls in {args.output}", flush=True)
    if args.execute:
        signal.signal(signal.SIGINT, stop_children)
        signal.signal(signal.SIGTERM, stop_children)
        outcomes=[]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for future in as_completed([pool.submit(run,row) for row in rows]):
                outcome=future.result()
                outcomes.append(outcome)
                print(json.dumps(outcome), flush=True)
        if STOP.is_set() or not all(row.get('complete') for row in outcomes):
            raise SystemExit(1)


if __name__ == '__main__':
    main()
