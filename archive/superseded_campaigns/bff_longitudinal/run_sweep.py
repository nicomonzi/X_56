#!/usr/bin/env python3
"""Esegue lo sweep longitudinale MBDyn variando esclusivamente V_INF.

Per default esegue 57.5, 60.0, 62.5 e 65.0 m/s. Gli output sono piatti nella
cartella Windows ``C:\\Users\\Utente\\Desktop\\bbf_longitudinal`` e ogni
prefisso contiene la velocita'. La destinazione puo' essere cambiata con la
variabile d'ambiente ``BFF_OUTPUT_DIR``.
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
from decimal import Decimal
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "main_bbf.mbd"
DEFAULT_OUTPUT = Path("/mnt/c/Users/Utente/Desktop/bbf_longitudinal")
OUTPUT = Path(os.environ.get("BFF_OUTPUT_DIR", str(DEFAULT_OUTPUT)))
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")
RHO_IPS = 9.7284e-8
EXPECTED_FINAL_TIME = 15.0
EXPECTED_TIME_STEP = 0.005
CASE_SUFFIXES = (".mbd", ".nc", ".out", ".log", ".stdout", ".mov", ".usr", ".aer")


def constant(text: str, name: str) -> float:
    match = re.search(
        rf"(?m)^\s*set:\s*const\s+(?:real|integer)\s+{re.escape(name)}\s*=\s*([^;]+);",
        text,
    )
    if not match:
        raise RuntimeError(f"costante MBDyn mancante: {name}")
    expression = match.group(1).strip()
    try:
        return float(expression)
    except ValueError:
        return float(eval(expression, {"__builtins__": {}}, {"deg2rad": math.pi / 180.0, "pi": math.pi}))


def replace_real_constant(text: str, name: str, value: float) -> str:
    rendered, count = re.subn(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*[^;]+;",
        f"set: const real {name} = {value:.8f};",
        text,
    )
    if count != 1:
        raise RuntimeError(f"{name}: attesa una definizione, trovate {count}")
    return rendered


def case_stem(velocity: float) -> str:
    return f"V_{velocity:08.4f}".replace(".", "p")


def render_case(velocity: float) -> str:
    if not math.isfinite(velocity) or velocity <= 0.0:
        raise ValueError("ogni velocita' deve essere positiva e finita")
    text = replace_real_constant(MODEL.read_text(), "V_INF", velocity)
    if not math.isclose(constant(text, "RHO_AIR"), RHO_IPS, rel_tol=0.0, abs_tol=1e-15):
        raise RuntimeError("RHO_AIR deve restare 9.7284e-8 in unita' IPS")
    if not math.isclose(constant(text, "FINAL_TIME"), EXPECTED_FINAL_TIME, abs_tol=1e-12):
        raise RuntimeError("FINAL_TIME deve restare 15 s")
    if not math.isclose(constant(text, "TIME_STEP"), EXPECTED_TIME_STEP, abs_tol=1e-12):
        raise RuntimeError("TIME_STEP deve restare 0.005 s con 25 modi elastici")
    joint = re.search(r"(?ms)^\s*joint:\s*1,\s*total pin joint,(.*?);", text)
    compact = re.sub(r"\s+", " ", joint.group(1)) if joint else ""
    if "position constraint, active, active, inactive, null" not in compact:
        raise RuntimeError("il joint 1 deve vincolare solo X e Y")
    if "orientation constraint, inactive, inactive, inactive, null" not in compact:
        raise RuntimeError("roll, pitch e yaw devono restare liberi")
    forbidden = ("AIRSPEED_PID", "TRIM_THRUST", "THRUST_TOTAL_DRIVE", "force: 9500")
    if any(item in text for item in forbidden):
        raise RuntimeError("trovato thrust/airspeed hold non ammesso")
    surfaces = (ROOT / "INCLUDE/control_surfaces.mbd").read_text()
    wf1l = re.search(r"(?ms)# WF1L.*?component,(.*?);", surfaces)
    bfr = re.search(r"(?ms)# BFR.*?component,(.*?);", surfaces)
    if not wf1l or "const, TRIM_SURFACE" not in wf1l.group(1):
        raise RuntimeError("WF1L deve usare TRIM_SURFACE")
    if not bfr or "const, BODY_TRIM_SURFACE" not in bfr.group(1):
        raise RuntimeError("BFR deve usare BODY_TRIM_SURFACE")
    driver = (ROOT / "INCLUDE/driver.mbd").read_text()
    if "reference, PITCH_BFF_FREE_DRIVE, reference, Q_BFF_FREE_DRIVE" not in driver:
        raise RuntimeError("la banca band-stop longitudinale non e' collegata")
    modal = (ROOT / "INCLUDE/modaljoint.mbd").read_text()
    modal_code = "\n".join(line.split("#", 1)[0] for line in modal.splitlines())
    selected_block = re.search(r"(?s)\b25\s*,\s*list\s*,(.*?)\binitial value", modal_code)
    selected_modes = [int(value) for value in re.findall(r"\b\d+\b", selected_block.group(1))] if selected_block else []
    initial_modes = [int(value) for value in re.findall(r"mode,\s*(\d+)", modal_code)]
    if selected_modes != list(range(7, 32)) or initial_modes != list(range(7, 32)):
        raise RuntimeError("modaljoint deve selezionare i 25 modi elastici FEM 7--31")
    with (ROOT / "INCLUDE/mbdyn_modal.fem").open(errors="ignore") as stream:
        fem_header = "".join(stream.readline() for _ in range(10))
    if not re.search(r"REV0\s+8527\s+60\s+0\s+0\s+0", fem_header):
        raise RuntimeError("mbdyn_modal.fem non e' il FEM SOL103 a 60 modi atteso")
    return text.replace('"./INCLUDE/', f'"{ROOT / "INCLUDE"}/')


def decimal_sweep(start: float, stop: float, step: float) -> list[float]:
    a, b, h = Decimal(str(start)), Decimal(str(stop)), Decimal(str(step))
    if h <= 0 or b < a:
        raise ValueError("richiesti start <= stop e step > 0")
    values: list[float] = []
    value = a
    while value <= b + Decimal("1e-12"):
        values.append(float(value))
        value += h
    return values


def clear_case(stem: str) -> None:
    for suffix in CASE_SUFFIXES:
        path = OUTPUT / f"{stem}{suffix}"
        if path.is_file():
            path.unlink()


def run_case(velocity: float, overwrite: bool) -> int:
    stem = case_stem(velocity)
    input_path = OUTPUT / f"{stem}.mbd"
    prefix = OUTPUT / stem
    nc_path = OUTPUT / f"{stem}.nc"
    rendered = render_case(velocity)
    if nc_path.exists() and not overwrite:
        if not input_path.exists() or input_path.read_text() != rendered:
            raise RuntimeError(f"{stem}: output obsoleto; ripetere con --overwrite")
        print(f"[reuse] {stem}")
        return 0
    if overwrite:
        clear_case(stem)
    input_path.write_text(rendered)
    command = [str(MBDYN), "-s", "-f", str(input_path), "-o", str(prefix)]
    q_pa = 0.5 * 1.039663910516137 * velocity**2
    print(f"[run] {stem}: V={velocity:.4f} m/s, q={q_pa:.2f} Pa")
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    (OUTPUT / f"{stem}.stdout").write_text(result.stdout + result.stderr)
    if not nc_path.exists():
        raise RuntimeError(f"{stem}: MBDyn non ha prodotto {nc_path.name}")
    print(f"[done] {stem}: returncode={result.returncode}")
    return result.returncode


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=float, default=57.5, help="prima velocita' [m/s]")
    parser.add_argument("--stop", type=float, default=65.0, help="ultima velocita' [m/s]")
    parser.add_argument("--step", type=float, default=2.5, help="incremento [m/s]")
    parser.add_argument("--velocities", nargs="+", type=float, help="lista esplicita; sostituisce start/stop/step")
    parser.add_argument("--overwrite", action="store_true", help="sovrascrive i casi omonimi")
    parser.add_argument("--stop-on-error", action="store_true", help="arresta lo sweep al primo returncode non nullo")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    velocities = args.velocities or decimal_sweep(args.start, args.stop, args.step)
    failures: list[tuple[float, int]] = []
    for velocity in velocities:
        returncode = run_case(velocity, args.overwrite)
        if returncode:
            failures.append((velocity, returncode))
            if args.stop_on_error:
                break
    print(f"Completati {len(velocities) - len(failures)}/{len(velocities)} casi. Output: {OUTPUT}")
    if failures:
        print("Casi con returncode non nullo:", ", ".join(f"{v:g} ({code})" for v, code in failures))


if __name__ == "__main__":
    main()
