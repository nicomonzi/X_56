#!/usr/bin/env python3
"""Prepare/run an excited-only, shadow-referenced fictitious-softening sweep.

The already completed physical-prestress shadow trajectory is used as the
moving reference of the additional force.  Consequently the added force is
zero on that trajectory, while its perturbation Jacobian adds

    (TOTAL_PRESTRESS_GAIN - 1) * Delta_n * K_h,n

to the physical prestress stiffness already present in the source deck.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/fictitious_softening_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent
BASELINE_DIR = (ROOT / "../bff_open_loop").resolve()
SOURCE_ROOT = Path("/mnt/c/Users/Utente/Desktop/BFF_PULLUP_STIFFNESS_SWEEP")
DEFAULT_OUTPUT = Path("/mnt/c/Users/Utente/Desktop/BFF_PULLUP_FAKE_SOFTENING_X4")
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")

VELOCITIES_MPS = (66.25, 66.75, 67.25)
NOMINAL_LOAD_FACTORS = (1.3, 1.6)
TOTAL_PRESTRESS_GAIN = 4.0
EXTRA_PRESTRESS_GAIN = TOTAL_PRESTRESS_GAIN - 1.0

META_RE = re.compile(r"(?m)^# STIFFNESS_STUDY_METADATA (\{.*\})$")
REAL_RE = r"(?m)^set:\s*const\s+real\s+{name}\s*=\s*([-+0-9.eE]+)\s*;"

sys.path.insert(0, str(ROOT))
from compare_paired_response import growth_metrics, paired_delta  # noqa: E402


def close(a: float, b: float) -> bool:
    return abs(a - b) < 1e-8


def metadata(path: Path) -> dict:
    match = META_RE.search(path.read_text())
    if not match:
        raise ValueError(f"metadata mancanti: {path}")
    return json.loads(match.group(1))


def constant(path: Path, name: str) -> float:
    match = re.search(REAL_RE.format(name=re.escape(name)), path.read_text())
    if not match:
        raise ValueError(f"{name} non trovato in {path}")
    return float(match.group(1))


def source_catalog() -> dict[tuple[float, float, bool], tuple[Path, Path, dict]]:
    result: dict[tuple[float, float, bool], tuple[Path, Path, dict]] = {}
    for mbd in sorted((SOURCE_ROOT / "cases").glob("*.mbd")):
        meta = metadata(mbd)
        if meta.get("campaign") != "prestress_rom":
            continue
        key = (
            float(meta["velocity_mps"]),
            float(meta["nominal_load_factor"]),
            bool(meta["excited"]),
        )
        if key in result:
            raise RuntimeError(f"caso sorgente duplicato: {key}")
        result[key] = (mbd, mbd.with_suffix(".nc"), meta)
    return result


def find_source(
    catalog: dict[tuple[float, float, bool], tuple[Path, Path, dict]],
    velocity: float,
    load: float,
    excited: bool,
) -> tuple[Path, Path, dict]:
    matches = [
        value
        for (v, n, e), value in catalog.items()
        if close(v, velocity) and close(n, load) and e == excited
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"atteso un caso prestress_rom V={velocity:g}, n={load:g}, "
            f"excited={excited}; trovati {len(matches)}"
        )
    mbd, nc, meta = matches[0]
    if not nc.is_file():
        raise FileNotFoundError(nc)
    return mbd, nc, meta


def tag(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}".replace("-", "m").replace(".", "p")


def write_shadow_history(shadow_nc: Path, output: Path) -> tuple[int, float, float]:
    with Dataset(shadow_nc) as data:
        time = np.asarray(data["time"][:], dtype=float).squeeze()
        modal = np.asarray(data["elem.joint.5.a"][:], dtype=float)
    if time.ndim != 1 or modal.shape != (len(time), 6):
        raise ValueError(f"storia modale non valida: {shadow_nc}")
    if len(time) < 2 or not np.all(np.diff(time) > 0.0):
        raise ValueError(f"tempo non strettamente crescente: {shadow_nc}")
    with output.open("w") as stream:
        stream.write("# time q7 q8 q9 q10 q11 q12\n")
        for t, row in zip(time, modal):
            stream.write(
                " ".join([f"{t:.16e}", *(f"{value:.16e}" for value in row)])
                + "\n"
            )
    return len(time), float(time[0]), float(time[-1])


def fictitious_force() -> str:
    modes = range(7, 13)
    lines = [
        "    # Fictitious sensitivity force centered on the completed physical shadow.",
        "    # It is zero on the shadow and changes only the excited-shadow tangent stiffness.",
        "    force: FICTITIOUS_SOFTENING_FORCE, modal, MODAL_JOINT,",
        "        list, 6, 7, 8, 9, 10, 11, 12,",
    ]
    for mode_i in modes:
        lines.append("        array, 6,")
        terms = []
        for mode_j in modes:
            terms.extend([
                f'            element, MODAL_JOINT, joint, string, "q[{mode_j}]",',
                (
                    '                string, "-FICTITIOUS_EXTRA_PRESTRESS_GAIN*'
                    'model::drive(PRESTRESS_DELTA_N_DRIVE,Time)*'
                    f'KHN_{mode_i:02d}_{mode_j:02d}*'
                    f'(Var-model::drive(SHADOW_Q{mode_j:02d}_DRIVE,Time))*'
                    '((Time>=SAS_OFF_START)&&(Time<SAS_ON_START))",'
                ),
            ])
        lines.extend(terms)
    lines.extend(["        output, no;", ""])
    return "\n".join(lines)


def render_case(source: Path, history: Path, shadow_nc: Path) -> tuple[str, dict]:
    text = source.read_text()
    source_meta = metadata(source)
    meta = dict(source_meta)
    meta.update({
        "campaign": "fictitious_softening_x4",
        "fictitious_softening": True,
        "fictitious_model": "full_Khn_centered_on_recorded_physical_shadow",
        "total_prestress_gain": TOTAL_PRESTRESS_GAIN,
        "extra_prestress_gain": EXTRA_PRESTRESS_GAIN,
        "shadow_reused": True,
        "shadow_reference_nc": str(shadow_nc.resolve()),
    })
    text, count = META_RE.subn(
        "# STIFFNESS_STUDY_METADATA " + json.dumps(meta, sort_keys=True), text, count=1
    )
    if count != 1:
        raise RuntimeError("metadata sorgente non univoci")

    constants_marker = "set: const integer MODAL_JOINT = 5;"
    constants = f"""

# Amplified fictitious prestress: total perturbation stiffness uses
# TOTAL_PRESTRESS_GAIN*K_h,n. The source deck already contains gain 1.
set: const integer FICTITIOUS_SOFTENING_FORCE = 9701;
set: const integer SHADOW_MODAL_FILE_DRIVER = 14001;
set: const real FICTITIOUS_TOTAL_PRESTRESS_GAIN = {TOTAL_PRESTRESS_GAIN:.16e};
set: const real FICTITIOUS_EXTRA_PRESTRESS_GAIN = {EXTRA_PRESTRESS_GAIN:.16e};
set: const integer SHADOW_Q07_DRIVE = 14011;
set: const integer SHADOW_Q08_DRIVE = 14012;
set: const integer SHADOW_Q09_DRIVE = 14013;
set: const integer SHADOW_Q10_DRIVE = 14014;
set: const integer SHADOW_Q11_DRIVE = 14015;
set: const integer SHADOW_Q12_DRIVE = 14016;
"""
    if text.count(constants_marker) != 1:
        raise RuntimeError("MODAL_JOINT non univoco")
    text = text.replace(constants_marker, constants_marker + constants, 1)

    text, count = re.subn(r"(?m)^(\s*forces:\s*)3(\s*;)", r"\g<1>4\g<2>", text, count=1)
    if count != 1:
        raise RuntimeError("conteggio forces=3 non trovato")
    control_marker = "    aerodynamic elements: 58;"
    if text.count(control_marker) != 1:
        raise RuntimeError("conteggio aerodynamic elements non univoco")
    text = text.replace(control_marker, control_marker + "\n    file drivers: 1;", 1)

    nodes_end = "end: nodes;"
    drivers = f'''end: nodes;

begin: drivers;
    file: SHADOW_MODAL_FILE_DRIVER, variable step, 6,
        interpolation, linear,
        pad zeros, no,
        "{history.resolve()}";
end: drivers;'''
    if text.count(nodes_end) != 1:
        raise RuntimeError("end: nodes non univoco")
    text = text.replace(nodes_end, drivers, 1)

    marker = "    # Idealized active modal damper"
    if text.count(marker) != 1:
        raise RuntimeError("marker del damper SAS non univoco")
    drives = "\n".join(
        f"    drive caller: SHADOW_Q{mode:02d}_DRIVE, file, SHADOW_MODAL_FILE_DRIVER, {mode - 6};"
        for mode in range(7, 13)
    )
    insertion = drives + "\n\n" + fictitious_force() + "\n"
    text = text.replace(marker, insertion + marker, 1)

    required = (
        "forces: 4;",
        "file drivers: 1;",
        "begin: drivers;",
        "force: PRESTRESS_ROM_FORCE, modal, MODAL_JOINT",
        "force: FICTITIOUS_SOFTENING_FORCE, modal, MODAL_JOINT",
        "FICTITIOUS_EXTRA_PRESTRESS_GAIN",
        "Var-model::drive(SHADOW_Q07_DRIVE,Time)",
        "((Time>=SAS_OFF_START)&&(Time<SAS_ON_START))",
    )
    if not all(item in text for item in required):
        raise RuntimeError("invarianti del caso fittizio non soddisfatte")
    return text, meta


def input_hash(text: str, history: Path) -> str:
    digest = hashlib.sha256(text.encode())
    digest.update(history.read_bytes())
    return digest.hexdigest()


def write_manifests(output: Path, rows: list[dict]) -> None:
    with (output / "manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "manifest.json").write_text(json.dumps({
        "description": "Excited-only fictitious full-matrix softening sweep",
        "equation": (
            "K_pert=K_fem+4*Delta_n*K_h_n; extra force="
            "-3*Delta_n*K_h_n*(q-q_shadow(t)) during SAS-off"
        ),
        "velocities_mps": VELOCITIES_MPS,
        "nominal_load_factors": NOMINAL_LOAD_FACTORS,
        "total_prestress_gain": TOTAL_PRESTRESS_GAIN,
        "source_campaign": str(SOURCE_ROOT),
        "cases": rows,
    }, indent=2) + "\n")


def prepare(output: Path, overwrite: bool) -> list[dict]:
    output.mkdir(parents=True, exist_ok=True)
    cases = output / "cases"
    histories = output / "shadow_histories"
    cases.mkdir(exist_ok=True)
    histories.mkdir(exist_ok=True)
    catalog = source_catalog()
    rows: list[dict] = []
    index = 0
    for load in NOMINAL_LOAD_FACTORS:
        for velocity in VELOCITIES_MPS:
            index += 1
            shadow_mbd, shadow_nc, _ = find_source(catalog, velocity, load, False)
            excited_mbd, excited_nc, _ = find_source(catalog, velocity, load, True)
            stem = f"fake_softening_x4__V_{tag(velocity)}__n_{tag(load)}__excited"
            history = histories / f"shadow__V_{tag(velocity)}__n_{tag(load)}.dat"
            samples, history_start, history_end = write_shadow_history(shadow_nc, history)
            rendered, _ = render_case(excited_mbd, history, shadow_nc)
            target = cases / f"{stem}.mbd"
            target_nc = target.with_suffix(".nc")
            digest = input_hash(rendered, history)
            if target_nc.is_file() and target.is_file() and target.read_text() != rendered and not overwrite:
                raise RuntimeError(f"input cambiato ma risultato esistente: {target_nc}")
            target.write_text(rendered)
            complete = result_is_complete(target_nc, target)
            rows.append({
                "case_index": index,
                "velocity_mps": velocity,
                "nominal_load_factor": load,
                "total_prestress_gain": TOTAL_PRESTRESS_GAIN,
                "extra_prestress_gain": EXTRA_PRESTRESS_GAIN,
                "stem": stem,
                "input_mbd": str(target.resolve()),
                "result_nc": str(target_nc.resolve()),
                "shadow_mbd": str(shadow_mbd.resolve()),
                "shadow_nc": str(shadow_nc.resolve()),
                "baseline_excited_nc": str(excited_nc.resolve()),
                "shadow_history": str(history.resolve()),
                "shadow_history_samples": samples,
                "shadow_history_start_s": history_start,
                "shadow_history_end_s": history_end,
                "input_and_history_sha256": digest,
                "status": "existing" if complete and not overwrite else "prepared",
                "message": "complete result retained" if complete and not overwrite else "excited input prepared",
            })
    write_manifests(output, rows)
    (output / "README.txt").write_text(
        "X-56 fictitious softening x4\n\n"
        "Six new excited-only cases: V=66.25/66.75/67.25 m/s and n=1.3/1.6.\n"
        "Existing physical-prestress shadows are reused as moving references.\n"
        "K_pert = K_fem + 4 Delta_n K_h,n during SAS-off.\n\n"
        "Run from /home/nicomonzi/X_56/workflows/maneuver_bff:\n"
        "python3 fictitious_softening_sweep.py --execute --jobs 2 --analyse\n"
    )
    return rows


def result_is_complete(nc: Path, mbd: Path) -> bool:
    if not nc.is_file():
        return False
    try:
        required = constant(mbd, "SAS_OFF_START") + constant(mbd, "SAS_OFF_DURATION")
        with Dataset(nc) as data:
            time = np.asarray(data["time"][:], dtype=float).squeeze()
        return time.ndim == 1 and len(time) >= 50 and float(time[-1]) >= required - 0.011
    except Exception:
        return False


def run_one(row: dict, overwrite: bool) -> tuple[str, str]:
    mbd = Path(row["input_mbd"])
    prefix = mbd.with_suffix("")
    if result_is_complete(prefix.with_suffix(".nc"), mbd) and not overwrite:
        return "reused", "existing complete result retained"
    if overwrite:
        for suffix in (".nc", ".out", ".log", ".stdout", ".mov"):
            prefix.with_suffix(suffix).unlink(missing_ok=True)
    result = subprocess.run(
        [str(MBDYN), "-s", "-f", str(mbd), "-o", str(prefix)],
        cwd=BASELINE_DIR,
        text=True,
        capture_output=True,
    )
    prefix.with_suffix(".stdout").write_text(result.stdout + result.stderr)
    if result.returncode != 0 or not result_is_complete(prefix.with_suffix(".nc"), mbd):
        return "failed", f"returncode={result.returncode}; vedere {prefix.with_suffix('.stdout')}"
    return "complete", "MBDyn completed"


def onset(points: list[dict]) -> dict:
    points = sorted(points, key=lambda row: row["velocity_mps"])
    for left, right in zip(points, points[1:]):
        sl, sr = left["sigma_per_s"], right["sigma_per_s"]
        if sl == 0.0 or sl * sr <= 0.0:
            vl, vr = left["velocity_mps"], right["velocity_mps"]
            value = vl if sl == 0.0 else vl - sl * (vr - vl) / (sr - sl)
            return {"status": "crossing_found", "bracket_mps": [vl, vr], "onset_mps": value}
    return {"status": "no_bracket", "bracket_mps": None, "onset_mps": None}


def analyse(output: Path, rows: list[dict]) -> None:
    analysis = output / "analysis"
    analysis.mkdir(exist_ok=True)
    results: list[dict] = []
    missing: list[str] = []
    for row in rows:
        mbd = Path(row["input_mbd"])
        shadow = Path(row["shadow_nc"])
        for gain, excited in (
            (1.0, Path(row["baseline_excited_nc"])),
            (TOTAL_PRESTRESS_GAIN, Path(row["result_nc"])),
        ):
            if not excited.is_file():
                missing.append(str(excited))
                continue
            time, delta = paired_delta(shadow, excited)
            velocity = float(row["velocity_mps"])
            frequency = float(np.interp(
                velocity,
                [65.0, 67.5],
                [2.0551998, 2.0584934],
            ))
            release = constant(mbd, "SAS_OFF_START")
            duration = constant(mbd, "SAS_OFF_DURATION")
            fit_start = release + 0.05 + 0.742 / frequency + 0.05
            fit_end = release + duration - 0.05
            metrics, _, _ = growth_metrics(time, delta, fit_start, fit_end, frequency)
            sigma = float(metrics["sigma_per_s"])
            results.append({
                "velocity_mps": velocity,
                "nominal_load_factor": float(row["nominal_load_factor"]),
                "total_prestress_gain": gain,
                "sigma_per_s": sigma,
                "verdict": "unstable" if sigma > 0.05 else "stable" if sigma < -0.05 else "near_onset",
                "fit_start_s": fit_start,
                "fit_end_s": fit_end,
                "energy_growth_ratio": metrics["energy_growth_ratio"],
                "shadow_nc": str(shadow),
                "excited_nc": str(excited),
            })
    if not results:
        raise RuntimeError("nessun risultato completo da analizzare")
    with (analysis / "growth_rates.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    onsets = []
    for load in NOMINAL_LOAD_FACTORS:
        for gain in (1.0, TOTAL_PRESTRESS_GAIN):
            points = [
                row for row in results
                if close(row["nominal_load_factor"], load)
                and close(row["total_prestress_gain"], gain)
            ]
            if points:
                onsets.append({
                    "nominal_load_factor": load,
                    "total_prestress_gain": gain,
                    "point_count": len(points),
                    **onset(points),
                })
    with (analysis / "onset.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(onsets[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(onsets)
    (analysis / "summary.json").write_text(json.dumps({
        "equation": "K_pert=K_fem+gain*Delta_n*K_h_n",
        "new_excited_case_count": len(rows),
        "available_growth_result_count": len(results),
        "missing_results": missing,
        "onsets": onsets,
    }, indent=2) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    colors = {1.3: "tab:blue", 1.6: "tab:red"}
    for load in NOMINAL_LOAD_FACTORS:
        for gain, style in ((1.0, "o--"), (TOTAL_PRESTRESS_GAIN, "o-")):
            points = sorted(
                (row for row in results if close(row["nominal_load_factor"], load)
                 and close(row["total_prestress_gain"], gain)),
                key=lambda row: row["velocity_mps"],
            )
            if points:
                axes[0].plot(
                    [row["velocity_mps"] for row in points],
                    [row["sigma_per_s"] for row in points],
                    style,
                    color=colors[load],
                    alpha=0.6 if gain == 1.0 else 1.0,
                    label=f"n={load:g}, gain Khn={gain:g}",
                )
    for gain, style in ((1.0, "o--"), (TOTAL_PRESTRESS_GAIN, "o-")):
        available = [row for row in onsets if row["total_prestress_gain"] == gain and row["onset_mps"] is not None]
        if available:
            axes[1].plot(
                [row["nominal_load_factor"] for row in available],
                [row["onset_mps"] for row in available],
                style,
                label=f"gain Khn={gain:g}",
            )
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set(xlabel="TAS [m/s]", ylabel="sigma [1/s]", title="BFF growth: physical vs fictitious softening")
    axes[1].set(xlabel="nominal load factor", ylabel="onset TAS [m/s]", title="Flutter onset shift")
    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(analysis / "summary.png", dpi=190)
    plt.close(fig)
    print(f"[analysis] {analysis}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--execute", action="store_true", help="run the six new excited cases")
    parser.add_argument("--analyse", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs deve essere almeno 1")
    output = args.output.expanduser().resolve()
    rows = prepare(output, args.overwrite)
    print(f"[prepared] {len(rows)} excited in {output / 'cases'}")
    print(f"[reused shadows] {len(rows)} from {SOURCE_ROOT / 'cases'}")
    if args.execute:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(run_one, row, args.overwrite): i for i, row in enumerate(rows)}
            for done, future in enumerate(as_completed(futures), start=1):
                index = futures[future]
                status, message = future.result()
                rows[index]["status"] = status
                rows[index]["message"] = message
                write_manifests(output, rows)
                print(f"[{done}/{len(rows)}] {rows[index]['stem']}: {status}", flush=True)
        if any(row["status"] == "failed" for row in rows):
            raise SystemExit(2)
    else:
        print("[no-run] usare --execute per avviare MBDyn")
    if args.analyse:
        analyse(output, rows)


if __name__ == "__main__":
    main()
