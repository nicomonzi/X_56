#!/usr/bin/env python3
"""Run a compact velocity sweep around the observed BBF-sensitive region."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parent
MAIN = ROOT / "main_bbf.mbd"
ANALYZER = ROOT / "analyze_sweep.py"
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")

# MODIFICA: griglia compatta. La risoluzione e' 0.1 m/s nella vecchia
# transizione 34.0--34.25 m/s; pochi punti esterni verificano che l'andamento
# sia fisico e monotono invece di una singola perdita di controllo.
DEFAULT_VELOCITIES = (33.75, 34.0, 34.1, 34.2, 34.3, 34.5, 35.0, 36.0)
DEFAULT_OUTPUT = Path(
    "/mnt/c/Users/Utente/Desktop/RESULTS_BBF_DT002_MODES_7_23_COMPACT"
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the stabilized-trim model near the BBF-sensitive velocity "
            "range and then generate the V-g/V-f analysis."
        )
    )
    parser.add_argument(
        "--velocities",
        type=float,
        nargs="+",
        default=DEFAULT_VELOCITIES,
        help="Velocity samples in m/s; default: 33.75 34 34.1 34.2 34.3 34.5 35 36.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Result root; default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--no-analysis",
        action="store_true",
        help="Run MBDyn only and skip analyze_sweep.py.",
    )
    return parser.parse_args()


def case_directory(output: Path, velocity: float) -> Path:
    token = f"{velocity:07.3f}".replace(".", "p")
    return output / f"V_{token}_mps"


def input_at_velocity(source: str, velocity: float) -> str:
    pattern = r"(?m)^set:\s*const\s+real\s+V_INF\s*=\s*[^;]+;"
    generated, replacements = re.subn(
        pattern,
        f"set: const real V_INF = {velocity:.6f};",
        source,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("Could not replace the unique V_INF definition")
    return generated


def preserve_input(case: Path, generated: str) -> Path:
    case.mkdir(parents=True, exist_ok=True)
    permanent = case / "case_input.mbd"
    if permanent.exists():
        if permanent.read_text(encoding="utf-8") != generated:
            raise RuntimeError(
                "Existing case was produced by a different model and will "
                f"not be overwritten: {case}"
            )
    else:
        permanent.write_text(generated, encoding="utf-8")
    return permanent


def result_is_complete(case: Path) -> bool:
    """Accept a NetCDF result only after the complete BBF window is present."""
    result = case / "case.nc"
    log = case / "case.log"
    if not result.is_file() or not log.is_file():
        return False
    match = re.search(
        r"(?m)^\s*const real BFF_WINDOW_END\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*$",
        log.read_text(encoding="utf-8", errors="replace"),
    )
    if match is None:
        return False
    required_time = float(match.group(1).replace("D", "E"))
    try:
        with Dataset(result) as dataset:
            return bool(
                len(dataset.variables["time"]) > 0
                and float(dataset.variables["time"][-1]) >= required_time
            )
    except (OSError, KeyError, IndexError):
        return False


def run_case(output: Path, source: str, velocity: float) -> None:
    case = case_directory(output, velocity)
    result = case / "case.nc"
    generated = input_at_velocity(source, velocity)

    if result_is_complete(case):
        preserve_input(case, generated)
        print(f"Reuse complete case at {velocity:g} m/s: {result}", flush=True)
        return
    if case.exists() and any(case.iterdir()):
        raise RuntimeError(
            "Incomplete non-empty case is preserved; move it before retrying: "
            f"{case}"
        )

    preserve_input(case, generated)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=".run_bbf_boundary_",
        suffix=".mbd",
        dir=ROOT,
        delete=False,
    ) as temporary:
        temporary.write(generated)
        active_input = Path(temporary.name)

    print(f"\n=== BBF boundary test: V_INF = {velocity:g} m/s ===", flush=True)
    try:
        subprocess.run(
            [
                str(MBDYN),
                "-f",
                str(active_input),
                "-o",
                str(case / "case"),
            ],
            cwd=ROOT,
            check=True,
        )
    finally:
        active_input.unlink(missing_ok=True)
    if not result_is_complete(case):
        raise RuntimeError(
            "MBDyn stopped before completing the BBF observation window; "
            f"the partial files were preserved in {case}"
        )


def main() -> None:
    args = arguments()
    output = args.output.expanduser().resolve()
    velocities = sorted(set(args.velocities))

    if not MAIN.is_file() or not MBDYN.is_file():
        raise SystemExit("main_bbf.mbd or the MBDyn executable is missing")
    if any(value <= 0.0 for value in velocities):
        raise SystemExit("Every velocity must be positive")

    source = MAIN.read_text(encoding="utf-8")
    output.mkdir(parents=True, exist_ok=True)

    # MODIFICA: questa copia documenta il modello comune a tutte le nuove run.
    model_copy = output / "main_bbf_model_used.mbd"
    if model_copy.exists() and model_copy.read_text(encoding="utf-8") != source:
        raise RuntimeError(
            "The output directory already belongs to a different model: "
            f"{output}"
        )
    if not model_copy.exists():
        model_copy.write_text(source, encoding="utf-8")

    print(
        "Velocity grid [m/s]: " + ", ".join(f"{v:g}" for v in velocities),
        flush=True,
    )
    print(f"Results: {output}", flush=True)
    for velocity in velocities:
        run_case(output, source, velocity)

    if not args.no_analysis:
        subprocess.run(
            [
                sys.executable,
                str(ANALYZER),
                "--results",
                str(output),
                "--output",
                str(output / "ANALYSIS"),
            ],
            cwd=ROOT,
            check=True,
        )


if __name__ == "__main__":
    main()
