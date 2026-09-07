#!/usr/bin/env python3
"""Prepare the two X-56 SOL103/STATSUB prestressed modal exports.

The script never launches Nastran. With ``--overwrite`` it removes only the
obsolete QRMETH/SOL106 campaign and known outputs of the generated SOL103
jobs, so a stale result cannot be mistaken for a new one.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOL103_STEM = "prestress_sol103_statsub"
ALTER_NAME = "MBDyn_NASTRAN_alter_SOL103_2024.nas"
SOLVER_OUTPUT_SUFFIXES = (".f04", ".f06", ".log", ".op2", ".fem")
OBSOLETE_CASE_FILES = (
    "MBDyn_NASTRAN_alter_SOL106_2024.nas",
    "prestress_sol106_diagnostic.bdf",
    "prestress_sol106_diagnostic.f04",
    "prestress_sol106_diagnostic.f06",
    "prestress_sol106_diagnostic.log",
    "prestress_sol106_diagnostic.op2",
)


def case_tag(load_factor: float) -> str:
    return f"n{load_factor:.3f}".replace(".", "p")


def has_param(bulk: str, name: str) -> bool:
    """Accept both free-field and fixed-field PARAM entries."""
    return re.search(
        rf"(?mi)^\s*PARAM(?:\s*,\s*|\s+){re.escape(name)}(?:\s*,|\s+)",
        bulk,
    ) is not None


def parameter_names(bulk: str) -> list[str]:
    names = []
    for line in bulk.splitlines():
        if line.lstrip().startswith("$"):
            continue
        match = re.match(r"(?i)^\s*PARAM(?:\s*,\s*|\s+)([A-Z0-9_]+)(?:\s*,|\s+)", line)
        if match:
            names.append(match.group(1).upper())
    return names


def convert_deck(text: str, load_factor: float, achieved_nz: float) -> str:
    """Add the MBDyn SOL103 ALTER without changing the validated bulk model."""
    if not re.search(r"(?mi)^\s*SOL\s+103\s*$", text):
        raise ValueError("SOL 103 non trovato nel deck sorgente")
    if not re.search(r"(?mi)^\s*STATSUB\s*=\s*1\s*$", text):
        raise ValueError("STATSUB=1 non trovato nel deck sorgente")
    if "BEGIN BULK" not in text:
        raise ValueError("BEGIN BULK non trovato nel deck sorgente")

    # The ALTER contains SOL 103 itself and must replace the executive SOL line.
    deck = re.sub(
        r"(?mi)^\s*SOL\s+103\s*$",
        f"INCLUDE '{ALTER_NAME}'",
        text,
        count=1,
    )
    deck = re.sub(
        r"(?mi)^TITLE=.*$",
        f"TITLE=X56 PRESTRESS SOL103 STATSUB NZ={load_factor:g}",
        deck,
        count=1,
    )
    deck = deck.replace(
        "SUBTITLE=CG-SUPPORTED NORMAL MODES ABOUT PRELOADED STATE",
        "SUBTITLE=FREE-FREE NORMAL MODES ABOUT PRELOADED STATE",
        1,
    )
    deck, vector_count = re.subn(
        r"(?mi)^\s*VECTOR\([^\n]*\)\s*=\s*901\s*$",
        "    VECTOR(SORT1,REAL,PLOT)=ALL",
        deck,
        count=1,
    )
    if vector_count != 1:
        raise ValueError("selezione VECTOR=901 del subcase modale non trovata")
    case_control_end = deck.index("BEGIN BULK")
    subcase2 = re.search(r"(?mi)^\s*SUBCASE\s+2\s*$", deck[:case_control_end])
    if subcase2 is None:
        raise ValueError("SUBCASE 2 modale non trovato")
    modal_control = deck[subcase2.start():case_control_end]
    modal_control, spc_count = re.subn(
        r"(?mi)^\s*SPC\s*=\s*620\s*\n",
        "    $ No SPC here: prestressed modes must remain free-free\n",
        modal_control,
        count=1,
    )
    if spc_count != 1:
        raise ValueError("SPC=620 del subcase modale non trovato")
    deck = deck[:subcase2.start()] + modal_control + deck[case_control_end:]
    header = (
        "$ Prepared by maneuver_stiffness_rom/nastran/prepare_cases.py\n"
        f"$ Nominal n_z={load_factor:g}; recovered mean body-z specific n={achieved_nz:.10g}.\n"
        "$ Statically supported preload followed by free-free SOL103 modes with STATSUB=1.\n"
        "$ Original CQUADR/CTRIAR formulation retained; all mode-shape grids exported.\n"
        "$ Follower stiffness excluded.\n"
    )
    first_newline = deck.find("\n")
    deck = deck[: first_newline + 1] + header + deck[first_newline + 1 :]

    begin = deck.index("BEGIN BULK") + len("BEGIN BULK")
    bulk = deck[begin:]
    if not has_param(bulk, "FOLLOWK"):
        deck = deck[:begin] + "\nPARAM,FOLLOWK,NO" + deck[begin:]
    names = parameter_names(deck[begin:])
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"PARAM duplicati nel deck generato: {', '.join(duplicates)}")
    if "QRMETH=3" in deck.upper() or re.search(r"(?mi)^\s*SOL\s+106\s*$", deck):
        raise ValueError("il deck SOL103 contiene impostazioni obsolete SOL106/QRMETH")
    return deck


def remove_known_outputs(target: Path) -> None:
    for suffix in SOLVER_OUTPUT_SUFFIXES:
        (target / f"{SOL103_STEM}{suffix}").unlink(missing_ok=True)
    (target / "mbdyn_modal.mat").unlink(missing_ok=True)


def runner_text() -> str:
    return f'''#!/usr/bin/env bash
set -euo pipefail
solver="${{NASTRAN_CMD:-nast}}"
stem={SOL103_STEM}
rm -f -- "$stem.f04" "$stem.f06" "$stem.log" "$stem.op2" mbdyn_modal.mat
"$solver" old=no "$stem.bdf"
test -s "$stem.f06" || {{ echo "F06 mancante o vuoto" >&2; exit 2; }}
if grep -qE "USER FATAL MESSAGE|SYSTEM FATAL MESSAGE|FATAL ERROR" "$stem.f06"; then
    echo "Nastran ha scritto un errore fatale in $stem.f06" >&2
    exit 3
fi
grep -q "END OF JOB" "$stem.f06" || {{ echo "END OF JOB non trovato" >&2; exit 4; }}
test -s "$stem.op2" || {{ echo "OP2 mancante o vuoto" >&2; exit 5; }}
test -s mbdyn_modal.mat || {{ echo "mbdyn_modal.mat mancante o vuoto" >&2; exit 6; }}
echo "SOL103/STATSUB completato: F06, OP2 e mbdyn_modal.mat presenti."
'''


def prepare(config_path: Path, output: Path, overwrite: bool) -> list[Path]:
    config_path = config_path.expanduser().resolve()
    config = json.loads(config_path.read_text())
    alter = Path(config["nastran"]["alter_sol103"]).expanduser().resolve()
    if not alter.is_file():
        raise FileNotFoundError(alter)
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    obsolete_linear = output / "linear_qrmeth3"
    if overwrite and obsolete_linear.exists():
        shutil.rmtree(obsolete_linear)

    made = []
    for entry in config["nastran"]["source_cases"]:
        load_factor = float(entry["load_factor"])
        source = Path(entry["directory"]).expanduser().resolve()
        load_summary = json.loads((source / "load_summary.json").read_text())
        achieved_nz = float(load_summary["specific_acceleration_body_g_mean"][2])
        target = output / case_tag(load_factor)
        if target.exists() and any(target.iterdir()) and not overwrite:
            raise FileExistsError(f"{target} non vuota; usare --overwrite")
        target.mkdir(parents=True, exist_ok=True)
        if overwrite:
            remove_known_outputs(target)
            for name in OBSOLETE_CASE_FILES:
                (target / name).unlink(missing_ok=True)

        for name in ("preload_loads.bdf", "rbe3s.bdf", "load_summary.json"):
            shutil.copy2(source / name, target / name)
        if target.joinpath("BULK").exists() and overwrite:
            shutil.rmtree(target / "BULK")
        if not target.joinpath("BULK").exists():
            shutil.copytree(source / "BULK", target / "BULK")
        shutil.copy2(alter, target / ALTER_NAME)

        source_deck = (source / "prestressed_modes_nodlm.bdf").read_text()
        deck = target / f"{SOL103_STEM}.bdf"
        deck.write_text(convert_deck(source_deck, load_factor, achieved_nz))
        runner = target / "RUN_BY_USER.sh"
        runner.write_text(runner_text())
        runner.chmod(0o755)
        manifest = {
            "purpose": "fixed-basis geometric-stiffness snapshot",
            "load_factor": load_factor,
            "load_factor_definition": "nominal_maneuver_class",
            "achieved_body_z_specific_load_factor": achieved_nz,
            "source": str(source / "prestressed_modes_nodlm.bdf"),
            "deck": deck.name,
            "solution": "supported static preload; free-free SOL103 cold start with STATSUB=1",
            "shell_formulation": "original CQUADR/CTRIAR; no QRMETH conversion",
            "follower_force_stiffness": False,
            "run_by_codex": False,
            "expected_outputs": [
                f"{SOL103_STEM}.f06",
                f"{SOL103_STEM}.op2",
                "mbdyn_modal.mat",
            ],
            "warning": "KHH is in the snapshot eigenbasis and must be projected into fixed Phi0",
        }
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        made.append(target)
    (output / "README.md").write_text(
        "# X-56 SOL103/STATSUB exports\n\n"
        "Run only `n1p000/RUN_BY_USER.sh` and `n1p600/RUN_BY_USER.sh`. "
        "Each job must produce a non-empty F06, OP2 and `mbdyn_modal.mat`.\n"
    )
    return made


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config.json")
    parser.add_argument("--output", type=Path, help="override della destinazione configurata in ZENO")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.expanduser().resolve().read_text())
    output = args.output or Path(config["nastran"]["case_output_directory"])
    for path in prepare(args.config, output, args.overwrite):
        print(path)


if __name__ == "__main__":
    main()
