#!/usr/bin/env python3
"""Costruisce la repository X_56 riorganizzata partendo da TESI.

La sorgente non viene modificata. Per sicurezza la destinazione deve contenere
soltanto ``.git``; lo script non sovrascrive una migrazione esistente.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import stat
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1]
DEFAULT_TARGET = Path("/home/nicomonzi/X_56")
COMMON_SKIP_DIRS = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".matplotlib",
    ".venv_x56",
}
COMMON_SKIP_SUFFIXES = {".pyc", ".pyo"}
RUNTIME_SUFFIXES = {
    ".nc",
    ".aer",
    ".mov",
    ".mod",
    ".jnt",
    ".usr",
    ".stdout",
    ".bylog",
}
TEXT_SUFFIXES = {
    ".py",
    ".sh",
    ".bat",
    ".md",
    ".txt",
    ".json",
    ".csv",
    ".xml",
    ".mbd",
    ".bdf",
    ".dat",
    ".in",
    ".c81",
    ".patch",
    ".env",
    ".example",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def ensure_empty_repository(target: Path) -> None:
    if not target.is_dir() or not (target / ".git").is_dir():
        raise SystemExit(f"La destinazione non e' una repository Git: {target}")
    unexpected = sorted(path.name for path in target.iterdir() if path.name != ".git")
    if unexpected:
        raise SystemExit(
            "La destinazione non e' vuota; nessun file e' stato modificato: "
            + ", ".join(unexpected)
        )


def should_skip(
    relative: Path,
    *,
    skip_dirs: set[str],
    skip_suffixes: set[str],
    skip_names: set[str],
) -> bool:
    if any(part in COMMON_SKIP_DIRS or part in skip_dirs for part in relative.parts):
        return True
    if relative.name in skip_names:
        return True
    if relative.suffix.lower() in COMMON_SKIP_SUFFIXES | skip_suffixes:
        return True
    return False


def copy_filtered(
    source: Path,
    target: Path,
    *,
    skip_dirs: set[str] | None = None,
    skip_suffixes: set[str] | None = None,
    skip_names: set[str] | None = None,
    fem_assets: dict[str, Path] | None = None,
) -> None:
    skip_dirs = skip_dirs or set()
    skip_suffixes = skip_suffixes or set()
    skip_names = skip_names or set()
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        if should_skip(
            relative,
            skip_dirs=skip_dirs,
            skip_suffixes=skip_suffixes,
            skip_names=skip_names,
        ):
            continue
        destination = target / relative
        if path.is_symlink():
            continue
        if path.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
            continue
        if not path.is_file():
            continue
        if path.suffix.lower() == ".fem" and fem_assets:
            digest = sha256(path)
            asset = fem_assets.get(digest)
            if asset:
                destination.parent.mkdir(parents=True, exist_ok=True)
                relative_asset = os.path.relpath(asset, destination.parent)
                destination.symlink_to(relative_asset)
                continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)


def copy_files(source: Path, target: Path, names: list[str]) -> None:
    target.mkdir(parents=True, exist_ok=True)
    for name in names:
        path = source / name
        if path.is_file():
            shutil.copy2(path, target / path.name)


def symlink_bulk(link: Path, canonical: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists() or link.is_symlink():
        raise RuntimeError(f"destinazione BULK gia' presente: {link}")
    link.symlink_to(os.path.relpath(canonical, link.parent), target_is_directory=True)


def replace(path: Path, old: str, new: str, required: bool = True) -> None:
    text = path.read_text(encoding="utf-8", errors="replace")
    count = text.count(old)
    if required and count == 0:
        raise RuntimeError(f"pattern non trovato in {path}: {old!r}")
    if count:
        path.write_text(text.replace(old, new), encoding="utf-8")


def rewrite_old_absolute_paths(target: Path) -> None:
    mapping = {
        "/home/nicomonzi/TESI/BFF_open_loop": str(target / "workflows/bff_open_loop"),
        "/home/nicomonzi/TESI/BFF_maneuver_envelope": str(target / "workflows/maneuver_bff"),
        "/home/nicomonzi/TESI/MANOUVER_STIFNESS": str(target / "workflows/maneuver_bff"),
        "/home/nicomonzi/TESI/BFF_DUST_55": str(target / "workflows/coupled_dust"),
        "/home/nicomonzi/TESI/_trim_coupled_stage": str(target / "workflows/coupled_dust/trim_55_stage"),
        "/home/nicomonzi/TESI/NASTRAN/REALASED_MODEL": str(target / "models/nastran/x56_r11"),
        "/home/nicomonzi/TESI/NASTRAN_SIMULATIONS/01_SOL103_60_MODES": str(target / "validation/modal/sol103_60_modes"),
        "/home/nicomonzi/TESI/NASTRAN_SIMULATIONS/02_COUPMASS_STUDY": str(target / "validation/modal/coupmass_study"),
        "/home/nicomonzi/TESI/NASTRAN_SIMULATIONS/03_GRAVITY_5G": str(target / "validation/static_5g"),
        "/home/nicomonzi/TESI/X56_AERO_POLAR": str(target / "validation/aero_polar"),
        "/home/nicomonzi/TESI/TRIM": str(target / "workflows/trim"),
        "/home/nicomonzi/TESI/DUST": str(target / "models/dust/x56"),
        "/home/nicomonzi/TESI/MBDYN": str(target / "archive/legacy_mbdyn"),
        "/home/nicomonzi/TESI/TEST": str(target / "archive/prototypes_test"),
        "/home/nicomonzi/TESI/bbf_manouver": str(target / "archive/superseded_campaigns/bbf_manouver"),
        "/home/nicomonzi/TESI/bff_eigenalaisy": str(target / "archive/superseded_campaigns/bff_eigenalaisy"),
        "/home/nicomonzi/TESI/bff_longitudinal": str(target / "archive/superseded_campaigns/bff_longitudinal"),
    }
    roots = [target / name for name in ("workflows", "models", "validation", "tools")]
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.is_symlink() or not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            updated = text
            for old, new in mapping.items():
                updated = updated.replace(old, new)
            if updated != text:
                path.write_text(updated, encoding="utf-8")


def make_portable(target: Path) -> None:
    bff = target / "workflows/bff_open_loop"
    maneuver = target / "workflows/maneuver_bff"
    trim = target / "workflows/trim"
    aero = target / "validation/aero_polar"

    replace(
        bff / "INCLUDE/modaljoint.mbd",
        '"mbdyn_modal.fem"',
        '"../../../assets/fem/mbdyn_modal_60.fem"',
    )
    replace(
        trim / "INCLUDE/modaljoint.mbd",
        '"../../MBDYN/BBF/INCLUDE/mbdyn_modal.fem"',
        '"../../../assets/fem/mbdyn_modal_40.fem"',
    )
    replace(
        trim / "INCLUDE/aerobody.mbd",
        '"../../MBDYN/BBF/INCLUDE/naca0012.c81"',
        '"../../../models/aerodynamics/c81/baseline/naca0012.c81"',
    )

    config_path = maneuver / "study_config.json"
    config = config_path.read_text(encoding="utf-8")
    absolute_baseline = str(target / "workflows/bff_open_loop")
    config_path.write_text(
        config.replace(f'"baseline_directory": "{absolute_baseline}"', '"baseline_directory": "../bff_open_loop"'),
        encoding="utf-8",
    )
    replace(
        maneuver / "maneuver_case.py",
        'BASELINE = Path(CONFIG["baseline_directory"])',
        'BASELINE_VALUE = Path(CONFIG["baseline_directory"])\nBASELINE = BASELINE_VALUE if BASELINE_VALUE.is_absolute() else (ROOT / BASELINE_VALUE).resolve()',
    )
    replace(
        maneuver / "analyse_time_domain_pairs.py",
        'BASELINE = Path(CONFIG["baseline_directory"])',
        'BASELINE_VALUE = Path(CONFIG["baseline_directory"])\nBASELINE = BASELINE_VALUE if BASELINE_VALUE.is_absolute() else (ROOT / BASELINE_VALUE).resolve()',
    )

    campaign_config_path = maneuver / "campaign_config.json"
    campaign_config = campaign_config_path.read_text(encoding="utf-8")
    replacements = {
        str(target / "workflows/maneuver_bff"): ".",
        str(target / "workflows/bff_open_loop/INCLUDE/mbdyn_modal.fem"): "../bff_open_loop/INCLUDE/mbdyn_modal.fem",
        str(target / "workflows/bff_open_loop/INCLUDE/modaljoint.mbd"): "../bff_open_loop/INCLUDE/modaljoint.mbd",
        str(target / "workflows/bff_open_loop/main_bff_open_loop.mbd"): "../bff_open_loop/main_bff_open_loop.mbd",
        str(target / "workflows/bff_open_loop"): "../bff_open_loop",
    }
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        campaign_config = campaign_config.replace(old, new)
    campaign_config_path.write_text(campaign_config, encoding="utf-8")
    replace(
        maneuver / "campaign.py",
        'MANEUVER_DIR = Path(MODEL["maneuver_directory"])\nBASELINE_DIR = Path(MODEL["baseline_directory"])',
        'def _model_path(value: str) -> Path:\n'
        '    path = Path(value)\n'
        '    return path if path.is_absolute() else (ROOT / path).resolve()\n\n\n'
        'MANEUVER_DIR = _model_path(MODEL["maneuver_directory"])\n'
        'BASELINE_DIR = _model_path(MODEL["baseline_directory"])',
    )
    replace(
        maneuver / "run_pullup_causal_controls.py",
        'BASE = ROOT.parent / "BFF_open_loop"',
        'BASE = ROOT.parent / "bff_open_loop"',
    )
    replace(
        maneuver / "analyse_prestressed_modes.py",
        'BASELINE_F06 = ROOT.parent / "NASTRAN_SIMULATIONS/01_SOL103_60_MODES/MAIN/sol103_60_modes.f06"\n'
        'BASELINE_FEM = ROOT.parent / "BFF_open_loop/INCLUDE/mbdyn_modal.fem"',
        'REPO_ROOT = ROOT.parents[1]\n'
        'BASELINE_F06 = REPO_ROOT / "validation/modal/sol103_60_modes/MAIN/sol103_60_modes.f06"\n'
        'BASELINE_FEM = REPO_ROOT / "assets/fem/mbdyn_modal_60.fem"',
    )

    replace(
        bff / "aero_static_validation.py",
        'ROOT = Path(__file__).resolve().parent\nMODEL = ROOT / "main_bff_open_loop.mbd"',
        'ROOT = Path(__file__).resolve().parent\nREPO_ROOT = ROOT.parents[1]\nMODEL = ROOT / "main_bff_open_loop.mbd"',
    )
    replace(
        bff / "aero_static_validation.py",
        f'NASTRAN_JSON = Path("{target}/validation/aero_polar/test/nastran_coefficients.json")',
        'NASTRAN_JSON = REPO_ROOT / "validation/aero_polar/test/nastran_coefficients.json"',
    )
    replace(
        bff / "dust_static_validation.py",
        'ROOT = Path(__file__).resolve().parent\nSOURCE = Path("' + str(target / "models/dust/x56") + '")',
        'ROOT = Path(__file__).resolve().parent\nREPO_ROOT = ROOT.parents[1]\nSOURCE = REPO_ROOT / "models/dust/x56"',
    )
    replace(
        bff / "nastran_flutter_reference.py",
        'RELEASED_F06 = Path(\n    "' + str(target / "models/nastran/x56_r11/FLUTTER_TEST") + '/"\n    "nsfluttr_cfg24611_10lb.f06"\n)',
        'REPO_ROOT = Path(__file__).resolve().parents[2]\nRELEASED_F06 = REPO_ROOT / "models/nastran/x56_r11/FLUTTER_TEST/nsfluttr_cfg24611_10lb.f06"',
    )
    replace(
        bff / "run_sweep.sh",
        "cd " + str(target / "workflows/bff_open_loop"),
        'cd "$(dirname "$0")"',
    )
    replace(
        target / "models/nastran/x56_r11/MAIN/export_x56_paraview.py",
        'DEFAULT_OUTPUT = Path("/home/nicomonzi/TESI/X56_PARAVIEW")',
        'DEFAULT_OUTPUT = Path(__file__).resolve().parents[3] / "results/figures/paraview"',
    )

    dynamic_replacements = {
        aero / "test/run_mbdyn_simulation.py": [
            (f'test_dir = Path("{target}/validation/aero_polar/test")', 'test_dir = Path(__file__).resolve().parent'),
            (f'mbdyn_dir = Path("{target}/validation/aero_polar/mbdyn")', 'mbdyn_dir = test_dir.parent / "mbdyn"'),
        ],
        aero / "test/compare_results.py": [
            (f'test_dir = Path("{target}/validation/aero_polar/test")', 'test_dir = Path(__file__).resolve().parent'),
        ],
        aero / "test/create_c81_file.py": [
            (f'test_dir = Path("{target}/validation/aero_polar/test")', 'test_dir = Path(__file__).resolve().parent'),
            (f'mbdyn_include = Path("{target}/validation/aero_polar/mbdyn/INCLUDE")', 'mbdyn_include = test_dir.parent / "mbdyn/INCLUDE"'),
        ],
        aero / "test/create_spanwise_c81_files.py": [
            (f'test_dir = Path("{target}/validation/aero_polar/test")', 'test_dir = Path(__file__).resolve().parent'),
            (f'mbdyn_include = Path("{target}/validation/aero_polar/mbdyn/INCLUDE")', 'mbdyn_include = test_dir.parent / "mbdyn/INCLUDE"'),
        ],
        aero / "test/extract_nastran_aero.py": [
            (f'nastran_dir = Path("{target}/validation/aero_polar/nastran")', 'test_dir = Path(__file__).resolve().parent\n    nastran_dir = test_dir.parent / "nastran"'),
            (f'test_dir = Path("{target}/validation/aero_polar/test")\n    test_dir.mkdir(exist_ok=True)', 'test_dir.mkdir(exist_ok=True)'),
            (f'c81_file = Path("{target}/validation/aero_polar/mbdyn/INCLUDE/x56_nastran.c81")', 'c81_file = test_dir.parent / "mbdyn/INCLUDE/x56_nastran.c81"'),
        ],
        aero / "test/create_mbdyn_models.py": [
            (f'main_file = Path("{target}/validation/aero_polar/mbdyn/main_x56.mbd")', 'aero_root = Path(__file__).resolve().parent.parent\n    main_file = aero_root / "mbdyn/main_x56.mbd"'),
            (f'output_file = Path("{target}/validation/aero_polar/mbdyn/main_x56_comprehensive.mbd")', 'output_file = aero_root / "mbdyn/main_x56_comprehensive.mbd"'),
            (f'quick_test_file = Path("{target}/validation/aero_polar/mbdyn/main_x56_quicktest.mbd")', 'quick_test_file = aero_root / "mbdyn/main_x56_quicktest.mbd"'),
        ],
        aero / "test/create_aerobody_files.py": [
            (f'x56_aero_dir = Path("{target}/validation/aero_polar/mbdyn/INCLUDE")', 'x56_aero_dir = Path(__file__).resolve().parent.parent / "mbdyn/INCLUDE"'),
            (f'    bff_aero_dir = Path("{target}/archive/superseded_campaigns/bff_eigenalaisy/INCLUDE")\n', ''),
        ],
    }
    for path, pairs in dynamic_replacements.items():
        for old, new in pairs:
            replace(path, old, new)


def asset_manifest(target: Path) -> None:
    rows: list[list[object]] = []
    assets = target / "assets"
    for path in sorted(assets.rglob("*")):
        if path.is_file() and not path.is_symlink() and path.name != "SHA256SUMS.csv":
            rows.append([path.relative_to(target).as_posix(), path.stat().st_size, sha256(path)])
    stream = []
    output = target / "assets/SHA256SUMS.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "bytes", "sha256"])
        writer.writerows(rows)


def create_project_files(target: Path) -> None:
    write(
        target / ".gitignore",
        """# Python e ambienti locali
__pycache__/
*.py[cod]
.venv*/
.pytest_cache/
.mypy_cache/
.matplotlib/

# Build e runtime dei solver
build/
build-*/
build_*/
precice-run/
work/
**/output/
**/outputs/
**/mbdyn_output/
**/incomplete_attempt_*/

# Risultati binari grezzi: conservarli fuori Git con un manifest
*.nc
*.aer
*.mov
*.mod
*.jnt
*.usr
*.stdout
*.bylog

# File locali macchina
machine.env
.DS_Store
Thumbs.db
""",
    )
    write(
        target / "requirements.txt",
        """numpy
scipy
matplotlib
netCDF4
h5py
cycler
pyNastran
""",
    )
    write(
        target / "README.md",
        """# X-56 aeroelastic research workspace

Repository riorganizzata a partire da `TESI` il 5 settembre 2026. La sorgente
originale non è stata modificata. Il rapporto storico completo è in
`docs/report/RAPPORTO_GENERALE.md`; le decisioni di migrazione sono in
`docs/decision_log/MIGRAZIONE_DA_TESI.md`.

## Workflow attivi

- `workflows/bff_open_loop`: baseline BFF steady con SAS-off;
- `workflows/maneuver_bff`: envelope, campagne paired e studio prestress;
- `workflows/trim`: trim longitudinale MBDyn;
- `workflows/coupled_dust`: coupling MBDyn–DUST e staging trim 55 m/s.

I modelli e i dati condivisi sono sotto `models/` e `assets/`; le verifiche
scientifiche sono sotto `validation/`. `archive/` non è una sorgente attiva.

## Verifica dopo il clone

```bash
python3 tools/verify_dependencies.py
python3 workflows/bff_open_loop/run_case.py --help
python3 workflows/maneuver_bff/run_sweep.py --help
python3 workflows/trim/run_trim_sweep.py --help
python3 workflows/coupled_dust/run_case.py --help
```

Le simulazioni non partono con questi comandi. DUST/preCICE richiedono inoltre
binari e binding configurati in `workflows/coupled_dust/config/machine.env`.
Gli output pesanti devono restare fuori Git e vanno associati a manifest e hash.
""",
    )
    write(
        target / "docs/decision_log/MIGRAZIONE_DA_TESI.md",
        f"""# Migrazione da TESI

Data: 5 settembre 2026. Sorgente: `{SOURCE}`. Destinazione: `{target}`.

## Decisioni

- Il contenuto della sorgente è stato copiato, mai spostato.
- `BFF_maneuver_envelope` e `MANOUVER_STIFNESS` sono stati uniti in
  `workflows/maneuver_bff` perché il secondo importa direttamente il primo.
- I FEM identici sono conservati una sola volta in `assets/fem`; i vecchi nomi
  sono link simbolici relativi, compatibili con Linux/WSL.
- I BULK R11 canonici sono in `assets/nastran_bulk/x56_r11`; i casi Nastran
  attivi mantengono il nome `BULK` tramite link relativo.
- Cache, virtualenv, build, tentativi incompleti e output grezzi NC/AER/MOV/MOD/
  JNT/USR non sono stati copiati. Manifest, CSV, JSON, report, figure e input
  sorgente sono stati conservati.
- `archive` è genealogia scientifica: non è incluso nel gate di dipendenze dei
  workflow attivi.
- I percorsi esterni ZENO e Desktop restano configurabili perché contengono
  risultati non inclusi nel repository; non sono dipendenze sorgente interne.

## Mappa principale

| TESI | X_56 |
|---|---|
| `BFF_open_loop` | `workflows/bff_open_loop` |
| `BFF_maneuver_envelope` + `MANOUVER_STIFNESS` | `workflows/maneuver_bff` |
| `TRIM` | `workflows/trim` |
| `BFF_DUST_55` | `workflows/coupled_dust` |
| `DUST` | `models/dust/x56` |
| `NASTRAN/REALASED_MODEL` | `models/nastran/x56_r11` |
| `NASTRAN_SIMULATIONS` | `validation/modal` e `validation/static_5g` |
| `X56_AERO_POLAR` | `validation/aero_polar` |
| `MBDYN`, `TEST`, vecchi `bff_*` | `archive` |

## Limiti mantenuti espliciti

La produzione DUST resta bloccata; la validazione statica 5 g non soddisfa la
soglia esterna; il gate fisico del prestress è rettificato e i controlli causali
sono ancora necessari. La riorganizzazione non cambia questi verdetti.
""",
    )
    write(
        target / "docs/metodologia/README.md",
        """# Metodologia

Questa directory raccoglie documenti trasversali. I dettagli esecutivi restano
accanto ai workflow per evitare che documentazione e codice divergano.

- `FILTERS_AND_ACTUATORS.md`: filtri digitali, notch e attuatori MBDyn;
- `BFF_OPEN_LOOP_TECHNICAL_REPORT.md`: correzioni della baseline;
- `MANEUVER_MODAL_NASTRAN_DLM_REVIEW.md`: rettifica prestress/DLM;
- `NASTRAN_VALIDATION.md`: base modale, COUPMASS e convergenza 5 g.
""",
    )
    write(
        target / "models/mbdyn/x56_modal/README.md",
        """# Modello modale X-56

Le basi FEM binarie sono centralizzate in `assets/fem`. I file MBDyn che
definiscono nodi, joint e aerodinamica restano nei workflow perché cambiano con
il protocollo sperimentale. `mbdyn_modal_60.fem` è la base recente a 60 modi;
`mbdyn_modal_40.fem` conserva la base storica usata dal trim e dai prototipi.
""",
    )
    write(
        target / "assets/README.md",
        """# Asset condivisi

File grandi canonici, identificati anche in `SHA256SUMS.csv`. Non modificare un
asset in-place: aggiungere una nuova versione e aggiornare esplicitamente il
workflow consumatore. I link simbolici nel repository devono restare relativi.
""",
    )
    write(
        target / "results/README.md",
        """# Risultati

`manifests` conserva la provenienza delle campagne, `processed` i risultati
ridotti e `figures` le tavole. I risultati grezzi pesanti restano fuori Git e
devono essere identificati tramite percorso, commit e SHA-256 nel manifest.
""",
    )
    write(
        target / "archive/README.md",
        """# Archivio

Rami storici e prototipi conservati per provenienza. Non sono la base dei
workflow correnti e non rientrano nel gate automatico delle dipendenze. I file
runtime pesanti sono stati omessi; i FEM condivisi sono link agli asset.
""",
    )

    verifier = r'''#!/usr/bin/env python3
"""Verifica statica della repository riorganizzata senza lanciare solver."""
from __future__ import annotations
import ast
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ACTIVE = [ROOT / name for name in ("workflows", "models", "validation", "tools")]
TEXT = {".py", ".json", ".xml", ".mbd", ".bdf", ".dat", ".in", ".c81", ".sh", ".bat"}
OLD_TESIS_ROOT = "/" + "home/nicomonzi/TESI"
errors = []

for path in ROOT.rglob("*"):
    if path.is_symlink() and not path.exists():
        errors.append(f"link rotto: {path.relative_to(ROOT)} -> {os.readlink(path)}")

for base in ACTIVE:
    for path in base.rglob("*"):
        if path.is_symlink() or not path.is_file() or path.suffix.lower() not in TEXT:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if OLD_TESIS_ROOT in text:
            errors.append(f"vecchio percorso TESI: {path.relative_to(ROOT)}")
        if path.suffix == ".py":
            try:
                ast.parse(text, filename=str(path))
            except SyntaxError as error:
                errors.append(f"Python non compilabile: {path.relative_to(ROOT)}: {error}")
        elif path.suffix == ".json":
            try:
                json.loads(text)
            except json.JSONDecodeError as error:
                errors.append(f"JSON non valido: {path.relative_to(ROOT)}: {error}")

for path in (ROOT / "workflows").rglob("*.mbd"):
    text = path.read_text(encoding="utf-8", errors="replace")
    for value in re.findall(r'include:\s*"([^"]+)"', text, flags=re.IGNORECASE):
        if any(token in value for token in ("$", "{", "}")):
            continue
        candidate = Path(value)
        if candidate.is_absolute():
            candidates = [candidate]
        else:
            workflow_relative = path.relative_to(ROOT / "workflows")
            workflow_root = ROOT / "workflows" / workflow_relative.parts[0]
            candidates = [path.parent / candidate, workflow_root / candidate]
        if not any(item.exists() for item in candidates):
            errors.append(f"include MBDyn mancante: {path.relative_to(ROOT)} -> {value}")
    for value in re.findall(r'"([^"]+\.(?:fem|c81))"', text, flags=re.IGNORECASE):
        candidate = Path(value)
        candidate = candidate if candidate.is_absolute() else path.parent / candidate
        if not candidate.exists():
            errors.append(f"asset MBDyn mancante: {path.relative_to(ROOT)} -> {value}")

for base in (ROOT / "models/nastran", ROOT / "validation"):
    for path in base.rglob("*.bdf"):
        text = path.read_text(encoding="utf-8", errors="replace")
        for value in re.findall(r"(?im)^\s*INCLUDE\s+['\"]([^'\"]+)['\"]", text):
            candidate = path.parent / value
            if not candidate.exists():
                errors.append(f"include Nastran mancante: {path.relative_to(ROOT)} -> {value}")

for path in ROOT.rglob("*"):
    if path.is_file() and not path.is_symlink() and path.stat().st_size >= 100_000_000:
        errors.append(f"file >=100 MB incompatibile con GitHub: {path.relative_to(ROOT)}")

if errors:
    print("VERIFICA FALLITA")
    for error in errors:
        print("-", error)
    raise SystemExit(1)
print("VERIFICA OK: Python/JSON, link, path TESI, include MBDyn/Nastran e limite 100 MB")
'''
    write(target / "tools/verify_dependencies.py", verifier)
    (target / "tools/verify_dependencies.py").chmod(
        (target / "tools/verify_dependencies.py").stat().st_mode | stat.S_IXUSR
    )


def migrate(target: Path) -> None:
    ensure_empty_repository(target)

    # Directory richieste, incluse quelle inizialmente solo documentali.
    for relative in (
        "docs/report", "docs/metodologia", "docs/decision_log",
        "models/nastran/x56_r11", "models/mbdyn/x56_modal",
        "models/aerodynamics/c81", "models/dust/x56",
        "workflows/trim", "workflows/bff_open_loop", "workflows/maneuver_bff",
        "workflows/coupled_dust", "validation/modal", "validation/static_5g",
        "validation/aero_polar", "validation/coupling_mesh",
        "results/manifests", "results/processed", "results/figures",
        "assets/fem", "assets/nastran_bulk/x56_r11", "tools",
        "archive/legacy_mbdyn", "archive/prototypes_test",
        "archive/superseded_campaigns",
    ):
        (target / relative).mkdir(parents=True, exist_ok=True)

    # Asset canonici prima di ogni altro albero: le copie FEM diventano link.
    fem_sources = {
        "mbdyn_modal_40.fem": SOURCE / "MBDYN/BBF/INCLUDE/mbdyn_modal.fem",
        "mbdyn_modal_60.fem": SOURCE / "BFF_open_loop/INCLUDE/mbdyn_modal.fem",
        "nsvibe_test.fem": SOURCE / "NASTRAN/FEMParser_OK/nsvibe_test.fem",
    }
    fem_assets: dict[str, Path] = {}
    for name, source in fem_sources.items():
        destination = target / "assets/fem" / name
        shutil.copy2(source, destination)
        fem_assets[sha256(source)] = destination

    bulk = target / "assets/nastran_bulk/x56_r11"
    copy_filtered(SOURCE / "NASTRAN/REALASED_MODEL/BULK", bulk)

    # Documentazione e modelli canonici.
    copy_filtered(SOURCE / "report", target / "docs/report")
    copy_filtered(
        SOURCE / "NASTRAN/REALASED_MODEL",
        target / "models/nastran/x56_r11",
        skip_dirs={"BULK", ".venv_x56", ".matplotlib"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "models/nastran/x56_r11/BULK", bulk)
    copy_filtered(SOURCE / "DUST", target / "models/dust/x56")
    copy_filtered(SOURCE / "BFF_open_loop/INCLUDE", target / "models/mbdyn/x56_modal/includes", fem_assets=fem_assets)
    for source_dir, name in (
        (SOURCE / "BFF_open_loop/INCLUDE", "baseline"),
        (SOURCE / "X56_AERO_POLAR/mbdyn/INCLUDE", "polar_study"),
    ):
        copy_filtered(
            source_dir,
            target / "models/aerodynamics/c81" / name,
            skip_names={path.name for path in source_dir.iterdir() if path.is_file() and path.suffix.lower() != ".c81"},
            fem_assets=fem_assets,
        )

    # Workflow attivi.
    runtime = set(RUNTIME_SUFFIXES)
    copy_filtered(
        SOURCE / "BFF_open_loop", target / "workflows/bff_open_loop",
        skip_suffixes=runtime, fem_assets=fem_assets,
    )
    copy_filtered(
        SOURCE / "TRIM", target / "workflows/trim",
        skip_suffixes=runtime | {".log", ".out"}, fem_assets=fem_assets,
    )
    copy_filtered(
        SOURCE / "BFF_maneuver_envelope", target / "workflows/maneuver_bff",
        skip_dirs={"output_n1_verification"}, skip_suffixes=runtime,
        fem_assets=fem_assets,
    )
    envelope_readme = target / "workflows/maneuver_bff/README.md"
    if envelope_readme.exists():
        envelope_readme.rename(target / "workflows/maneuver_bff/README_ENVELOPE.md")
    copy_filtered(
        SOURCE / "MANOUVER_STIFNESS", target / "workflows/maneuver_bff",
        skip_dirs={"runs", "runs_superseded_full_matrix"},
        skip_suffixes=runtime,
        fem_assets=fem_assets,
    )
    copy_filtered(
        SOURCE / "BFF_DUST_55", target / "workflows/coupled_dust",
        skip_dirs={
            "output", "precice-run", "AILERON_ANIMATION", "DEFLECTION_CASES",
            "SURFACE_SEQUENCE_4231",
        },
        skip_suffixes=runtime | {".log", ".out", ".h5", ".vtu", ".pvd"},
        fem_assets=fem_assets,
    )
    copy_filtered(
        SOURCE / "_trim_coupled_stage", target / "workflows/coupled_dust/trim_55_stage",
        skip_suffixes=runtime | {".log", ".out", ".h5", ".vtu", ".pvd"},
        fem_assets=fem_assets,
    )

    # Validazioni.
    copy_filtered(
        SOURCE / "NASTRAN_SIMULATIONS/01_SOL103_60_MODES",
        target / "validation/modal/sol103_60_modes",
        skip_dirs={"BULK"}, skip_suffixes=runtime | {".log", ".out"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "validation/modal/sol103_60_modes/BULK", bulk)
    copy_filtered(
        SOURCE / "NASTRAN_SIMULATIONS/02_COUPMASS_STUDY",
        target / "validation/modal/coupmass_study",
        skip_dirs={"BULK"}, skip_suffixes=runtime | {".log", ".out"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "validation/modal/coupmass_study/BULK", bulk)
    copy_filtered(
        SOURCE / "GRAFICI/MODAL_COMPARISON_10LB",
        target / "validation/modal/mode_shape_comparison",
        skip_dirs={"BULK"}, skip_suffixes=runtime | {".log", ".out"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "validation/modal/mode_shape_comparison/nastran/BULK", bulk)
    copy_filtered(
        SOURCE / "NASTRAN_SIMULATIONS/03_GRAVITY_5G",
        target / "validation/static_5g",
        skip_dirs={"BULK"}, skip_suffixes=runtime | {".log", ".out"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "validation/static_5g/nastran/BULK", bulk)
    copy_filtered(
        SOURCE / "GRAFICI/SOL101",
        target / "validation/static_5g/shape_comparison",
        skip_dirs={"BULK"}, skip_suffixes=runtime | {".log", ".out"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "validation/static_5g/shape_comparison/nastran/BULK", bulk)
    copy_filtered(
        SOURCE / "X56_AERO_POLAR", target / "validation/aero_polar",
        skip_dirs={"BULK", "mbdyn_output"}, skip_suffixes=runtime | {".log", ".out"},
        fem_assets=fem_assets,
    )
    symlink_bulk(target / "validation/aero_polar/nastran/BULK", bulk)
    copy_filtered(SOURCE / "BFF_DUST_55/meshes", target / "validation/coupling_mesh/meshes")
    copy_filtered(SOURCE / "BFF_DUST_55/reports", target / "validation/coupling_mesh/reports")
    copy_files(
        SOURCE / "BFF_DUST_55/tools", target / "validation/coupling_mesh/tools",
        ["audit_mesh.py", "configure_meshes.py"],
    )

    # Risultati elaborati e figure.
    copy_filtered(SOURCE / "results", target / "results/processed/global_comparisons")
    copy_filtered(SOURCE / "MANOUVER_STIFNESS/results", target / "results/processed/maneuver_bff")
    copy_filtered(SOURCE / "GRAFICI/X56_BEAM_GEOMETRY_PLOTS", target / "results/figures/beam_geometry")
    copy_filtered(SOURCE / "GRAFICI/X56_PARAVIEW", target / "results/figures/paraview")
    for path in sorted((SOURCE / "MANOUVER_STIFNESS").rglob("manifest.*")):
        relative = path.relative_to(SOURCE / "MANOUVER_STIFNESS")
        destination = target / "results/manifests/maneuver_bff" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)

    # Archivio leggero: sorgenti, input e risultati elaborati, non runtime grezzo.
    copy_filtered(
        SOURCE / "MBDYN", target / "archive/legacy_mbdyn",
        skip_dirs={"output"}, skip_suffixes=runtime | {".log", ".out", ".m"},
        fem_assets=fem_assets,
    )
    copy_filtered(
        SOURCE / "TEST", target / "archive/prototypes_test",
        skip_dirs={"output", "trim_sweep_results", "build_dust"},
        skip_suffixes=runtime | {".log", ".out"}, fem_assets=fem_assets,
    )
    for name in ("bbf_manouver", "bff_eigenalaisy", "bff_longitudinal"):
        copy_filtered(
            SOURCE / name, target / "archive/superseded_campaigns" / name,
            skip_dirs={"output"}, skip_suffixes=runtime | {".log", ".out"},
            fem_assets=fem_assets,
        )
    copy_filtered(
        SOURCE / "MANOUVER_STIFNESS/runs_superseded_full_matrix",
        target / "archive/superseded_campaigns/maneuver_full_matrix",
        skip_suffixes={".mbd"},
    )

    create_project_files(target)
    copy_files(
        SOURCE / "TEST", target / "docs/metodologia", ["FILTERS_AND_ACTUATORS.md"]
    )
    copy_files(
        SOURCE / "BFF_open_loop", target / "docs/metodologia", ["TECHNICAL_REPORT.md"]
    )
    source = target / "docs/metodologia/TECHNICAL_REPORT.md"
    if source.exists():
        source.rename(target / "docs/metodologia/BFF_OPEN_LOOP_TECHNICAL_REPORT.md")
    shutil.copy2(
        SOURCE / "MANOUVER_STIFNESS/REVIEW_MODAL_NASTRAN_DLM.md",
        target / "docs/metodologia/MANEUVER_MODAL_NASTRAN_DLM_REVIEW.md",
    )
    shutil.copy2(
        SOURCE / "NASTRAN_SIMULATIONS/README.md",
        target / "docs/metodologia/NASTRAN_VALIDATION.md",
    )

    rewrite_old_absolute_paths(target)
    make_portable(target)
    asset_manifest(target)

    print(f"Migrazione completata in {target}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", nargs="?", type=Path, default=DEFAULT_TARGET)
    args = parser.parse_args()
    migrate(args.target.expanduser().resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
