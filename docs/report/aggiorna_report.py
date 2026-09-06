#!/usr/bin/env python3
"""Aggiorna inventario, duplicati e registro progressi della cartella TESI.

La scansione e' in sola lettura fuori da ``report``. I prodotti di questo
script sono volutamente separati dai rapporti interpretativi scritti a mano.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import os
import subprocess
import time
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


REPORT_DIR = Path(__file__).resolve().parent
ROOT = REPORT_DIR.parent
DATA_DIR = REPORT_DIR / "dati"
SNAPSHOT_PATH = DATA_DIR / "snapshot.json"
PROGRESS_PATH = REPORT_DIR / "PROGRESSI.md"
EXCLUDED_TOP_LEVEL = {".git", "report"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def run_git(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, text=True, capture_output=True, check=False
    )
    return result.stdout.strip()


def regular_files() -> tuple[list[Path], int]:
    files: list[Path] = []
    symlinks = 0
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if relative.parts and relative.parts[0] in EXCLUDED_TOP_LEVEL:
            continue
        if path.is_symlink():
            symlinks += 1
        elif path.is_file():
            files.append(path)
    return sorted(files), symlinks


def project_directories() -> list[Path]:
    directories: list[Path] = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if relative.parts and relative.parts[0] in EXCLUDED_TOP_LEVEL:
            continue
        if not path.is_symlink() and path.is_dir():
            directories.append(path)
    return sorted(directories)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_text(header: list[str], rows: list[list[object]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue()


def git_history_rows() -> list[list[object]]:
    raw = run_git("log", "--date=iso-strict", "--pretty=format:%H%x1f%ad%x1f%an%x1f%s")
    rows: list[list[object]] = []
    for index, line in enumerate(raw.splitlines(), start=1):
        fields = line.split("\x1f", 3)
        if len(fields) == 4:
            rows.append([index, *fields])
    return rows


def load_previous() -> dict[str, object] | None:
    if not SNAPSHOT_PATH.exists():
        return None
    try:
        value = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def format_bytes(value: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    number = float(value)
    for unit in units:
        if abs(number) < 1024.0 or unit == units[-1]:
            return f"{number:.2f} {unit}"
        number /= 1024.0
    raise AssertionError("unreachable")


def append_progress(
    snapshot: dict[str, object], previous: dict[str, object] | None, note: str
) -> None:
    if not PROGRESS_PATH.exists():
        atomic_write(
            PROGRESS_PATH,
            "# Registro progressi di TESI\n\n"
            "Registro append-only generato da `python3 report/aggiorna_report.py`. "
            "La cartella `report` e `.git` sono escluse dalla scansione.\n\n",
        )

    current_files = snapshot["files"]
    assert isinstance(current_files, dict)
    previous_files = previous.get("files", {}) if previous else {}
    if not isinstance(previous_files, dict):
        previous_files = {}

    current_names = set(current_files)
    previous_names = set(previous_files)
    added = sorted(current_names - previous_names)
    removed = sorted(previous_names - current_names)
    modified = sorted(
        name
        for name in current_names & previous_names
        if current_files[name].get("sha256") != previous_files[name].get("sha256")
    )

    lines = [
        f"## Snapshot {snapshot['generated_at_utc']}",
        "",
        f"- Git HEAD: `{snapshot['git_head']}`; stato: "
        f"{snapshot['git_status_entries']} voci non pulite.",
        f"- File regolari: {snapshot['regular_file_count']}; dimensione: "
        f"{format_bytes(int(snapshot['regular_file_bytes']))}; symlink: "
        f"{snapshot['symlink_count']}.",
        f"- Duplicati: {snapshot['duplicate_groups']} gruppi, "
        f"{snapshot['redundant_extra_copies']} copie oltre la prima, massimo "
        f"teorico recuperabile {format_bytes(int(snapshot['potential_reclaim_bytes']))}.",
    ]
    if previous:
        lines.append(
            f"- Variazioni dal precedente snapshot: +{len(added)} aggiunti, "
            f"-{len(removed)} rimossi, {len(modified)} modificati."
        )
    else:
        lines.append("- Snapshot iniziale: non esisteva un confronto precedente.")
    if note:
        lines.append(f"- Nota: {note}")

    def add_examples(label: str, names: list[str]) -> None:
        if names:
            examples = ", ".join(f"`{name}`" for name in names[:12])
            suffix = f" (+{len(names) - 12} altri)" if len(names) > 12 else ""
            lines.append(f"- {label}: {examples}{suffix}.")

    if previous:
        add_examples("Aggiunti", added)
        add_examples("Rimossi", removed)
        add_examples("Modificati", modified)
    lines.extend(["", ""])
    with PROGRESS_PATH.open("a", encoding="utf-8") as stream:
        stream.write("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", default="", help="nota da aggiungere al registro")
    parser.add_argument(
        "--force-progress",
        action="store_true",
        help="aggiunge una voce anche se il contenuto del progetto non e' cambiato",
    )
    args = parser.parse_args()

    started = time.monotonic()
    files, symlink_count = regular_files()
    tracked = set(run_git("ls-files", "-z").split("\0"))
    tracked.discard("")
    status_lines = [
        line
        for line in run_git("status", "--short", "--untracked-files=all").splitlines()
        if not line[3:].startswith("report/")
    ]

    records: dict[str, dict[str, object]] = {}
    folder = defaultdict(lambda: {"files": 0, "bytes": 0, "ext": Counter()})
    hashes = defaultdict(list)
    zero_byte_files = 0

    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        stat = path.stat()
        digest = sha256(path)
        top = relative.split("/", 1)[0]
        extension = path.suffix.lower() or "[senza_estensione]"
        mtime = datetime.fromtimestamp(stat.st_mtime, timezone.utc).replace(
            microsecond=0
        ).isoformat()
        records[relative] = {
            "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": digest,
        }
        folder[top]["files"] += 1
        folder[top]["bytes"] += stat.st_size
        folder[top]["ext"][extension] += 1
        if stat.st_size == 0:
            zero_byte_files += 1
        else:
            hashes[(stat.st_size, digest)].append(relative)

    duplicate_groups = [
        (size, digest, sorted(paths))
        for (size, digest), paths in hashes.items()
        if len(paths) > 1
    ]
    duplicate_groups.sort(
        key=lambda item: (-(len(item[2]) - 1) * item[0], item[2][0])
    )
    duplicate_occurrences = sum(len(paths) for _, _, paths in duplicate_groups)
    redundant_copies = sum(len(paths) - 1 for _, _, paths in duplicate_groups)
    reclaim = sum((len(paths) - 1) * size for size, _, paths in duplicate_groups)

    inventory_rows: list[list[object]] = []
    for relative, data in records.items():
        path = Path(relative)
        inventory_rows.append(
            [
                relative,
                path.parts[0],
                path.suffix.lower() or "[senza_estensione]",
                data["bytes"],
                datetime.fromtimestamp(
                    int(data["mtime_ns"]) / 1_000_000_000, timezone.utc
                ).replace(microsecond=0).isoformat(),
                data["sha256"],
                "si" if relative in tracked else "no",
            ]
        )

    duplicate_rows: list[list[object]] = []
    for group_id, (size, digest, paths) in enumerate(duplicate_groups, start=1):
        for position, relative in enumerate(paths, start=1):
            duplicate_rows.append(
                [
                    group_id,
                    digest,
                    size,
                    len(paths),
                    (len(paths) - 1) * size,
                    position,
                    relative.split("/", 1)[0],
                    relative,
                ]
            )

    internal = defaultdict(lambda: {"groups": 0, "extra": 0, "reclaim": 0})
    for size, _digest, paths in duplicate_groups:
        by_top = Counter(path.split("/", 1)[0] for path in paths)
        for top, count in by_top.items():
            if count > 1:
                internal[top]["groups"] += 1
                internal[top]["extra"] += count - 1
                internal[top]["reclaim"] += (count - 1) * size

    folder_rows: list[list[object]] = []
    for top in sorted(folder):
        extensions = "; ".join(
            f"{key}:{value}" for key, value in folder[top]["ext"].most_common(12)
        )
        top_files = [name for name in records if name.split("/", 1)[0] == top]
        folder_rows.append(
            [
                top,
                folder[top]["files"],
                folder[top]["bytes"],
                sum(name in tracked for name in top_files),
                sum(name not in tracked for name in top_files),
                internal[top]["groups"],
                internal[top]["extra"],
                internal[top]["reclaim"],
                extensions,
            ]
        )

    directory_stats = defaultdict(
        lambda: {
            "direct_files": 0,
            "recursive_files": 0,
            "recursive_bytes": 0,
            "latest_mtime_utc": "",
            "ext": Counter(),
        }
    )
    for relative, data in records.items():
        file_path = Path(relative)
        ancestors = list(file_path.parents)[:-1]
        mtime = datetime.fromtimestamp(
            int(data["mtime_ns"]) / 1_000_000_000, timezone.utc
        ).replace(microsecond=0).isoformat()
        extension = file_path.suffix.lower() or "[senza_estensione]"
        for ancestor in ancestors:
            key = ancestor.as_posix()
            directory_stats[key]["recursive_files"] += 1
            directory_stats[key]["recursive_bytes"] += int(data["bytes"])
            directory_stats[key]["ext"][extension] += 1
            if mtime > directory_stats[key]["latest_mtime_utc"]:
                directory_stats[key]["latest_mtime_utc"] = mtime
        directory_stats[file_path.parent.as_posix()]["direct_files"] += 1

    directory_rows: list[list[object]] = []
    for directory in project_directories():
        relative = directory.relative_to(ROOT).as_posix()
        stats = directory_stats[relative]
        extensions = "; ".join(
            f"{key}:{value}" for key, value in stats["ext"].most_common(12)
        )
        directory_rows.append(
            [
                relative,
                len(Path(relative).parts),
                stats["direct_files"],
                stats["recursive_files"],
                stats["recursive_bytes"],
                stats["latest_mtime_utc"],
                extensions,
            ]
        )

    code_rows: list[list[object]] = []
    code_exclusions = {".venv_x56", "site-packages", "dist-packages", "build_dust", "__pycache__"}
    for relative in sorted(name for name in records if name.endswith(".py")):
        path = Path(relative)
        if any(part in code_exclusions for part in path.parts):
            continue
        source = (ROOT / path).read_text(encoding="utf-8", errors="replace")
        parse_status = "ok"
        classes: list[str] = []
        functions: list[str] = []
        imports: list[str] = []
        description = ""
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", SyntaxWarning)
                tree = ast.parse(source)
            description = ((ast.get_docstring(tree) or "").strip().splitlines() or [""])[0]
            for node in tree.body:
                if isinstance(node, ast.ClassDef):
                    classes.append(node.name)
                elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    functions.append(node.name)
                elif isinstance(node, ast.Import):
                    imports.extend(alias.name for alias in node.names)
                elif isinstance(node, ast.ImportFrom):
                    imports.append(node.module or ".")
        except (SyntaxError, ValueError) as error:
            parse_status = f"errore: {type(error).__name__}: {error}"
        code_rows.append(
            [
                relative,
                source.count("\n") + 1,
                parse_status,
                description,
                "; ".join(classes),
                "; ".join(functions),
                "; ".join(dict.fromkeys(imports)),
            ]
        )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    atomic_write(
        DATA_DIR / "inventario_file.csv",
        csv_text(
            ["percorso", "cartella", "estensione", "byte", "mtime_utc", "sha256", "tracciato_git"],
            inventory_rows,
        ),
    )
    atomic_write(
        DATA_DIR / "duplicati_sha256.csv",
        csv_text(
            [
                "gruppo",
                "sha256",
                "byte_file",
                "copie",
                "byte_ridondanti_massimi",
                "posizione_nel_gruppo",
                "cartella",
                "percorso",
            ],
            duplicate_rows,
        ),
    )
    atomic_write(
        DATA_DIR / "riepilogo_cartelle.csv",
        csv_text(
            [
                "cartella",
                "file_regolari",
                "byte",
                "file_tracciati_git",
                "file_non_tracciati_o_ignorati",
                "gruppi_duplicati_interni",
                "copie_extra_interne",
                "byte_recuperabili_interni_massimi",
                "estensioni_principali",
            ],
            folder_rows,
        ),
    )
    atomic_write(
        DATA_DIR / "inventario_sottocartelle.csv",
        csv_text(
            [
                "sottocartella",
                "profondita",
                "file_diretti",
                "file_ricorsivi",
                "byte_ricorsivi",
                "ultimo_mtime_utc",
                "estensioni_principali_ricorsive",
            ],
            directory_rows,
        ),
    )
    atomic_write(
        DATA_DIR / "catalogo_codice_python.csv",
        csv_text(
            [
                "percorso",
                "righe",
                "parsing_ast",
                "descrizione_docstring",
                "classi_top_level",
                "funzioni_top_level",
                "import_top_level",
            ],
            code_rows,
        ),
    )
    atomic_write(
        DATA_DIR / "cronologia_git.csv",
        csv_text(
            ["ordine_dal_piu_recente", "commit", "data", "autore", "messaggio"],
            git_history_rows(),
        ),
    )

    digest = hashlib.sha256()
    for relative, data in records.items():
        digest.update(relative.encode("utf-8", errors="surrogateescape"))
        digest.update(b"\0")
        digest.update(str(data["bytes"]).encode("ascii"))
        digest.update(b"\0")
        digest.update(str(data["sha256"]).encode("ascii"))
        digest.update(b"\n")

    previous = load_previous()
    snapshot: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "root": str(ROOT),
        "excluded_top_level": sorted(EXCLUDED_TOP_LEVEL),
        "regular_file_count": len(files),
        "regular_file_bytes": sum(int(data["bytes"]) for data in records.values()),
        "symlink_count": symlink_count,
        "zero_byte_files": zero_byte_files,
        "duplicate_groups": len(duplicate_groups),
        "duplicate_occurrences": duplicate_occurrences,
        "redundant_extra_copies": redundant_copies,
        "potential_reclaim_bytes": reclaim,
        "project_digest": digest.hexdigest(),
        "git_head": run_git("rev-parse", "--short=12", "HEAD") or "non disponibile",
        "git_commit_count": int(run_git("rev-list", "--all", "--count") or 0),
        "git_status_entries": len(status_lines),
        "duration_seconds": round(time.monotonic() - started, 3),
        "files": records,
    }
    changed = not previous or previous.get("project_digest") != snapshot["project_digest"]
    if changed or args.force_progress or args.note:
        append_progress(snapshot, previous, args.note)
    atomic_write(SNAPSHOT_PATH, json.dumps(snapshot, indent=2, ensure_ascii=False) + "\n")

    print(
        f"Scansione completata: {len(files)} file, {len(duplicate_groups)} gruppi "
        f"duplicati, {redundant_copies} copie extra."
    )
    print(f"Dati aggiornati in {DATA_DIR}")
    print("Registro progressi aggiornato." if changed or args.force_progress or args.note else "Nessuna variazione: registro invariato.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
