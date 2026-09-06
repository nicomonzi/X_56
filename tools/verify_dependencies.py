#!/usr/bin/env python3
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
