#!/usr/bin/env python3
"""Rebuild the three committed discretizations from the common topology."""
from __future__ import annotations

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "meshes/COARSE/parametric_mesh.in"
LEVELS = {
    "COARSE": (5, [2, 2, 3, 3, 5, 30], 900),
    "MEDIUM": (10, [2, 3, 4, 5, 9, 34], 2280),
    "FINE": (17, [2, 3, 4, 5, 9, 40], 4284),
}


def main() -> None:
    text = BASE.read_text()
    text = re.sub(r"^coupling_node_file\s*=.*$",
                  "coupling_node_file = model/dust/coupling_nodes.in",
                  text, flags=re.MULTILINE)
    text = re.sub(r"^airfoil\s*=.*airfoilsection/([^/\s]+)$",
                  r"airfoil = model/dust/airfoilsection/\1",
                  text, flags=re.MULTILINE)
    for level, (n_chord, spans, expected) in LEVELS.items():
        configured = re.sub(r"^nelem_chord\s*=\s*\d+$",
                            f"nelem_chord = {n_chord}", text,
                            flags=re.MULTILINE)
        iterator = iter(spans)
        configured, count = re.subn(
            r"^nelem_span\s*=\s*\d+$",
            lambda _: f"nelem_span = {next(iterator)}",
            configured, flags=re.MULTILINE,
        )
        if count != 6:
            raise RuntimeError(f"Expected six span regions, found {count}")
        panels = 4 * n_chord * sum(spans)
        if panels != expected:
            raise RuntimeError(f"{level}: computed {panels}, expected {expected}")
        (ROOT / f"meshes/{level}/parametric_mesh.in").write_text(configured)
        print(f"{level}: chord={n_chord}, spans={spans}, panels={panels}")


if __name__ == "__main__":
    main()
