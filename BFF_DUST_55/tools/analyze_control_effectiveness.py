#!/usr/bin/env python3
"""Compute central-difference DUST effectiveness after five rigid cases exist."""
from __future__ import annotations
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/control_effectiveness/loads.csv"


def main() -> None:
    if not INPUT.exists():
        raise SystemExit(
            "Rigid DUST cases have not been run. Expected columns "
            "case,Fz_lbf,My_lbfin in results/control_effectiveness/loads.csv"
        )
    with INPUT.open(newline="") as stream:
        values = {row["case"]: row for row in csv.DictReader(stream)}
    required = {"theta_minus", "theta_plus", "bf_minus", "bf_plus"}
    if not required <= values.keys():
        raise SystemExit("Missing cases: " + ", ".join(sorted(required-values.keys())))
    def derivative(plus: str, minus: str, field: str) -> float:
        return (float(values[plus][field])-float(values[minus][field]))/2.0
    result = {
        "dFz_dtheta_lbf_per_deg": derivative("theta_plus", "theta_minus", "Fz_lbf"),
        "dMy_dtheta_lbfin_per_deg": derivative("theta_plus", "theta_minus", "My_lbfin"),
        "dFz_ddelta_BF_lbf_per_deg": derivative("bf_plus", "bf_minus", "Fz_lbf"),
        "dMy_ddelta_BF_lbfin_per_deg": derivative("bf_plus", "bf_minus", "My_lbfin"),
    }
    output = INPUT.with_name("effectiveness.json")
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__": main()
