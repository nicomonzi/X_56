#!/usr/bin/env python3
"""Consolidate accepted MBDyn text states and preCICE diagnostics into CSV."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


SURFACES = {
    1004: "bfl", 1008: "wf1l", 1011: "wf2l", 1014: "wf3l", 1017: "wf4l",
    2004: "bfr", 2008: "wf1r", 2011: "wf2r", 2014: "wf3r", 2017: "wf4r",
}


def blocks(path: Path) -> list[dict[str, list[float]]]:
    """Split an MBDyn text output whenever its first record label repeats."""
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields:
            rows.append((fields[0], [float(value) for value in fields[1:]]))
    if not rows:
        raise RuntimeError(f"No records in {path}")
    first = rows[0][0]
    starts = [index for index, row in enumerate(rows) if row[0] == first]
    starts.append(len(rows))
    result = []
    for begin, end in zip(starts, starts[1:]):
        block = dict(rows[begin:end])
        if len(block) != end - begin:
            raise RuntimeError(f"Duplicate labels inside a state block in {path}")
        result.append(block)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="SMOKE or PRODUCTION output directory")
    parser.add_argument("--prefix", default="smoke", help="MBDyn output basename")
    args = parser.parse_args()
    output = args.output.resolve()

    response_path = output / "coupled_response.csv"
    with response_path.open(newline="") as stream:
        response = list(csv.DictReader(stream))
    mov = blocks(output / f"{args.prefix}.mov")
    mod = blocks(output / f"{args.prefix}.mod")
    jnt = blocks(output / f"{args.prefix}.jnt")
    structural_steps = len(mov) - 1
    if not (len(mov) == len(mod) == len(jnt)) or len(response) not in {
        structural_steps, structural_steps - 1
    }:
        raise RuntimeError(
            f"State count mismatch: response={len(response)}, mov={len(mov)}, "
            f"mod={len(mod)}, jnt={len(jnt)}; expected equal structural outputs "
            "and at most one missing terminal coupling diagnostic"
        )

    old_moment = bool(response) and "sum_nodal_my_lbfin" in response[0]
    dt = float(response[0]["time_s"]) if response else 0.0
    fieldnames = [
        "time_s", "state", "coupling_iterations", "theta_rad", "q_radps",
        "vz_inps", "z_in", "swb1_modal_q", "swb1_modal_qdot",
        "global_fz_lbf", "global_my_cg_lbfin", "moment_definition",
        "left_tip_dz_in", "right_tip_dz_in", "symmetric_tip_dz_in",
        "antisymmetric_tip_dz_in", "dust_particle_count",
        *(f"{name}_rad" for name in SURFACES.values()),
    ]
    destination = output / "diagnostics.csv"
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(1, structural_steps + 1):
            coupled = response[index - 1] if index <= len(response) else {}
            base = mov[index]["990000"]
            swb1 = mod[index]["5.7"]
            row = {
                "time_s": coupled.get("time_s", f"{index * dt:.9f}"),
                "state": coupled.get("state", response[-1]["state"] if response else ""),
                "coupling_iterations": coupled.get("coupling_iterations", ""),
                "theta_rad": f"{base[4]:.12e}", "q_radps": f"{base[10]:.12e}",
                "vz_inps": f"{base[8]:.12e}", "z_in": f"{base[2]:.12e}",
                "swb1_modal_q": f"{swb1[0]:.12e}",
                "swb1_modal_qdot": f"{swb1[1]:.12e}",
                "global_fz_lbf": coupled.get("total_fz_lbf", ""),
                "global_my_cg_lbfin": coupled.get("global_my_cg_lbfin",
                                                    coupled.get("sum_nodal_my_lbfin", "")),
                "moment_definition": ("sum_of_nodal_moments_only_legacy_smoke"
                                      if old_moment else "resultant_about_CG"),
                "left_tip_dz_in": coupled.get("left_tip_dz_in", ""),
                "right_tip_dz_in": coupled.get("right_tip_dz_in", ""),
                "symmetric_tip_dz_in": coupled.get("symmetric_tip_dz_in", ""),
                "antisymmetric_tip_dz_in": coupled.get("antisymmetric_tip_dz_in", ""),
                "dust_particle_count": "",
            }
            for label, name in SURFACES.items():
                # Total-joint output: relative rotation vector Phi occupies fields 16:18.
                row[f"{name}_rad"] = f"{jnt[index][str(label)][15]:.12e}"
            writer.writerow(row)
    print(f"Wrote {structural_steps} accepted states to {destination}")
    if len(response) != structural_steps:
        print("NOTE: terminal structural state has no aerodynamic diagnostic row; "
              "its load/tip fields are intentionally blank")


if __name__ == "__main__":
    main()
