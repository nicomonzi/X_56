#!/usr/bin/env python3
"""Run a wake-developed rigid X-56 DUST polar on the available FINE mesh.

This is an independent sign/trend check.  The available topology has no
winglets and the repository mesh audit reports residual moment convergence
uncertainty, so DUST is not used as the sole calibration target.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[1]
SOURCE = REPO_ROOT / "models/dust/x56"
DEFAULT_DUST = Path("/home/nicomonzi/src/dust-master/build/bin/dust")
DEFAULT_PRE = Path("/home/nicomonzi/src/dust-master/build/bin/dust_pre")
DEFAULT_POST = Path("/home/nicomonzi/src/dust-master/build/bin/dust_post")
CG_LOCAL_IN = (63.187383385809, 0.110529571088, 1.079793658848)
S_REF_IN2 = 8064.0
C_REF_IN = 24.0
RHO = 9.7284e-8


def fine_rigid_mesh() -> str:
    text = (SOURCE / "parametric_mesh.in").read_text()
    # Use the vortex-lattice surface for a rigid polar.  The available closed
    # panel topology is not watertight (as documented by its audit) and drives
    # the pressure-panel solve to NaN from the first step.
    text = re.sub(r"(?m)^el_type\s*=\s*\w+", "el_type = v", text)
    text = re.sub(r"(?m)^nelem_chord\s*=\s*\d+", "nelem_chord = 17", text)
    divisions = iter((2, 3, 4, 5, 9, 40))
    text = re.sub(r"(?m)^nelem_span\s*=\s*\d+", lambda _m: f"nelem_span = {next(divisions)}", text)
    text = re.sub(r"(?m)^(amplitude\s*=\s*)[-+]?\d+(?:\.\d+)?", r"\g<1>0.0", text)
    if len(re.findall(r"(?m)^amplitude\s*=\s*0\.0", text)) != 10:
        raise RuntimeError("expected ten fixed zero-angle hinges")
    return text


def write_case(case: Path, angle_deg: float, velocity_mps: float, end_time: float) -> int:
    case.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE / "airfoilsection", case / "airfoilsection", dirs_exist_ok=True)
    (case / "parametric_mesh.in").write_text(fine_rigid_mesh())
    (case / "dust_pre.in").write_text(
        "comp_name = X_56\n"
        "geo_file = parametric_mesh.in\n"
        "ref_tag = centerbody\n"
        "file_name = geo_input.h5\n"
    )
    x, y, z = CG_LOCAL_IN
    (case / "References.in").write_text(
        "reference_tag = centerbody\nparent_tag = 0\n"
        "origin = (/0., 0., 0./)\n"
        "orientation = (/1.,0.,0., 0.,1.,0., 0.,0.,1./)\n"
        "multiple = F\nmoving = F\n\n"
        "reference_tag = CG\nparent_tag = centerbody\n"
        f"origin = (/{x:.12f}, {y:.12f}, {z:.12f}/)\n"
        "orientation = (/1.,0.,0., 0.,1.,0., 0.,0.,1./)\n"
        "multiple = F\nmoving = F\n"
    )
    speed = velocity_mps / 0.0254
    alpha = math.radians(angle_deg)
    frames = int(round(end_time / 0.05)) + 1
    (case / "dust.in").write_text(f"""basename = Output/case
geometry_file = geo_input.h5
reference_file = References.in

tstart = 0.0
tend = {end_time:.6f}
dt = 0.01
dt_out = 0.05
output_start = T
ndt_update_wake = 1

rho_inf = {RHO:.12e}
a_inf = 13385.8267716535
P_inf = 14.6959487755
u_inf = (/{speed * math.cos(alpha):.12f}, 0.0, {speed * math.sin(alpha):.12f}/)

particles_box_min = (/-500.0, -1000.0, -800.0/)
particles_box_max = (/3500.0, 1000.0, 800.0/)
fmm = T
box_length = 250.0
n_box = (/16, 8, 7/)
octree_origin = (/-500.0, -1000.0, -800.0/)
n_octree_levels = 4
min_octree_part = 10
multipole_degree = 2
vortstretch = T
diffusion = T
penetration_avoidance = T
n_wake_panels = 40
n_wake_particles = 500000
Kvortex_rad = 1.0
""")
    # Never average the t=0 initialization frame; it has no developed wake and
    # its pressure-load fields may legitimately be NaN.
    first_average = max(2, frames - 4)
    (case / "dust_post.in").write_text(f"""basename = Postpro/x56
data_basename = Output/case

analysis = {{
  type = integral_loads
  name = loads_cg
  start_res = {first_average}
  end_res = {frames}
  step_res = 1
  format = dat
  average = T
  component = X_56
  reference_tag = CG
}}
""")
    return frames


def run(command: list[str], cwd: Path, log_name: str, environment: dict[str, str]) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, env=environment)
    log = result.stdout + result.stderr
    (cwd / log_name).write_text(log)
    if result.returncode or "\nERROR in" in log:
        raise RuntimeError(f"{' '.join(command)} failed; see {cwd / log_name}")


def parse_average(case: Path, velocity_mps: float) -> dict[str, float]:
    files = sorted((case / "Postpro").glob("*loads_cg*.dat"))
    if not files:
        files = sorted((case / "Postpro").glob("*.dat"))
    if len(files) != 1:
        raise RuntimeError(f"expected one DUST load file, found {files}")
    rows = []
    for line in files[0].read_text().splitlines():
        if line.lstrip().startswith("#") or not line.strip():
            continue
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            continue
    if not rows:
        raise RuntimeError(f"no numeric loads in {files[0]}")
    values = rows[-1]
    if len(values) >= 18:
        values = values[:6]
    elif len(values) == 7:
        values = values[1:]
    if len(values) != 6:
        raise RuntimeError(f"unexpected integral-load columns in {files[0]}: {values}")
    fx, fy, fz, mx, my, mz = values
    speed = velocity_mps / 0.0254
    q_dyn = 0.5 * RHO * speed**2
    return {
        "CX": fx / (q_dyn * S_REF_IN2), "CY": fy / (q_dyn * S_REF_IN2), "CZ": fz / (q_dyn * S_REF_IN2),
        "CMX_CG": mx / (q_dyn * S_REF_IN2 * C_REF_IN),
        "CMY_CG": my / (q_dyn * S_REF_IN2 * C_REF_IN),
        "CMZ_CG": mz / (q_dyn * S_REF_IN2 * C_REF_IN),
        "load_file": str(files[0]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--angles", type=float, nargs="+", default=[-1.0, 0.0, 1.0])
    parser.add_argument("--velocity", type=float, default=60.8421)
    parser.add_argument("--end-time", type=float, default=0.80,
                        help="wake-development time; 0.8 s is about 5.9 spans")
    parser.add_argument("--output", type=Path, default=ROOT / "aero_validation/dust_fine")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dust", type=Path, default=DEFAULT_DUST)
    parser.add_argument("--dust-pre", type=Path, default=DEFAULT_PRE)
    parser.add_argument("--dust-post", type=Path, default=DEFAULT_POST)
    args = parser.parse_args()
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "8")
    results = {}
    for angle in args.angles:
        stem = f"alpha_{angle:+05.1f}".replace("+", "p").replace("-", "m").replace(".", "p")
        case = args.output / stem
        if args.overwrite and case.exists():
            shutil.rmtree(case)
        frames = write_case(case, angle, args.velocity, args.end_time)
        final_result = case / "Output" / f"case_res_{frames:04d}.h5"
        if not final_result.exists():
            (case / "Output").mkdir(exist_ok=True)
            run([str(args.dust_pre), "dust_pre.in"], case, "dust_pre.log", environment)
            run([str(args.dust), "dust.in"], case, "dust.log", environment)
        (case / "Postpro").mkdir(exist_ok=True)
        run([str(args.dust_post), "dust_post.in"], case, "dust_post.log", environment)
        results[angle] = parse_average(case, args.velocity)
        print(f"alpha={angle:+.1f}: CZ={results[angle]['CZ']:+.6f}, CMY_CG={results[angle]['CMY_CG']:+.6f}")
    angles = np.array(sorted(results))
    summary = {
        "mesh": "FINE vortex lattice, 17 chordwise divisions, 2142 panels nominal, no winglets",
        "wake_development_s": args.end_time,
        "wake_convection_spans": args.end_time * args.velocity / (336.0 * 0.0254),
        "limitations": ["available topology has no winglets", "closed pressure-panel mesh produced NaN and was replaced by the corresponding vortex lattice", "mesh audit reports 14.59% MEDIUM-to-FINE My difference"],
        "cases": results,
        "slopes_per_deg": {
            "CZ_alpha": float(np.polyfit(angles, [results[a]["CZ"] for a in angles], 1)[0]),
            "CMY_CG_alpha": float(np.polyfit(angles, [results[a]["CMY_CG"] for a in angles], 1)[0]),
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "dust_polar_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["slopes_per_deg"], indent=2))


if __name__ == "__main__":
    main()
