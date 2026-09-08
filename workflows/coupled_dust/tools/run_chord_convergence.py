#!/usr/bin/env python3
"""DUST chordwise convergence study using the TRIM_DUST_1 X-56 geometry.

The span discretization is deliberately fixed at [2, 3, 4, 3, 9, 40].
Only ``nelem_chord`` changes.  The standalone study uses the vortex-lattice
surface because the closed pressure-panel topology is intended for the
coupled moving-grid run and is not robust as a rigid standalone polar.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
SOURCE = REPO / "workflows/TRIM_DUST_1/case"
CHORD_LEVELS = (10, 20, 30, 35, 40)
SPAN_LEVELS = (2, 3, 4, 3, 9, 40)
CG_LOCAL_IN = (63.187383385809, 0.110529571088, 1.079793658848)
S_REF_IN2 = 8064.0
C_REF_IN = 24.0
RHO = 9.7284e-8


def executable(value: str | None, fallback: str) -> Path:
    candidate = Path(value).expanduser() if value else Path(shutil.which(fallback) or "")
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise RuntimeError(f"eseguibile non disponibile: {fallback}")
    return candidate.resolve()


def render_mesh(n_chord: int) -> str:
    text = (SOURCE / "mesh_canonical.in").read_text()
    text = re.sub(r"(?m)^el_type\s*=\s*\w+", "el_type = v", text)
    text = re.sub(r"(?m)^nelem_chord\s*=\s*\d+", f"nelem_chord = {n_chord}", text)
    # The standalone 0.8.2 build uses the legacy spelling.  It implements the
    # same leading-edge cosine distribution called ``cosine_le`` by the
    # coupled development build.
    text = re.sub(r"(?m)^type_chord\s*=.*$", "type_chord = cosineLE", text)
    spans = iter(SPAN_LEVELS)
    text, count = re.subn(
        r"(?m)^nelem_span\s*=\s*\d+",
        lambda _match: f"nelem_span = {next(spans)}",
        text,
    )
    if count != len(SPAN_LEVELS):
        raise RuntimeError(f"attese sei regioni spanwise, trovate {count}")
    text, count = re.subn(
        r"(?m)^(amplitude\s*=\s*)[-+]?\d+(?:\.\d+)?",
        r"\g<1>0.0",
        text,
    )
    if count != 10:
        raise RuntimeError(f"attese dieci hinge, trovate {count}")
    return text


def write_case(case: Path, n_chord: int, angle_deg: float,
               velocity_mps: float, end_time: float) -> int:
    case.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE / "airfoilsection", case / "airfoilsection", dirs_exist_ok=True)
    (case / "mesh.in").write_text(render_mesh(n_chord))
    (case / "dust_pre.in").write_text(
        "comp_name = X56\ngeo_file = mesh.in\nref_tag = centerbody\n"
        "file_name = geo_input.h5\n"
    )
    x, y, z = CG_LOCAL_IN
    (case / "References.in").write_text(
        "reference_tag = centerbody\nparent_tag = 0\norigin = (/0.,0.,0./)\n"
        "orientation = (/1.,0.,0., 0.,1.,0., 0.,0.,1./)\nmultiple = F\nmoving = F\n\n"
        "reference_tag = CG\nparent_tag = centerbody\n"
        f"origin = (/{x:.12f},{y:.12f},{z:.12f}/)\n"
        "orientation = (/1.,0.,0., 0.,1.,0., 0.,0.,1./)\nmultiple = F\nmoving = F\n"
    )
    speed = velocity_mps / 0.0254
    alpha = math.radians(angle_deg)
    dt_out = 0.05
    frames = int(round(end_time / dt_out)) + 1
    (case / "dust.in").write_text(f"""basename = Output/case
geometry_file = geo_input.h5
reference_file = References.in
tstart = 0.0
tend = {end_time:.9g}
dt = 0.01
dt_out = {dt_out}
output_start = T
ndt_update_wake = 1
rho_inf = {RHO:.12e}
a_inf = 13385.8267716535
P_inf = 14.6959487755
u_inf = (/{speed*math.cos(alpha):.12f},0.0,{speed*math.sin(alpha):.12f}/)
particles_box_min = (/-500.0,-1000.0,-800.0/)
particles_box_max = (/3500.0,1000.0,800.0/)
fmm = T
box_length = 250.0
n_box = (/16,8,7/)
octree_origin = (/-500.0,-1000.0,-800.0/)
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
  component = X56
  reference_tag = CG
}}
""")
    return frames


def run(command: list[str], cwd: Path, log_name: str,
        environment: dict[str, str]) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True,
                            env=environment)
    log = result.stdout + result.stderr
    (cwd / log_name).write_text(log)
    if result.returncode or "\nERROR in" in log:
        raise RuntimeError(f"{' '.join(command)} fallito; vedere {cwd/log_name}")


def parse_loads(case: Path, velocity_mps: float) -> dict[str, float]:
    files = sorted((case / "Postpro").glob("*loads_cg*.dat"))
    if len(files) != 1:
        raise RuntimeError(f"file carichi inattesi: {files}")
    rows: list[list[float]] = []
    for line in files[0].read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            pass
    values = rows[-1]
    if len(values) >= 18:
        values = values[:6]
    elif len(values) == 7:
        values = values[1:]
    if len(values) != 6 or not np.isfinite(values).all():
        raise RuntimeError(f"carichi non validi: {values}")
    fx, fy, fz, mx, my, mz = values
    speed = velocity_mps / 0.0254
    q_dyn = 0.5 * RHO * speed**2
    return {
        "Fx_lbf": fx, "Fy_lbf": fy, "Fz_lbf": fz,
        "Mx_lbfin": mx, "My_lbfin": my, "Mz_lbfin": mz,
        "CX": fx/(q_dyn*S_REF_IN2), "CZ": fz/(q_dyn*S_REF_IN2),
        "CMY_CG": my/(q_dyn*S_REF_IN2*C_REF_IN),
    }


def relative_change(new: float, old: float) -> float:
    scale = max(abs(new), abs(old), 1.e-14)
    return 100.0 * abs(new - old) / scale


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--angle", type=float, default=1.0)
    parser.add_argument("--velocity", type=float, default=60.8421)
    parser.add_argument("--end-time", type=float, default=0.8)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/chord_convergence")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dust", default=os.environ.get("DUST_BIN"))
    parser.add_argument("--dust-pre", default=os.environ.get("DUST_PRE_BIN"))
    parser.add_argument("--dust-post", default=os.environ.get("DUST_POST_BIN"))
    args = parser.parse_args()
    bins = {
        "dust": executable(args.dust, "dust"),
        "pre": executable(args.dust_pre, "dust_pre"),
        "post": executable(args.dust_post, "dust_post"),
    }
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1")
    results: dict[str, dict[str, object]] = {}
    for n_chord in CHORD_LEVELS:
        case = args.output / f"chord_{n_chord:02d}"
        if args.overwrite and case.exists():
            shutil.rmtree(case)
        frames = write_case(case, n_chord, args.angle, args.velocity, args.end_time)
        final = case / "Output" / f"case_res_{frames:04d}.h5"
        if not final.exists():
            (case / "Output").mkdir(exist_ok=True)
            run([str(bins["pre"]), "dust_pre.in"], case, "dust_pre.log", env)
            run([str(bins["dust"]), "dust.in"], case, "dust.log", env)
        (case / "Postpro").mkdir(exist_ok=True)
        run([str(bins["post"]), "dust_post.in"], case, "dust_post.log", env)
        try:
            loads = parse_loads(case, args.velocity)
        except RuntimeError as error:
            results[str(n_chord)] = {
                "status": "invalid",
                "reason": str(error),
            }
            print(f"chord={n_chord}: INVALID ({error})", flush=True)
            continue
        results[str(n_chord)] = {"status": "ok", **loads}
        print(f"chord={n_chord}: Fz={loads['Fz_lbf']:+.6f} lbf, "
              f"My={loads['My_lbfin']:+.6f} lbf in", flush=True)
    changes = {}
    for previous, current in zip(CHORD_LEVELS, CHORD_LEVELS[1:]):
        previous_result = results[str(previous)]
        current_result = results[str(current)]
        if previous_result["status"] != "ok" or current_result["status"] != "ok":
            changes[f"{previous}_to_{current}"] = {"status": "not_available"}
            continue
        changes[f"{previous}_to_{current}"] = {
            key: relative_change(float(current_result[key]), float(previous_result[key]))
            for key in ("Fz_lbf", "My_lbfin", "CZ", "CMY_CG")
        }
    summary = {
        "geometry_source": "workflows/TRIM_DUST_1/case/mesh_canonical.in",
        "span_elements_by_region": list(SPAN_LEVELS),
        "chord_elements": list(CHORD_LEVELS),
        "production_chord_elements": 30,
        "selection": (
            "30 elementi di corda: 40 produce carichi NaN; 35 aumenta ancora "
            "l'aspect ratio dei pannelli e non mostra convergenza monotona del momento."
        ),
        "standalone_element_type": "vortex lattice",
        "angle_deg": args.angle,
        "velocity_mps": args.velocity,
        "wake_development_s": args.end_time,
        "results": results,
        "relative_change_percent": changes,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(changes, indent=2))


if __name__ == "__main__":
    main()
