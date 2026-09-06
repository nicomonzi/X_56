#!/usr/bin/env python3
"""Generate DUST meshes, audit geometry, and export the FINE mesh to VTU."""
from __future__ import annotations

from collections import Counter
import json
import os
from pathlib import Path
import shutil
import subprocess

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LEVELS = {"COARSE": 900, "MEDIUM": 2280, "FINE": 4284}
WORK = ROOT / "reports/generated_meshes"
HINGE_Y = np.array([16.83, 45.93, 54.2175, 79.5, 104.785, 130.07, 155.35])
REGION_Y = np.array([0.0, 3.038, 9.35, 19.021, 25.665, 50.0, 168.0])


def executable() -> Path:
    value = os.environ.get("DUST_PRE_BIN") or shutil.which("dust_pre")
    if not value:
        raise SystemExit("DUST_PRE_BIN is not set and dust_pre is not on PATH")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise SystemExit(f"DUST preprocessor not found: {path}")
    return path


def build(level: str, dust_pre: Path) -> Path:
    WORK.mkdir(parents=True, exist_ok=True)
    config = WORK / f"dust_pre_{level.lower()}.in"
    output = WORK / f"{level.lower()}_geo.h5"
    config.write_text(
        "comp_name = X56\n"
        f"geo_file = meshes/{level}/parametric_mesh.in\n"
        "ref_tag = centerbody\n"
        f"file_name = {output.relative_to(ROOT)}\n"
    )
    subprocess.run([str(dust_pre), str(config.relative_to(ROOT))],
                   cwd=ROOT, check=True, capture_output=True, text=True)
    return output


def polygon_metrics(points: np.ndarray, elements: np.ndarray):
    quads = points[elements]
    cross1 = np.cross(quads[:, 1] - quads[:, 0], quads[:, 2] - quads[:, 0])
    cross2 = np.cross(quads[:, 2] - quads[:, 0], quads[:, 3] - quads[:, 0])
    areas = 0.5 * (np.linalg.norm(cross1, axis=1) + np.linalg.norm(cross2, axis=1))
    normals = cross1 + cross2
    lengths = np.stack([
        np.linalg.norm(quads[:, (i + 1) % 4] - quads[:, i], axis=1)
        for i in range(4)
    ], axis=1)
    aspect = lengths.max(axis=1) / np.maximum(lengths.min(axis=1), 1.e-14)
    centroids = quads.mean(axis=1)
    return areas, aspect, normals, centroids


def quantized_rows(values: np.ndarray, tolerance: float = 1.e-7) -> set[tuple[int, ...]]:
    return {tuple(row) for row in np.rint(values / tolerance).astype(np.int64)}


def audit(path: Path, expected: int) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path) as h5:
        geometry = h5["Components/Comp001/Geometry"]
        points = np.asarray(geometry["rr"], dtype=float)
        elements = np.asarray(geometry["ee"], dtype=int) - 1
        coupling = np.asarray(geometry["CouplingNodes"], dtype=float)
    areas, aspect, normals, centroids = polygon_metrics(points, elements)
    mirrored = points.copy(); mirrored[:, 1] *= -1
    symmetry_missing = len(quantized_rows(mirrored) - quantized_rows(points))
    element_keys = [tuple(sorted(row.tolist())) for row in elements]
    duplicate_elements = len(element_keys) - len(set(element_keys))
    edges = Counter(tuple(sorted((int(row[i]), int(row[(i + 1) % 4]))))
                    for row in elements for i in range(4))
    nonmanifold_edges = sum(count > 2 for count in edges.values())
    boundary_edges = sum(count == 1 for count in edges.values())
    y_values = np.abs(points[:, 1])
    hinge_errors = {f"{value:g}": float(np.min(np.abs(y_values - value)))
                    for value in HINGE_Y}
    region_errors = {f"{value:g}": float(np.min(np.abs(y_values - value)))
                     for value in REGION_Y}
    normal_norm = np.linalg.norm(normals, axis=1)
    # For the symmetric zero-camber sections, outward upper/lower normals must
    # agree with the sign of the local thickness coordinate.
    outward = []
    rounded_y = np.round(np.abs(centroids[:, 1]), 7)
    for value in np.unique(rounded_y):
        indices = np.flatnonzero(rounded_y == value)
        mid_z = np.median(centroids[indices, 2])
        useful = indices[np.abs(centroids[indices, 2] - mid_z) > 1.e-9]
        outward.extend((normals[useful, 2] * (centroids[useful, 2] - mid_z) > 0).tolist())
    hinge_alignment_error = max(hinge_errors.values())
    result = {
        "elements": int(len(elements)), "expected_elements": expected,
        "points": int(len(points)), "coupling_nodes": int(len(coupling)),
        "finite": bool(np.isfinite(points).all() and np.isfinite(normals).all()),
        "left_right_symmetry_missing_points": symmetry_missing,
        "duplicate_elements": duplicate_elements,
        "boundary_edges": boundary_edges,
        "nonmanifold_edges": nonmanifold_edges,
        "minimum_panel_area_in2": float(areas.min()),
        "maximum_panel_area_in2": float(areas.max()),
        "minimum_aspect_ratio": float(aspect.min()),
        "mean_aspect_ratio": float(aspect.mean()),
        "maximum_aspect_ratio": float(aspect.max()),
        "zero_normal_panels": int(np.count_nonzero(normal_norm < 1.e-10)),
        "outward_normal_fraction": float(np.mean(outward)),
        "hinge_endpoint_nearest_y_error_in": hinge_errors,
        "region_boundary_nearest_y_error_in": region_errors,
        "free_stream_axis": "+X",
        "winglet_present": False,
        "hinge_boundaries_aligned": bool(hinge_alignment_error < 1.e-6),
        "topology_source": "reconstructed discretization of available no-winglet parametric geometry",
    }
    result["basic_geometry_pass"] = bool(
        len(elements) == expected and result["finite"] and symmetry_missing == 0
        and duplicate_elements == 0 and nonmanifold_edges == 0
        and result["zero_normal_panels"] == 0
        and result["outward_normal_fraction"] > 0.98
        and result["hinge_boundaries_aligned"]
    )
    return result, points, elements, aspect, coupling


def write_vtu(path: Path, points: np.ndarray, elements: np.ndarray,
              aspect: np.ndarray) -> None:
    connectivity = " ".join(str(int(v)) for v in elements.ravel())
    offsets = " ".join(str(4 * (i + 1)) for i in range(len(elements)))
    coordinates = " ".join(f"{v:.10g}" for v in points.ravel())
    aspects = " ".join(f"{v:.8g}" for v in aspect)
    path.write_text(f'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid><Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(elements)}">
    <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{coordinates}</DataArray></Points>
    <Cells>
      <DataArray type="Int32" Name="connectivity" format="ascii">{connectivity}</DataArray>
      <DataArray type="Int32" Name="offsets" format="ascii">{offsets}</DataArray>
      <DataArray type="UInt8" Name="types" format="ascii">{' '.join(['9'] * len(elements))}</DataArray>
    </Cells>
    <CellData Scalars="panel_aspect_ratio"><DataArray type="Float64" Name="panel_aspect_ratio" format="ascii">{aspects}</DataArray></CellData>
  </Piece></UnstructuredGrid>
</VTKFile>
''')


def write_nodes_vtu(path: Path, points: np.ndarray) -> None:
    coordinates = " ".join(f"{v:.10g}" for v in points.ravel())
    ids = " ".join(str(i) for i in range(len(points)))
    offsets = " ".join(str(i + 1) for i in range(len(points)))
    path.write_text(f'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid><Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(points)}">
    <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{coordinates}</DataArray></Points>
    <Cells>
      <DataArray type="Int32" Name="connectivity" format="ascii">{ids}</DataArray>
      <DataArray type="Int32" Name="offsets" format="ascii">{offsets}</DataArray>
      <DataArray type="UInt8" Name="types" format="ascii">{' '.join(['1'] * len(points))}</DataArray>
    </Cells>
  </Piece></UnstructuredGrid>
</VTKFile>
''')


def main() -> None:
    dust_pre = executable()
    report = {"classification": "production diagnostic mesh with residual convergence uncertainty",
              "production_blocked": True, "levels": {}}
    fine_data = None
    for level, expected in LEVELS.items():
        result = audit(build(level, dust_pre), expected)
        report["levels"][level] = result[0]
        if level == "FINE": fine_data = result[1:]
    report["blocking_reasons"] = [
        "Available parametric geometry has no winglet component, so the wing/winglet junction cannot be audited.",
        "Control-surface hinge endpoints do not all coincide with spanwise mesh lines in the reconstructed FINE mesh.",
        "The original mesh-study input files were absent; the three discretizations were reconstructed from the available topology.",
        "Given MEDIUM-to-FINE changes (Fz 3.16%, My 14.59%, spanwise 26.33%), FINE is not formally converged.",
        "The requested modal count and requested highest frequency are inconsistent with the available FEM basis.",
        "Rigid DUST control-effectiveness derivatives have not yet been executed.",
    ]
    (ROOT / "reports/fine_mesh_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    assert fine_data is not None
    write_vtu(ROOT / "reports/fine_mesh_audit.vtu", *fine_data[:3])
    write_nodes_vtu(ROOT / "reports/fine_coupling_nodes.vtu", fine_data[3])
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
