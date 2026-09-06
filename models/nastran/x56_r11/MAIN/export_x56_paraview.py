#!/usr/bin/env python3
"""Esporta geometria e modi X-56 in VTU/PVD per ParaView.

Non richiede il pacchetto Python VTK: scrive direttamente VTK XML ASCII.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
from pyNastran.op2.op2 import read_op2

from plot_x56_modes import (
    DEFAULT_BDF,
    DEFAULT_OP2,
    eigenvector_result,
    mode_displacements,
    node_positions,
    read_x56_bdf,
    result_modes_and_frequencies,
)


DEFAULT_OUTPUT = Path(__file__).resolve().parents[3] / "results/figures/paraview"

# Tipo Nastran: (numero nodi usati, VTK cell type)
SHELL_CELL_TYPES = {
    "CTRIA3": (3, 5),   # VTK_TRIANGLE
    "CTRIAR": (3, 5),
    "CTRIA6": (6, 22),  # VTK_QUADRATIC_TRIANGLE
    "CQUAD4": (4, 9),   # VTK_QUAD
    "CQUADR": (4, 9),
    "CSHEAR": (4, 9),
    "CQUAD8": (8, 23),  # VTK_QUADRATIC_QUAD
}


def fmt_floats(values: np.ndarray) -> str:
    return " ".join(f"{value:.9e}" for value in values.ravel())


def fmt_ints(values: np.ndarray) -> str:
    return " ".join(str(int(value)) for value in values.ravel())


def shell_cells(bdf, positions: dict[int, np.ndarray]):
    cells: list[tuple[list[int], int, int, int]] = []
    used_nodes: set[int] = set()
    for eid, element in sorted(bdf.elements.items()):
        definition = SHELL_CELL_TYPES.get(element.type)
        if definition is None:
            continue
        node_count, vtk_type = definition
        try:
            node_ids = [int(nid) for nid in element.node_ids[:node_count]]
        except (TypeError, ValueError):
            continue
        if len(node_ids) != node_count or any(nid not in positions for nid in node_ids):
            continue
        cells.append((node_ids, vtk_type, int(eid), int(element.Pid())))
        used_nodes.update(node_ids)
    if not cells:
        raise RuntimeError("Nessun elemento shell esportabile trovato.")
    return cells, sorted(used_nodes)


def write_vtu(
    output: Path,
    points: np.ndarray,
    node_ids: np.ndarray,
    connectivity: np.ndarray,
    offsets: np.ndarray,
    vtk_types: np.ndarray,
    element_ids: np.ndarray,
    property_ids: np.ndarray,
    displacement: np.ndarray,
    normalized: np.ndarray,
    mode: int,
    frequency: float,
    eigenvalue: float,
) -> None:
    magnitude = np.linalg.norm(normalized, axis=1)
    with output.open("w", encoding="utf-8") as stream:
        stream.write('<?xml version="1.0"?>\n')
        stream.write(
            '<VTKFile type="UnstructuredGrid" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        stream.write("  <UnstructuredGrid>\n")
        stream.write("    <FieldData>\n")
        stream.write(
            f'      <DataArray type="Int32" Name="ModeNumber" '
            f'NumberOfTuples="1" format="ascii">{mode}</DataArray>\n'
        )
        stream.write(
            f'      <DataArray type="Float64" Name="Frequency_Hz" '
            f'NumberOfTuples="1" format="ascii">{frequency:.12e}</DataArray>\n'
        )
        stream.write(
            f'      <DataArray type="Float64" Name="Eigenvalue" '
            f'NumberOfTuples="1" format="ascii">{eigenvalue:.12e}</DataArray>\n'
        )
        stream.write("    </FieldData>\n")
        stream.write(
            f'    <Piece NumberOfPoints="{len(points)}" '
            f'NumberOfCells="{len(offsets)}">\n'
        )
        stream.write('      <PointData Vectors="ModeShape_Normalized">\n')
        stream.write(
            '        <DataArray type="Int64" Name="NodeID" format="ascii">\n'
            f"          {fmt_ints(node_ids)}\n"
            "        </DataArray>\n"
        )
        stream.write(
            '        <DataArray type="Float64" Name="ModeShape" '
            'NumberOfComponents="3" format="ascii">\n'
            f"          {fmt_floats(displacement)}\n"
            "        </DataArray>\n"
        )
        stream.write(
            '        <DataArray type="Float64" Name="ModeShape_Normalized" '
            'NumberOfComponents="3" format="ascii">\n'
            f"          {fmt_floats(normalized)}\n"
            "        </DataArray>\n"
        )
        stream.write(
            '        <DataArray type="Float64" Name="NormalizedMagnitude" '
            'format="ascii">\n'
            f"          {fmt_floats(magnitude)}\n"
            "        </DataArray>\n"
        )
        stream.write("      </PointData>\n")
        stream.write("      <CellData>\n")
        stream.write(
            '        <DataArray type="Int64" Name="ElementID" format="ascii">\n'
            f"          {fmt_ints(element_ids)}\n"
            "        </DataArray>\n"
        )
        stream.write(
            '        <DataArray type="Int64" Name="PropertyID" format="ascii">\n'
            f"          {fmt_ints(property_ids)}\n"
            "        </DataArray>\n"
        )
        stream.write("      </CellData>\n")
        stream.write("      <Points>\n")
        stream.write(
            '        <DataArray type="Float64" Name="Coordinates_in" '
            'NumberOfComponents="3" format="ascii">\n'
            f"          {fmt_floats(points)}\n"
            "        </DataArray>\n"
        )
        stream.write("      </Points>\n")
        stream.write("      <Cells>\n")
        for name, vtk_name, values, vtk_dtype in (
            ("connectivity", "connectivity", connectivity, "Int64"),
            ("offsets", "offsets", offsets, "Int64"),
            ("types", "types", vtk_types, "UInt8"),
        ):
            stream.write(
                f'        <DataArray type="{vtk_dtype}" Name="{vtk_name}" '
                f'format="ascii">\n          {fmt_ints(values)}\n'
                "        </DataArray>\n"
            )
        stream.write("      </Cells>\n")
        stream.write("    </Piece>\n")
        stream.write("  </UnstructuredGrid>\n")
        stream.write("</VTKFile>\n")


def write_pvd(
    output: Path, entries: list[tuple[int, float, str]]
) -> None:
    with output.open("w", encoding="utf-8") as stream:
        stream.write('<?xml version="1.0"?>\n')
        stream.write(
            '<VTKFile type="Collection" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        stream.write("  <Collection>\n")
        for mode, frequency, filename in entries:
            name = escape(f"Mode {mode} - {frequency:.3f} Hz")
            stream.write(
                f'    <DataSet timestep="{mode}" group="" part="0" '
                f'name="{name}" file="{escape(filename)}"/>\n'
            )
        stream.write("  </Collection>\n")
        stream.write("</VTKFile>\n")


def write_readme(output: Path) -> None:
    text = """# X-56 modal shapes for ParaView

1. Open `x56_modes.pvd` in ParaView and click **Apply**.
2. Select a mode with the time controls; timestep equals the Nastran mode.
3. Apply **Filters > Alphabetical > Warp By Vector**.
4. Set **Vectors** to `ModeShape_Normalized`.
5. Set **Scale Factor** in inches (for example 30 or 50), then click **Apply**.
6. Set **Representation** to `Wireframe` or `Surface With Edges`.

Modes 1-6 are the numerical free-free rigid-body modes from Nastran.
Elastic modes start at mode 7 (3.217 Hz).

Arrays:
- `ModeShape`: original mass-normalized Nastran eigenvector.
- `ModeShape_Normalized`: same shape normalized to a maximum of 1 inch.
- `NormalizedMagnitude`: magnitude useful for coloring.
- `NodeID`, `ElementID`, `PropertyID`: original Nastran identifiers.

All coordinates and Warp scale factors are in inches.
"""
    (output / "README.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bdf", type=Path, default=DEFAULT_BDF)
    parser.add_argument("--op2", type=Path, default=DEFAULT_OP2)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Leggo BDF: {args.bdf}")
    bdf = read_x56_bdf(args.bdf)
    print(f"Leggo OP2: {args.op2}")
    op2 = read_op2(
        str(args.op2), combine=True, build_dataframe=False, debug=False
    )
    result = eigenvector_result(op2)
    modes, frequencies = result_modes_and_frequencies(result)
    eigenvalues = np.asarray(result.eigns, dtype=float)

    positions = node_positions(bdf)
    cells, used_node_ids = shell_cells(bdf, positions)
    node_index = {nid: index for index, nid in enumerate(used_node_ids)}
    points = np.vstack([positions[nid] for nid in used_node_ids])

    connectivity_list: list[int] = []
    offsets: list[int] = []
    vtk_types: list[int] = []
    element_ids: list[int] = []
    property_ids: list[int] = []
    for nids, vtk_type, eid, pid in cells:
        connectivity_list.extend(node_index[nid] for nid in nids)
        offsets.append(len(connectivity_list))
        vtk_types.append(vtk_type)
        element_ids.append(eid)
        property_ids.append(pid)

    connectivity = np.asarray(connectivity_list, dtype=np.int64)
    offsets_array = np.asarray(offsets, dtype=np.int64)
    vtk_types_array = np.asarray(vtk_types, dtype=np.uint8)
    node_ids_array = np.asarray(used_node_ids, dtype=np.int64)
    element_ids_array = np.asarray(element_ids, dtype=np.int64)
    property_ids_array = np.asarray(property_ids, dtype=np.int64)

    print(
        f"Esporto {len(modes)} modi, {len(points)} nodi shell, "
        f"{len(cells)} elementi shell in {args.output_dir}"
    )
    pvd_entries: list[tuple[int, float, str]] = []
    for index, (mode, frequency) in enumerate(zip(modes, frequencies)):
        displacement_by_node = mode_displacements(result, bdf, index)
        displacement = np.vstack(
            [displacement_by_node.get(nid, np.zeros(3)) for nid in used_node_ids]
        )
        max_norm = float(np.max(np.linalg.norm(displacement, axis=1)))
        normalized = displacement / max_norm if max_norm > 0.0 else displacement
        filename = f"x56_mode_{int(mode):02d}.vtu"
        write_vtu(
            output=args.output_dir / filename,
            points=points,
            node_ids=node_ids_array,
            connectivity=connectivity,
            offsets=offsets_array,
            vtk_types=vtk_types_array,
            element_ids=element_ids_array,
            property_ids=property_ids_array,
            displacement=displacement,
            normalized=normalized,
            mode=int(mode),
            frequency=float(frequency),
            eigenvalue=float(eigenvalues[index]),
        )
        pvd_entries.append((int(mode), float(frequency), filename))
        print(f"  modo {int(mode):2d}: {float(frequency):8.3f} Hz -> {filename}")

    write_pvd(args.output_dir / "x56_modes.pvd", pvd_entries)
    write_readme(args.output_dir)
    print(f"Creato: {args.output_dir / 'x56_modes.pvd'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
