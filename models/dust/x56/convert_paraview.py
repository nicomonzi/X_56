#!/usr/bin/env python3
import re
import struct
from pathlib import Path

SOURCE = Path("Postpro")
TARGET = SOURCE / "Paraview"


def offset(header, name):
    match = re.search(
        rb'<DataArray[^>]*Name="' + re.escape(name.encode()) + rb'"[^>]*offset="(\d+)"',
        header,
        re.IGNORECASE,
    )
    if not match:
        raise RuntimeError(f"DataArray {name} non trovato")
    return int(match.group(1))


def block(raw, position):
    length = struct.unpack_from("<I", raw, position)[0]
    return raw[position + 4 : position + 4 + length]


def values(raw, position, code):
    data = block(raw, position)
    size = struct.calcsize(code)
    return struct.unpack(f"<{len(data) // size}{code}", data)


def rows(data, floating=False):
    lines = []
    for start in range(0, len(data), 12):
        row = data[start : start + 12]
        if floating:
            row = (f"{value:.9g}" for value in row)
        else:
            row = (str(value) for value in row)
        lines.append("     " + " ".join(row))
    return "\n".join(lines)


def convert(source, target):
    document = source.read_bytes()
    marker = re.search(rb'<AppendedData\s+encoding="raw">\s*_', document)
    if not marker:
        raise RuntimeError(f"{source} non contiene VTK appended raw")
    closing = document.rfind(b"</AppendedData>")
    header = document[: marker.start()]
    raw = document[marker.end() : closing]

    dimensions = re.search(rb'NumberOfPoints="(\d+)" NumberOfCells="(\d+)"', header)
    if not dimensions:
        raise RuntimeError(f"{source} non contiene le dimensioni della mesh")
    number_of_points, number_of_cells = map(int, dimensions.groups())

    points = values(raw, offset(header, "Coordinates"), "f")
    connectivity = values(raw, offset(header, "connectivity"), "i")
    cell_offsets = values(raw, offset(header, "offsets"), "i")
    cell_types = values(raw, offset(header, "types"), "i")
    intensity = values(raw, offset(header, "Singularity_Intensity"), "f")

    xml = f'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
 <UnstructuredGrid>
  <Piece NumberOfPoints="{number_of_points}" NumberOfCells="{number_of_cells}">
   <Points>
    <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="ascii">
{rows(points, True)}
    </DataArray>
   </Points>
   <Cells>
    <DataArray type="Int32" Name="connectivity" format="ascii">
{rows(connectivity)}
    </DataArray>
    <DataArray type="Int32" Name="offsets" format="ascii">
{rows(cell_offsets)}
    </DataArray>
    <DataArray type="UInt8" Name="types" format="ascii">
{rows(cell_types)}
    </DataArray>
   </Cells>
   <CellData Scalars="scalars">
    <DataArray type="Float32" Name="Singularity_Intensity" format="ascii">
{rows(intensity, True)}
    </DataArray>
   </CellData>
  </Piece>
 </UnstructuredGrid>
</VTKFile>
'''
    target.write_text(xml)


def read_dt():
    settings = Path("dust.in").read_text()
    match = re.search(r"^\s*dt_out\s*=\s*([0-9.eE+-]+)", settings, re.MULTILINE)
    return float(match.group(1)) if match else 1.0


def main():
    files = sorted(SOURCE.glob("sim_paraview-*.vtu"))
    if not files:
        raise SystemExit("Nessun file VTK prodotto da dust_post")
    TARGET.mkdir(parents=True, exist_ok=True)
    for source in files:
        convert(source, TARGET / source.name)

    datasets = "\n".join(
        f'  <DataSet timestep="{index * read_dt():.9g}" group="" part="0" file="Paraview/{source.name}"/>'
        for index, source in enumerate(files)
    )
    collection = f'''<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
 <Collection>
{datasets}
 </Collection>
</VTKFile>
'''
    (SOURCE / "sim.pvd").write_text(collection)
    print(f"Creato Postpro/sim.pvd con {len(files)} risultati ParaView")


if __name__ == "__main__":
    main()
