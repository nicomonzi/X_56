#!/usr/bin/env python3
"""Repair DUST raw-appended VTU XML and create the trim animation PVD."""
from __future__ import annotations
import base64
from pathlib import Path
import re
import struct

ROOT = Path(__file__).resolve().parent
POST = ROOT / "work" / "current" / "paraview"
RAW_TAG = b'<AppendedData encoding="raw">'
RAW_END = b"  </AppendedData>"
ARRAY = re.compile(r'<DataArray(?P<attrs>[^>]*?)\sformat="appended"\s+'
                   r'offset="(?P<offset>\d+)"\s*/>', re.IGNORECASE)


def portable(path: Path) -> None:
    source = path.read_bytes()
    tag = source.find(RAW_TAG)
    if tag < 0:
        return
    start = source.find(b"_", tag + len(RAW_TAG)) + 1
    stop = source.rfind(RAW_END)
    if start <= 0 or stop <= start:
        raise RuntimeError(f"VTU appended non valido: {path}")
    header = source[:tag].decode("ascii")
    payload = source[start:stop]

    def replace(match: re.Match[str]) -> str:
        offset = int(match.group("offset"))
        size = struct.unpack_from("<I", payload, offset)[0]
        block = payload[offset:offset + 4 + size]
        return (f'<DataArray{match.group("attrs")} format="binary">'
                f'{base64.b64encode(block).decode("ascii")}</DataArray>')

    xml, count = ARRAY.subn(replace, header)
    if not count:
        raise RuntimeError(f"nessun array appended in {path}")
    path.write_text(xml.rstrip() + "\n</VTKFile>\n")


def main() -> None:
    files = sorted(POST.glob("trim55_coupled*.vtu"))
    if not files:
        raise SystemExit("VTU mancanti: esegui prima dust_post case/dust_post.in da work/current")
    for path in files:
        portable(path)
    stems = ("trim55_coupled", "trim55_coupled_wpan", "trim55_coupled_wpart")
    sequences = [sorted(POST.glob(f"{stem}-[0-9][0-9][0-9][0-9].vtu")) for stem in stems]
    rows = []
    for part, sequence in enumerate(sequences):
        for path in sequence:
            index = int(path.stem.rsplit("-", 1)[1])
            rows.append(((index - 1) * 0.04, part, path.name))
    rows.sort()
    target = POST / "trim55_all.pvd"
    target.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
        '  <Collection>\n' + "\n".join(
            f'    <DataSet timestep="{t:.9g}" group="" part="{part}" file="{name}"/>'
            for t, part, name in rows) + '\n  </Collection>\n</VTKFile>\n')
    print(f"Convertiti {len(files)} VTU")
    print(f"Apri in ParaView: {target}")


if __name__ == "__main__":
    main()
