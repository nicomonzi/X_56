#!/usr/bin/env python3
"""Replace dust_post's fixed surface coordinates with saved coupled DUST geometry."""
from __future__ import annotations

import argparse
import base64
from pathlib import Path
import re
import struct

import h5py
import numpy as np


COORDINATES = re.compile(
    rb'(<DataArray\s+type="Float32"\s+Name="Coordinates"\s+'
    rb'NumberOfComponents="3"\s+format="binary">)([^<]*)(</DataArray>)',
    re.DOTALL,
)


def inject(vtu: Path, result: Path) -> None:
    with h5py.File(result, "r") as source:
        coordinates = np.asarray(
            source["Components/Comp001/Geometry/rr"], dtype="<f4"
        )
    payload = coordinates.tobytes(order="C")
    block = struct.pack("<I", len(payload)) + payload
    encoded = base64.b64encode(block)
    document = vtu.read_bytes()
    document, count = COORDINATES.subn(
        lambda match: match.group(1) + encoded + match.group(3), document, count=1
    )
    if count != 1:
        raise RuntimeError(f"Coordinates array not found exactly once: {vtu}")
    vtu.write_bytes(document)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paraview_directory", type=Path)
    parser.add_argument("results_directory", type=Path)
    args = parser.parse_args()
    files = sorted(args.paraview_directory.glob("x56_coupled-[0-9][0-9][0-9][0-9].vtu"))
    if not files:
        raise SystemExit("No coupled surface VTUs found")
    for vtu in files:
        index = vtu.stem.rsplit("-", 1)[1]
        result = args.results_directory / f"case_res_{index}.h5"
        if not result.is_file():
            raise FileNotFoundError(result)
        inject(vtu, result)
    print(f"Injected coupled DUST surface coordinates: {len(files)} VTUs")


if __name__ == "__main__":
    main()
