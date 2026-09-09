#!/usr/bin/env python3
"""Convert DUST raw-appended VTU files to portable inline-Base64 VTU."""
from __future__ import annotations
import argparse
import base64
from pathlib import Path
import re
import struct

TAG = b'<AppendedData encoding="raw">'
END = b"  </AppendedData>"
PATTERN = re.compile(
    r'<DataArray(?P<a>[^>]*?)\sformat="appended"\s+offset="(?P<o>\d+)"\s*/>',
    re.IGNORECASE,
)


def convert(path: Path) -> bool:
    source = path.read_bytes(); tag = source.find(TAG)
    if tag < 0: return False
    start = source.find(b"_", tag + len(TAG)) + 1; stop = source.rfind(END)
    if start <= 0 or stop <= start: raise ValueError(f"Invalid appended section: {path}")
    header = source[:tag].decode("ascii"); payload = source[start:stop]
    def replacement(match):
        offset = int(match.group("o")); size = struct.unpack_from("<I", payload, offset)[0]
        block = payload[offset:offset+4+size]
        if len(block) != size+4: raise ValueError(f"Truncated block: {path}")
        return f'<DataArray{match.group("a")} format="binary">{base64.b64encode(block).decode("ascii")}</DataArray>'
    text, count = PATTERN.subn(replacement, header)
    if not count: raise ValueError(f"No appended arrays: {path}")
    temporary = path.with_suffix(".vtu.tmp"); temporary.write_text(text.rstrip()+"\n</VTKFile>\n")
    temporary.replace(path); return True


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("directory", type=Path)
    args = parser.parse_args(); files = sorted(args.directory.glob("*.vtu"))
    print(f"Portable VTU converted: {sum(convert(path) for path in files)}/{len(files)}")


if __name__ == "__main__": main()
