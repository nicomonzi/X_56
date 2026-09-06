#!/usr/bin/env bash
set -euo pipefail
solver="${NASTRAN_CMD:-nastran}"
"$solver" sol103_60_modes.bdf scr=yes old=no
