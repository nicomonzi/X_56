#!/usr/bin/env bash
set -euo pipefail
solver="${NASTRAN_CMD:-nastran}"
"$solver" sol101_gravity_5g.bdf scr=yes old=no
