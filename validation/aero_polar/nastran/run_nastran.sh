#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
"${NASTRAN_CMD:-nastran}" x56_polar.bdf scr=yes old=no
