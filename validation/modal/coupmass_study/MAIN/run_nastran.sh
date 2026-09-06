#!/usr/bin/env bash
set -euo pipefail
solver="${NASTRAN_CMD:-nastran}"
"$solver" sol103_coupmass_lumped.bdf scr=yes old=no
"$solver" sol103_coupmass_coupled.bdf scr=yes old=no
