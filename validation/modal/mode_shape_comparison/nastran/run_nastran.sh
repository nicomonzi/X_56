#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
NASTRAN_CMD="${NASTRAN_CMD:-nastran}"
"${NASTRAN_CMD}" sol103_10lb_f06.bdf old=no mode=i8

test -s sol103_10lb_f06.f06
if grep -Eq "FATAL MESSAGE|USER FATAL" sol103_10lb_f06.f06; then
    echo "Nastran reported a fatal error; inspect sol103_10lb_f06.f06." >&2
    exit 1
fi

echo "Created sol103_10lb_f06.f06 (printed modes) and sol103_10lb_f06.op2."
echo "The ALTER also created mbdyn_modal.mat for femgen."
