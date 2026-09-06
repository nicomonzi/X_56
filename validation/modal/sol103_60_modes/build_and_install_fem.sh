#!/usr/bin/env bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
femgen_cmd="${FEMGEN_CMD:-femgen}"
cd "$here/MAIN"
"$femgen_cmd" sol103_60_modes -o mbdyn_modal_60
if ! grep -q "RECORD GROUP 10" mbdyn_modal_60.fem; then
    echo "ERROR: incomplete FEM. The OP2 does not contain all modal vectors." >&2
    echo "Rerun sol103_60_modes.bdf with DISPLACEMENT(PLOT)=ALL and VECTOR(PLOT)=ALL." >&2
    exit 2
fi
cp mbdyn_modal_60.fem "$here/../03_GRAVITY_5G/mbdyn/mbdyn_modal.fem"
echo "Installed 60-mode FEM in 03_GRAVITY_5G/mbdyn/mbdyn_modal.fem"
