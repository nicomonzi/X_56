#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
MBDYN_CMD="${MBDYN_CMD:-/usr/local/mbdyn/bin/mbdyn}"
"${MBDYN_CMD}" -f modal_mode_check.mbd -o mode_check

echo "Created mode_check.mov and the other MBDyn output files."
echo "Edit MODE_TO_EXCITE in modal_mode_check.mbd to inspect modes 7 through 18."
