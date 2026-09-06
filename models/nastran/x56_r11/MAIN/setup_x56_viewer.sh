#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
venv_dir="${script_dir}/.venv_x56"

python3 -m venv "${venv_dir}"
"${venv_dir}/bin/python" -m pip install --upgrade pip
"${venv_dir}/bin/python" -m pip install pyNastran

echo
echo "Installazione completata."
echo "Avvia il viewer con: ${script_dir}/run_x56_viewer.sh"
