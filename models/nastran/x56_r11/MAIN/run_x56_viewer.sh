#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
venv_python="${script_dir}/.venv_x56/bin/python"

if [[ ! -x "${venv_python}" ]]; then
    echo "Ambiente Python X-56 non trovato."
    echo "Esegui prima: ${script_dir}/setup_x56_viewer.sh"
    exit 1
fi

export MPLCONFIGDIR="${script_dir}/.matplotlib"
mkdir -p -- "${MPLCONFIGDIR}"

exec "${venv_python}" "${script_dir}/view_x56_modes.py" "$@"
