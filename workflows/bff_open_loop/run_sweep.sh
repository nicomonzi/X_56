#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
exec python3 run_sweep.py \
    --start 50 \
    --stop 70 \
    --step 2.5 \
    --refine-tolerance 0.25 \
    --output /mnt/c/Users/Utente/Desktop/BFF_open_loop \
    --clean \
    --overwrite "$@"
