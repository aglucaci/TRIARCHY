#!/bin/bash -euo pipefail
ROOT="/Users/agl4001/Documents/TRIVIUM/workflow/.."
CONFIG="/Users/agl4001/Documents/TRIVIUM/config/test.config.yaml"
if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${ROOT}/${CONFIG}"
fi
OUT="${PWD}/concordance.tsv"
cd "${ROOT}"
python3 workflow/scripts/compute_concordance.py       --config "${CONFIG}"       --output "${OUT}"
