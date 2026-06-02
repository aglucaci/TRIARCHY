#!/bin/bash -euo pipefail
ROOT="/Users/agl4001/Documents/TRIVIUM/workflow/.."
CONFIG="/Users/agl4001/Documents/TRIVIUM/config/test.config.yaml"
if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${ROOT}/${CONFIG}"
fi
OUT="${PWD}/signal_classification.tsv"
cd "${ROOT}"
python3 workflow/scripts/classify_signals.py       --config "${CONFIG}"       --output "${OUT}"
