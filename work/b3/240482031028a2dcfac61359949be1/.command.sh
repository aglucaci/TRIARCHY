#!/bin/bash -euo pipefail
ROOT="/Users/agl4001/Documents/TRIVIUM/workflow/.."
CONFIG="/Users/agl4001/Documents/TRIVIUM/config/test.config.yaml"
if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${ROOT}/${CONFIG}"
fi
SIGNALS="${PWD}/signal_classification.tsv"
CONCORDANCE="${PWD}/concordance.tsv"
OUT="${PWD}/report.html"
cd "${ROOT}"
python3 workflow/scripts/build_report.py       --config "${CONFIG}"       --signals "${SIGNALS}"       --concordance "${CONCORDANCE}"       --output "${OUT}"
