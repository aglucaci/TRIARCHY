#!/bin/bash -euo pipefail
ROOT="/Users/agl4001/Documents/TRIVIUM/workflow/.."
SAMPLESHEET="test_data/metadata/samplesheet.csv"
if [[ "${SAMPLESHEET}" != /* ]]; then
  SAMPLESHEET="${ROOT}/${SAMPLESHEET}"
fi
python3 "${ROOT}/workflow/scripts/validate_samplesheet.py"       --samplesheet "${SAMPLESHEET}"       > samplesheet.validated.txt
