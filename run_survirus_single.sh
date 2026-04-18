#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_survirus_single.sh samples.tsv SAMPLE_ID results

SAMPLES="${1:?need samples.tsv}"
SAMPLE_ID="${2:?need sample_id}"
OUTDIR="${3:-results}"
THREADS="${4:-6}"

python3 run_survirus_pipeline.py \
  --samples "$SAMPLES" \
  --outdir "$OUTDIR" \
  --sample-id "$SAMPLE_ID" \
  --threads "$THREADS" \
  --python2 "python2" \
  --surveyor "surveyor.py" \
  --bwa "bwa" \
  --samtools "samtools" \
  --dust "dust"
