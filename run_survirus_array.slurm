#!/usr/bin/env bash
#SBATCH -J survirus
#SBATCH -p cpu1,cpu2
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH -e logs/slurm-%A_%a.err

set -euo pipefail

# Usage:
#   sbatch --array=1-N%K run_survirus_array.slurm samples.tsv results [--wgs]
#   (Recommended) use submit_survirus_array.sh for resource-aware auto tuning.
#
# Notes:
# - samples.tsv must have header with columns: sample_id, input_R1, input_R2
# - SLURM_ARRAY_TASK_ID is used as 1-based row index over data rows.

SAMPLES="${1:-samples.tsv}"
OUTDIR="${2:-results}"
WGS_FLAG="${3:-}"

if [[ -n "$WGS_FLAG" && "$WGS_FLAG" != "--wgs" ]]; then
  echo "[ERROR] Third argument must be --wgs (or omitted), got: $WGS_FLAG" >&2
  exit 1
fi

mkdir -p logs

python3 run_survirus_pipeline.py \
  --samples "$SAMPLES" \
  --outdir "$OUTDIR" \
  --slurm-array \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --python2 "python2" \
  --surveyor "surveyor.py" \
  --bwa "${BWA_MEM2:-bwa-mem2}" \
  --samtools "samtools" \
  --dust "${DUST_EXEC:-sdust}" \
  ${WGS_FLAG:+--wgs}
