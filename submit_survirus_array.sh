#!/usr/bin/env bash
set -euo pipefail

# Resource-aware submit helper for SurVirus array jobs.
# It reads current idle CPUs from Slurm and computes a safe array concurrency.

usage() {
  cat <<'EOF'
Usage:
  ./submit_survirus_array.sh --samples samples.tsv --outdir results [options]

Options:
  --samples FILE            samples.tsv with header (required)
  --outdir DIR              output root dir passed to array script (required)
  --partitions P1,P2        Slurm partitions to query/submit (default: cpu1,cpu2)
  --fraction N              use at most 1/N of currently idle CPUs (default: 3)
  --threads-per-task N      cpus per sample task (default: 6)
  --mem-gb N                memory per task in GB (default: 32)
  --max-parallel N          hard cap for array concurrency (optional)
  --dry-run                 print computed sbatch command, do not submit

Examples:
  ./submit_survirus_array.sh --samples samples.tsv --outdir results --fraction 3
  ./submit_survirus_array.sh --samples samples.tsv --outdir results --fraction 2 --threads-per-task 8
EOF
}

SAMPLES=""
OUTDIR=""
PARTITIONS="cpu1,cpu2"
FRACTION=3
THREADS=6
MEM_GB=32
MAX_PARALLEL=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --partitions) PARTITIONS="$2"; shift 2 ;;
    --fraction) FRACTION="$2"; shift 2 ;;
    --threads-per-task) THREADS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$SAMPLES" ]] || { echo "[ERROR] --samples is required" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "[ERROR] --outdir is required" >&2; exit 1; }
[[ -f "$SAMPLES" ]] || { echo "[ERROR] samples.tsv not found: $SAMPLES" >&2; exit 1; }

[[ "$FRACTION" =~ ^[0-9]+$ ]] || { echo "[ERROR] --fraction must be positive integer" >&2; exit 1; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "[ERROR] --threads-per-task must be positive integer" >&2; exit 1; }
[[ "$MEM_GB" =~ ^[0-9]+$ ]] || { echo "[ERROR] --mem-gb must be positive integer" >&2; exit 1; }
(( FRACTION >= 2 )) || { echo "[ERROR] --fraction should be >=2 (shared cluster safety)" >&2; exit 1; }
(( THREADS >= 1 )) || { echo "[ERROR] --threads-per-task must be >=1" >&2; exit 1; }
(( MEM_GB >= 1 )) || { echo "[ERROR] --mem-gb must be >=1" >&2; exit 1; }

NUM_SAMPLES=$(( $(wc -l < "$SAMPLES") - 1 ))
(( NUM_SAMPLES >= 1 )) || { echo "[ERROR] samples.tsv has no data rows" >&2; exit 1; }

if ! command -v sinfo >/dev/null 2>&1; then
  echo "[ERROR] sinfo not found; cannot auto-tune by live idle CPUs" >&2
  exit 1
fi

IDLE_CPUS=$(sinfo -p "$PARTITIONS" -h -o "%C" | awk -F'/' '{idle+=$2} END{print idle+0}')
if [[ -z "$IDLE_CPUS" || "$IDLE_CPUS" -le 0 ]]; then
  echo "[WARN] No idle CPUs found in partitions=$PARTITIONS, forcing minimal parallelism=1" >&2
  IDLE_CPUS=1
fi

SAFE_CPUS=$(( IDLE_CPUS / FRACTION ))
(( SAFE_CPUS >= 1 )) || SAFE_CPUS=1
PARALLEL=$(( SAFE_CPUS / THREADS ))
(( PARALLEL >= 1 )) || PARALLEL=1

if (( MAX_PARALLEL > 0 && PARALLEL > MAX_PARALLEL )); then
  PARALLEL=$MAX_PARALLEL
fi
if (( PARALLEL > NUM_SAMPLES )); then
  PARALLEL=$NUM_SAMPLES
fi

ARRAY_SPEC="1-${NUM_SAMPLES}%${PARALLEL}"

CMD=(
  sbatch
  --partition "$PARTITIONS"
  --cpus-per-task "$THREADS"
  --mem "${MEM_GB}G"
  --array "$ARRAY_SPEC"
  run_survirus_array.slurm
  "$SAMPLES"
  "$OUTDIR"
)

echo "[INFO] samples        : $NUM_SAMPLES"
echo "[INFO] partitions     : $PARTITIONS"
echo "[INFO] idle_cpus      : $IDLE_CPUS"
echo "[INFO] fraction       : 1/$FRACTION"
echo "[INFO] safe_cpus      : $SAFE_CPUS"
echo "[INFO] threads/task   : $THREADS"
echo "[INFO] mem/task       : ${MEM_GB}G"
echo "[INFO] array_spec     : $ARRAY_SPEC"
echo "[INFO] submit command : ${CMD[*]}"

if (( DRY_RUN == 1 )); then
  exit 0
fi

"${CMD[@]}"
