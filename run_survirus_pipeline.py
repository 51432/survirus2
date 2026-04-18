#!/usr/bin/env python3
"""HPC-friendly wrapper for running SurVirus from samples.tsv."""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import Counter

DEFAULT_HOST = "/data/person/wup/public/liusy_files/reference_genomes/hg38_no_alt/reference/hg38_no_alt.fa"
DEFAULT_VIRUS = "/data/person/wup/public/liusy_files/reference_genomes/virus/reference/virus.fa"
DEFAULT_HOST_VIRUS = "/data/person/wup/public/liusy_files/reference_genomes/host_virus/reference/host_virus.fa"

REQUIRED_COLUMNS = ("sample_id", "input_R1", "input_R2")


def fail(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def read_samples_tsv(samples_path):
    if not os.path.exists(samples_path):
        fail(f"samples.tsv not found: {samples_path}")

    with open(samples_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            fail("samples.tsv is empty or missing header")

        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            fail(
                "samples.tsv missing required columns: "
                + ", ".join(missing)
                + f". Found columns: {reader.fieldnames}"
            )

        rows = []
        for idx, row in enumerate(reader, start=2):
            sample = {k: (row.get(k, "") or "").strip() for k in REQUIRED_COLUMNS}

            for key, value in sample.items():
                if not value:
                    fail(f"Line {idx}: column '{key}' is empty")

            rows.append(sample)

    if not rows:
        fail("samples.tsv has header but no sample rows")

    counts = Counter(r["sample_id"] for r in rows)
    dup = [sid for sid, c in counts.items() if c > 1]
    if dup:
        fail("Duplicate sample_id found: " + ", ".join(sorted(dup)))

    for row in rows:
        for key in ("input_R1", "input_R2"):
            if not os.path.exists(row[key]):
                fail(f"Sample '{row['sample_id']}': file does not exist: {row[key]}")

    return rows


def select_sample(rows, sample_id=None, array_task_id=None):
    if sample_id and array_task_id is not None:
        fail("Please use either --sample-id or --slurm-array/--task-id, not both")

    if sample_id:
        matches = [r for r in rows if r["sample_id"] == sample_id]
        if not matches:
            fail(f"sample_id '{sample_id}' not found in samples.tsv")
        return matches[0]

    if array_task_id is not None:
        if array_task_id < 1 or array_task_id > len(rows):
            fail(
                f"Array task id {array_task_id} out of range. "
                f"Valid range: 1..{len(rows)}"
            )
        return rows[array_task_id - 1]

    if len(rows) == 1:
        return rows[0]

    fail(
        "Multiple samples found. Please specify one of: --sample-id, "
        "--slurm-array (or --task-id), or provide a single-row samples.tsv"
    )


def ensure_output_dir(path, force=False):
    if os.path.exists(path):
        if not force:
            fail(
                f"Output directory already exists: {path}. "
                "Use --force to reuse the directory."
            )
    else:
        os.makedirs(path)


def resolve_aligner_execs(user_bwa, user_bwa_mem2):
    bwa_exec = user_bwa
    if shutil.which(bwa_exec) is None:
        fail(
            f"Classic bwa not found: {bwa_exec}. "
            "SurVirus still needs bwa aln/samse steps. Please set --bwa."
        )

    if user_bwa_mem2:
        if shutil.which(user_bwa_mem2) is None:
            fail(f"--bwa-mem2 not found in PATH: {user_bwa_mem2}")
        return bwa_exec, user_bwa_mem2, "user-specified --bwa-mem2"

    auto_mem2 = shutil.which("bwa-mem2")
    if auto_mem2:
        return bwa_exec, auto_mem2, "auto-detected bwa-mem2"

    return bwa_exec, bwa_exec, "fallback to bwa mem"


def main():
    parser = argparse.ArgumentParser(
        description="Run SurVirus from samples.tsv (single sample or Slurm array task)."
    )
    parser.add_argument("--samples", required=True, help="Path to samples.tsv")
    parser.add_argument("--outdir", required=True, help="Base output directory")

    parser.add_argument("--sample-id", default="", help="Run one sample by sample_id")
    parser.add_argument(
        "--slurm-array",
        action="store_true",
        help="Use SLURM_ARRAY_TASK_ID to select the sample row (1-based)",
    )
    parser.add_argument(
        "--task-id",
        type=int,
        default=None,
        help="Manually set task id (1-based), for local debugging",
    )

    parser.add_argument("--threads", type=int, default=6, help="Threads per sample")
    parser.add_argument("--python2", default="python2", help="Python2 executable for surveyor.py")
    parser.add_argument("--surveyor", default="surveyor.py", help="Path to surveyor.py")

    parser.add_argument("--host", default=DEFAULT_HOST, help="Host reference fasta")
    parser.add_argument("--virus", default=DEFAULT_VIRUS, help="Virus reference fasta")
    parser.add_argument("--host-virus", default=DEFAULT_HOST_VIRUS, help="Host+virus reference fasta")

    parser.add_argument("--bwa", default="bwa", help="bwa executable path")
    parser.add_argument(
        "--bwa-mem2",
        default="",
        help="Optional bwa-mem2 path. If omitted, auto-detect bwa-mem2 first, else fallback to bwa mem.",
    )
    parser.add_argument("--samtools", default="samtools", help="samtools executable path")
    parser.add_argument("--dust", default="dust", help="dust executable path")

    parser.add_argument("--force", action="store_true", help="Reuse existing sample output directory")
    parser.add_argument("--dry-run", action="store_true", help="Only print command, do not run")

    args = parser.parse_args()

    rows = read_samples_tsv(args.samples)

    task_id = args.task_id
    if args.slurm_array:
        env_tid = os.environ.get("SLURM_ARRAY_TASK_ID")
        if env_tid is None:
            fail("--slurm-array is set but SLURM_ARRAY_TASK_ID is not defined")
        try:
            task_id = int(env_tid)
        except ValueError:
            fail(f"Invalid SLURM_ARRAY_TASK_ID value: {env_tid}")

    sample = select_sample(rows, sample_id=args.sample_id or None, array_task_id=task_id)

    os.makedirs(args.outdir, exist_ok=True)
    sample_outdir = os.path.join(args.outdir, sample["sample_id"])
    ensure_output_dir(sample_outdir, force=args.force)
    bwa_exec, bwa_mem_exec, bwa_mode = resolve_aligner_execs(args.bwa, args.bwa_mem2)

    logs_dir = os.path.join(args.outdir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, f"{sample['sample_id']}.survirus.log")

    cmd = [
        args.python2,
        args.surveyor,
        f"{sample['input_R1']},{sample['input_R2']}",
        sample_outdir,
        args.host,
        args.virus,
        args.host_virus,
        "--fq",
        "--threads",
        str(args.threads),
        "--bwa",
        bwa_exec,
        "--bwa-mem",
        bwa_mem_exec,
        "--samtools",
        args.samtools,
        "--dust",
        args.dust,
    ]

    print(f"[INFO] Selected sample: {sample['sample_id']}")
    print(f"[INFO] Output dir: {sample_outdir}")
    print(f"[INFO] Log file: {log_path}")
    print(f"[INFO] BWA mode: {bwa_mode}")
    print(f"[INFO] bwa (aln/samse): {bwa_exec}")
    print(f"[INFO] bwa-mem engine: {bwa_mem_exec}")
    print("[INFO] Command:")
    print(" ".join(cmd))

    if args.dry_run:
        return

    with open(log_path, "w") as logf:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            print(line, end="")
            logf.write(line)
        code = proc.wait()

    if code != 0:
        fail(f"SurVirus failed for sample '{sample['sample_id']}'. Check log: {log_path}")

    print(f"[INFO] Completed sample: {sample['sample_id']}")


if __name__ == "__main__":
    main()
