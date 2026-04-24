import logging
import os
from typing import Dict, List, Optional, Tuple

import pandas as pd


BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]


def setup_logging(outdir: str, level: str = "INFO") -> None:
    os.makedirs(outdir, exist_ok=True)
    log_file = os.path.join(outdir, "analysis.log")
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )


def read_manifest(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"sample_id", "survirus_results_t1"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")
    if "subtype" not in df.columns:
        df["subtype"] = ""
    return df


def file_exists(path: str) -> bool:
    return bool(path) and os.path.exists(path)


def ensure_columns(df: pd.DataFrame, cols: List[str], fill_value="NA") -> pd.DataFrame:
    for c in cols:
        if c not in df.columns:
            df[c] = fill_value
    return df


def parse_bed(path: Optional[str], min_cols: int = 3) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=BED_COLS)
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < min_cols:
                continue
            while len(toks) < 6:
                toks.append(".")
            rows.append(toks[:6])
    if not rows:
        return pd.DataFrame(columns=BED_COLS)
    bed = pd.DataFrame(rows, columns=BED_COLS)
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce").fillna(0).astype(int)
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce").fillna(0).astype(int)
    return bed


def read_gene_gtf(path: Optional[str]) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["chrom", "start", "end", "gene", "strand"])
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = toks
            if feature != "gene":
                continue
            gene = ""
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("gene_name") or field.startswith("gene_id"):
                    parts = field.split(" ")
                    if len(parts) > 1:
                        gene = parts[1].strip('"')
                        break
            if not gene:
                gene = "UNKNOWN"
            rows.append((chrom, int(start) - 1, int(end), gene, strand))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "strand"])


def build_interval_index(df: pd.DataFrame) -> Dict[str, List[Tuple[int, int, dict]]]:
    idx: Dict[str, List[Tuple[int, int, dict]]] = {}
    if df.empty:
        return idx
    for _, r in df.iterrows():
        chrom = str(r["chrom"])
        meta = r.to_dict()
        idx.setdefault(chrom, []).append((int(r["start"]), int(r["end"]), meta))
    for chrom in idx:
        idx[chrom].sort(key=lambda x: x[0])
    return idx


def query_overlaps(index: Dict[str, List[Tuple[int, int, dict]]], chrom: str, pos: int) -> List[dict]:
    out = []
    for start, end, meta in index.get(str(chrom), []):
        if start <= pos < end:
            out.append(meta)
        elif start > pos:
            break
    return out


def nearest_interval(index: Dict[str, List[Tuple[int, int, dict]]], chrom: str, pos: int) -> Tuple[Optional[dict], Optional[int]]:
    best = None
    best_dist = None
    for start, end, meta in index.get(str(chrom), []):
        if start <= pos < end:
            return meta, 0
        dist = min(abs(pos - start), abs(pos - end))
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best = meta
        if start > pos and best_dist is not None and (start - pos) > best_dist:
            break
    return best, best_dist
