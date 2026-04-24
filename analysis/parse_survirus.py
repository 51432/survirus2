import os
import re
from typing import Dict, Tuple

import pandas as pd

from .utils import file_exists

MAIN_RE = re.compile(
    r"^ID=(?P<event_id>\d+)\s+"
    r"(?P<host_chr>[^:]+):(?P<host_strand>[+-])(?P<host_pos>\d+)\s+"
    r"(?P<virus_contig>[^:]+):(?P<virus_strand>[+-])(?P<virus_pos>\d+)\s+"
    r"SUPPORTING_PAIRS=(?P<supporting_pairs>[-0-9.]+)\s+"
    r"SPLIT_READS=(?P<split_reads>[-0-9.]+)\s+"
    r"HOST_PBS=(?P<host_pbs>[-0-9.]+)\s+"
    r"COVERAGE=(?P<coverage>[-0-9.]+)"
)
ALT_RE = re.compile(r"^ID=(?P<event_id>\d+)\s+(?P<host_chr>[^:]+):(?P<host_strand>[+-])(?P<host_pos>\d+)")


def _virus_type(virus_contig: str) -> str:
    core = virus_contig.split("|")[0]
    m = re.match(r"([A-Za-z]+\d+)", core)
    return m.group(1) if m else core


def parse_main_results(path: str, sample_id: str, subtype: str = "") -> pd.DataFrame:
    rows = []
    if not file_exists(path):
        return pd.DataFrame()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            m = MAIN_RE.match(line)
            if not m:
                continue
            d = m.groupdict()
            d["sample_id"] = sample_id
            d["subtype"] = subtype
            d["event_uid"] = f"{sample_id}|{d['event_id']}"
            d["virus_type"] = _virus_type(d["virus_contig"])
            rows.append(d)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    int_cols = ["host_pos", "virus_pos", "supporting_pairs", "split_reads"]
    float_cols = ["host_pbs", "coverage"]
    for c in int_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    for c in float_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def parse_alt_results(path: str) -> pd.DataFrame:
    if not file_exists(path):
        return pd.DataFrame(columns=["event_id", "n_alt_total", "n_alt_extra", "multi_mapping_risk"])
    counts: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            m = ALT_RE.match(line)
            if not m:
                continue
            e = m.group("event_id")
            counts[e] = counts.get(e, 0) + 1
    rows = []
    for event_id, n in counts.items():
        n_extra = max(n - 1, 0)
        rows.append({
            "event_id": str(event_id),
            "n_alt_total": n,
            "n_alt_extra": n_extra,
            "multi_mapping_risk": 1 if n_extra > 0 else 0,
        })
    return pd.DataFrame(rows)


def parse_remapped(path: str) -> pd.DataFrame:
    df = parse_main_results(path, sample_id="", subtype="")
    if df.empty:
        return pd.DataFrame(columns=["event_id", "remap_confirmed", "remap_supporting_pairs", "remap_split_reads"])
    out = df[["event_id", "supporting_pairs", "split_reads"]].copy()
    out = out.rename(columns={
        "supporting_pairs": "remap_supporting_pairs",
        "split_reads": "remap_split_reads",
    })
    out["event_id"] = out["event_id"].astype(str)
    out["remap_confirmed"] = 1
    return out


def merge_survirus_tables(main_df: pd.DataFrame, remap_df: pd.DataFrame, alt_df: pd.DataFrame) -> pd.DataFrame:
    df = main_df.copy()
    df["event_id"] = df["event_id"].astype(str)
    if not remap_df.empty:
        df = df.merge(remap_df, on="event_id", how="left")
    else:
        df["remap_confirmed"] = 0
        df["remap_supporting_pairs"] = pd.NA
        df["remap_split_reads"] = pd.NA

    if not alt_df.empty:
        df = df.merge(alt_df, on="event_id", how="left")
    for c, v in [("n_alt_total", 0), ("n_alt_extra", 0), ("multi_mapping_risk", 0), ("remap_confirmed", 0)]:
        if c not in df.columns:
            df[c] = v
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(v).astype(int)
    return df
