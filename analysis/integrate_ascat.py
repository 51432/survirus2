import os
from typing import Optional

import pandas as pd


def _read_ascat(path: Optional[str]) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str)
    colmap = {c.lower(): c for c in df.columns}
    needed = ["chr", "chromosome", "start", "end", "major_cn", "minor_cn", "total_cn"]
    # normalize names heuristically
    for src in ["chr", "chromosome"]:
        if src in colmap:
            df = df.rename(columns={colmap[src]: "chrom"})
            break
    for src, dst in [("start", "start"), ("end", "end"), ("major_cn", "major_cn"), ("minor_cn", "minor_cn"), ("total_cn", "total_cn")]:
        if src in colmap:
            df = df.rename(columns={colmap[src]: dst})
    for c in ["start", "end", "major_cn", "minor_cn", "total_cn"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def annotate_local_cnv(events: pd.DataFrame, manifest_row: pd.Series) -> pd.DataFrame:
    out = events.copy()
    out["segment_start"] = pd.NA
    out["segment_end"] = pd.NA
    out["local_total_cn"] = pd.NA
    out["local_major_cn"] = pd.NA
    out["local_minor_cn"] = pd.NA
    out["sample_ploidy"] = manifest_row.get("sample_ploidy", "NA") or "NA"
    out["sample_wgd_status"] = manifest_row.get("sample_wgd_status", "NA") or "NA"

    ascat_df = _read_ascat(manifest_row.get("ascat_segments", ""))
    if ascat_df.empty or "chrom" not in ascat_df.columns:
        return out

    for idx, r in out.iterrows():
        chrom = str(r["host_chr"])
        pos = int(r["host_pos"])
        seg = ascat_df[(ascat_df["chrom"].astype(str) == chrom) & (ascat_df["start"] <= pos) & (ascat_df["end"] >= pos)]
        if seg.empty:
            continue
        s = seg.iloc[0]
        out.at[idx, "segment_start"] = s.get("start", pd.NA)
        out.at[idx, "segment_end"] = s.get("end", pd.NA)
        out.at[idx, "local_total_cn"] = s.get("total_cn", pd.NA)
        out.at[idx, "local_major_cn"] = s.get("major_cn", pd.NA)
        out.at[idx, "local_minor_cn"] = s.get("minor_cn", pd.NA)
    return out
