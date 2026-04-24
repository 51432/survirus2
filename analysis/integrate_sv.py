import os
from typing import Optional

import pandas as pd


def _parse_manta(vcf_path: Optional[str]) -> pd.DataFrame:
    rows = []
    if not vcf_path or not os.path.exists(vcf_path):
        return pd.DataFrame(columns=["chrom", "pos", "svtype"])
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 8:
                continue
            chrom, pos, _, _, _, _, _, info = toks[:8]
            infod = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    infod[k] = v
            svtype = infod.get("SVTYPE", "NA")
            rows.append({"chrom": chrom, "pos": int(pos), "svtype": svtype})
            if "END" in infod:
                rows.append({"chrom": chrom, "pos": int(infod["END"]), "svtype": svtype})
    return pd.DataFrame(rows)


def _parse_tiddit(path: Optional[str]) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["chrom", "pos", "svtype"])
    df = pd.read_csv(path, sep="\t", dtype=str)
    col = {c.lower(): c for c in df.columns}
    chr_col = col.get("chr") or col.get("chrom")
    pos_col = col.get("pos") or col.get("start")
    type_col = col.get("svtype") or col.get("type")
    if not chr_col or not pos_col:
        return pd.DataFrame(columns=["chrom", "pos", "svtype"])
    out = pd.DataFrame({
        "chrom": df[chr_col],
        "pos": pd.to_numeric(df[pos_col], errors="coerce"),
        "svtype": df[type_col] if type_col else "NA",
    }).dropna(subset=["pos"])
    out["pos"] = out["pos"].astype(int)
    return out


def annotate_local_sv(events: pd.DataFrame, manta_sv: Optional[str], tiddit_sv: Optional[str]) -> pd.DataFrame:
    out = events.copy()
    sv = pd.concat([_parse_manta(manta_sv), _parse_tiddit(tiddit_sv)], ignore_index=True)

    for c in ["n_sv_breakpoints_50kb", "n_sv_breakpoints_100kb", "n_sv_breakpoints_500kb", "has_local_del", "has_local_dup", "has_local_inv", "has_local_tra", "local_complex_sv_flag"]:
        out[c] = 0

    if sv.empty:
        out[["n_sv_breakpoints_50kb", "n_sv_breakpoints_100kb", "n_sv_breakpoints_500kb"]] = pd.NA
        out[["has_local_del", "has_local_dup", "has_local_inv", "has_local_tra", "local_complex_sv_flag"]] = pd.NA
        return out

    sv["chrom"] = sv["chrom"].astype(str)
    sv["svtype"] = sv["svtype"].fillna("NA").str.upper()

    for idx, r in out.iterrows():
        chrom = str(r["host_chr"])
        pos = int(r["host_pos"])
        local = sv[sv["chrom"] == chrom]
        if local.empty:
            continue
        d = (local["pos"] - pos).abs()
        out.at[idx, "n_sv_breakpoints_50kb"] = int((d <= 50000).sum())
        out.at[idx, "n_sv_breakpoints_100kb"] = int((d <= 100000).sum())
        out.at[idx, "n_sv_breakpoints_500kb"] = int((d <= 500000).sum())

        near = local[d <= 100000]
        types = set(near["svtype"].tolist())
        out.at[idx, "has_local_del"] = 1 if "DEL" in types else 0
        out.at[idx, "has_local_dup"] = 1 if "DUP" in types else 0
        out.at[idx, "has_local_inv"] = 1 if "INV" in types else 0
        out.at[idx, "has_local_tra"] = 1 if ("TRA" in types or "BND" in types) else 0
        flags = [out.at[idx, "has_local_del"], out.at[idx, "has_local_dup"], out.at[idx, "has_local_inv"], out.at[idx, "has_local_tra"]]
        out.at[idx, "local_complex_sv_flag"] = 1 if sum(flags) >= 2 else 0

    return out
