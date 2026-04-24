import os
from typing import Optional

import pandas as pd


def annotate_local_expression(
    events: pd.DataFrame,
    rna_expr_path: Optional[str],
    cis_distance_bp: int = 100000,
    cis_z_cutoff: float = 2.0,
) -> pd.DataFrame:
    out = events.copy()
    out["nearest_gene_expr"] = pd.NA
    out["nearest_gene_zscore"] = pd.NA
    out["cis_effect_flag"] = pd.NA

    if not rna_expr_path or not os.path.exists(rna_expr_path):
        return out

    expr = pd.read_csv(rna_expr_path, sep="\t", dtype=str)
    if expr.empty or expr.shape[1] < 2:
        return out

    gene_col = expr.columns[0]
    expr = expr.rename(columns={gene_col: "gene"})
    for c in expr.columns[1:]:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")

    sample_cols = [c for c in expr.columns if c != "gene"]
    if not sample_cols:
        return out

    vals = expr[sample_cols]
    means = vals.mean(axis=1)
    stds = vals.std(axis=1).replace(0, pd.NA)
    expr["_mean"] = means
    expr["_std"] = stds

    emap = expr.set_index("gene")
    for idx, r in out.iterrows():
        sid = r["sample_id"]
        gene = r.get("nearest_gene", "NA")
        if sid not in emap.columns or gene not in emap.index:
            continue
        val = emap.at[gene, sid]
        mu = emap.at[gene, "_mean"]
        sd = emap.at[gene, "_std"]
        z = pd.NA if pd.isna(val) or pd.isna(mu) or pd.isna(sd) else (val - mu) / sd
        out.at[idx, "nearest_gene_expr"] = val
        out.at[idx, "nearest_gene_zscore"] = z

        dist = r.get("distance_to_tss", pd.NA)
        if pd.isna(z) or pd.isna(dist):
            out.at[idx, "cis_effect_flag"] = pd.NA
        else:
            out.at[idx, "cis_effect_flag"] = 1 if (float(dist) < cis_distance_bp and float(z) > cis_z_cutoff) else 0

    return out
