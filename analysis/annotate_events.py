from typing import Optional

import pandas as pd

from .utils import build_interval_index, nearest_interval, parse_bed, query_overlaps, read_gene_gtf


def annotate_events(
    events: pd.DataFrame,
    gene_gtf: Optional[str],
    tss_bed: Optional[str],
    promoter_bed: Optional[str],
    enhancer_bed: Optional[str],
    hpv_bed: Optional[str],
) -> pd.DataFrame:
    df = events.copy()

    genes = read_gene_gtf(gene_gtf)
    tss = parse_bed(tss_bed)
    promoter = parse_bed(promoter_bed)
    enhancer = parse_bed(enhancer_bed)
    hpv = parse_bed(hpv_bed)

    gene_idx = build_interval_index(genes)
    tss_idx = build_interval_index(tss)
    promoter_idx = build_interval_index(promoter)
    enhancer_idx = build_interval_index(enhancer)
    hpv_idx = build_interval_index(hpv)

    nearest_gene = []
    dist_tss = []
    host_cls = []
    host_detail = []
    promoter_flag = []
    enhancer_flag = []
    noncoding_flag = []
    virus_region = []
    virus_gene_seg = []

    for _, r in df.iterrows():
        chrom = str(r["host_chr"])
        pos = int(r["host_pos"])

        g_meta, _ = nearest_interval(gene_idx, chrom, pos)
        nearest_gene.append(g_meta.get("gene", "NA") if g_meta else "NA")

        t_meta, t_dist = nearest_interval(tss_idx, chrom, pos)
        dist_tss.append(t_dist if t_dist is not None else pd.NA)

        p_hits = query_overlaps(promoter_idx, chrom, pos)
        e_hits = query_overlaps(enhancer_idx, chrom, pos)
        g_hits = query_overlaps(gene_idx, chrom, pos)

        pflag = 1 if p_hits else 0
        eflag = 1 if e_hits else 0
        promoter_flag.append(pflag)
        enhancer_flag.append(eflag)

        if pflag:
            cls = "promoter_proximal"
        elif g_hits:
            cls = "intronic"
        elif eflag:
            cls = "enhancer_proximal"
        else:
            cls = "intergenic"

        # exon-level detail unavailable from gene-only GTF; keep exonic branch for extensibility
        if cls == "intronic" and g_meta:
            detail = f"gene:{g_meta.get('gene', 'NA')}"
        elif cls == "promoter_proximal" and p_hits:
            detail = f"promoter:{p_hits[0].get('name', '.') }"
        elif cls == "enhancer_proximal" and e_hits:
            detail = f"enhancer:{e_hits[0].get('name', '.') }"
        else:
            detail = "intergenic"

        host_cls.append(cls)
        host_detail.append(detail)
        noncoding_flag.append(0 if cls in {"exonic", "intronic", "promoter_proximal"} else 1)

        vchrom = str(r["virus_contig"])
        vpos = int(r["virus_pos"])
        v_hits = query_overlaps(hpv_idx, vchrom, vpos)
        if v_hits:
            virus_region.append(v_hits[0].get("name", "NA"))
            virus_gene_seg.append(v_hits[0].get("name", "NA"))
        else:
            virus_region.append("NA")
            virus_gene_seg.append("NA")

    df["nearest_gene"] = nearest_gene
    df["distance_to_tss"] = dist_tss
    df["host_region_class"] = host_cls
    df["host_region_detail"] = host_detail
    df["promoter_flag"] = promoter_flag
    df["enhancer_flag"] = enhancer_flag
    df["noncoding_flag"] = noncoding_flag
    df["virus_region"] = virus_region
    df["virus_gene_segment"] = virus_gene_seg
    return df
