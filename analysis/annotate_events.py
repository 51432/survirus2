from typing import Optional

import pandas as pd

from .utils import build_interval_index, nearest_interval, parse_bed, query_overlaps, read_gene_gtf


def _get_gene_from_bed_meta(meta: dict) -> str:
    if not meta:
        return "NA"
    return meta.get("gene") or meta.get("gene_name") or meta.get("score") or meta.get("name") or "NA"


def _get_transcript_from_bed_meta(meta: dict) -> str:
    if not meta:
        return "NA"
    return meta.get("transcript_id") or meta.get("name") or "NA"


def annotate_events(
    events: pd.DataFrame,
    gene_gtf: Optional[str],
    gene_body_bed: Optional[str],
    exon_bed: Optional[str],
    intron_bed: Optional[str],
    cds_bed: Optional[str],
    tss_bed: Optional[str],
    promoter_bed: Optional[str],
    enhancer_bed: Optional[str],
    hpv_bed: Optional[str],
) -> pd.DataFrame:
    df = events.copy()

    gene_body = parse_bed(gene_body_bed) if gene_body_bed else pd.DataFrame()
    if gene_body.empty:
        gene_body = read_gene_gtf(gene_gtf)
    exons = parse_bed(exon_bed)
    introns = parse_bed(intron_bed)
    cds = parse_bed(cds_bed)
    tss = parse_bed(tss_bed)
    promoter = parse_bed(promoter_bed)
    enhancer = parse_bed(enhancer_bed)
    hpv = parse_bed(hpv_bed)

    gene_body_idx = build_interval_index(gene_body)
    exon_idx = build_interval_index(exons)
    intron_idx = build_interval_index(introns)
    cds_idx = build_interval_index(cds)
    tss_idx = build_interval_index(tss)
    promoter_idx = build_interval_index(promoter)
    enhancer_idx = build_interval_index(enhancer)
    hpv_idx = build_interval_index(hpv)

    host_gene = []
    nearest_transcript = []
    nearest_gene = []
    dist_tss = []
    host_cls = []
    host_detail = []
    promoter_flag = []
    enhancer_flag = []
    coding_flag = []
    exonic_flag = []
    intronic_flag = []
    noncoding_flag = []
    virus_region = []
    virus_gene_seg = []

    for _, r in df.iterrows():
        chrom = str(r["host_chr"])
        pos = int(r["host_pos"])

        gene_hits = query_overlaps(gene_body_idx, chrom, pos)
        exon_hits = query_overlaps(exon_idx, chrom, pos)
        intron_hits = query_overlaps(intron_idx, chrom, pos)
        cds_hits = query_overlaps(cds_idx, chrom, pos)
        p_hits = query_overlaps(promoter_idx, chrom, pos)
        e_hits = query_overlaps(enhancer_idx, chrom, pos)

        t_meta, t_dist = nearest_interval(tss_idx, chrom, pos)
        if t_meta:
            nearest_transcript.append(_get_transcript_from_bed_meta(t_meta))
            nearest_gene.append(_get_gene_from_bed_meta(t_meta))
            dist_tss.append(t_dist if t_dist is not None else pd.NA)
        else:
            nearest_transcript.append("NA")
            nearest_gene.append("NA")
            dist_tss.append(pd.NA)

        pflag = int(bool(p_hits))
        eflag = int(bool(e_hits))
        promoter_flag.append(pflag)
        enhancer_flag.append(eflag)

        if pflag:
            cls = "promoter_proximal"
            p_meta = p_hits[0]
            p_gene = _get_gene_from_bed_meta(p_meta)
            p_tx = _get_transcript_from_bed_meta(p_meta)
            detail = f"promoter:{p_gene}|{p_tx}"
        elif cds_hits:
            cls = "coding_exon"
            cds_meta = cds_hits[0]
            detail = f"CDS:{_get_gene_from_bed_meta(cds_meta)}|{cds_meta.get('name', '.')}"
        elif exon_hits:
            cls = "noncoding_exon"
            exon_meta = exon_hits[0]
            detail = f"exon:{_get_gene_from_bed_meta(exon_meta)}|{exon_meta.get('name', '.')}"
        elif intron_hits:
            cls = "intronic"
            intron_meta = intron_hits[0]
            detail = f"intron:{_get_gene_from_bed_meta(intron_meta)}|{intron_meta.get('name', '.')}"
        elif eflag:
            cls = "enhancer_proximal"
            detail = f"enhancer:{e_hits[0].get('name', '.')}"
        else:
            cls = "intergenic"
            detail = "intergenic"

        if cds_hits:
            hgene = _get_gene_from_bed_meta(cds_hits[0])
        elif exon_hits:
            hgene = _get_gene_from_bed_meta(exon_hits[0])
        elif intron_hits:
            hgene = _get_gene_from_bed_meta(intron_hits[0])
        elif gene_hits:
            hgene = _get_gene_from_bed_meta(gene_hits[0])
        else:
            hgene = "NA"

        coding = int(bool(cds_hits))
        exonic = int(bool(exon_hits))
        intronic = int(bool(intron_hits))

        host_gene.append(hgene)
        host_cls.append(cls)
        host_detail.append(detail)
        coding_flag.append(coding)
        exonic_flag.append(exonic)
        intronic_flag.append(intronic)
        noncoding_flag.append(0 if cls == "coding_exon" else 1)

        vchrom = str(r["virus_contig"])
        vpos = int(r["virus_pos"])
        v_hits = query_overlaps(hpv_idx, vchrom, vpos)
        if v_hits:
            virus_region.append(v_hits[0].get("name", "NA"))
            virus_gene_seg.append(v_hits[0].get("name", "NA"))
        else:
            virus_region.append("NA")
            virus_gene_seg.append("NA")

    df["host_gene"] = host_gene
    df["nearest_transcript"] = nearest_transcript
    df["nearest_gene"] = nearest_gene
    df["distance_to_tss"] = dist_tss
    df["host_region_class"] = host_cls
    df["host_region_detail"] = host_detail
    df["promoter_flag"] = promoter_flag
    df["enhancer_flag"] = enhancer_flag
    df["coding_flag"] = coding_flag
    df["exonic_flag"] = exonic_flag
    df["intronic_flag"] = intronic_flag
    df["noncoding_flag"] = noncoding_flag
    df["virus_region"] = virus_region
    df["virus_gene_segment"] = virus_gene_seg
    return df
