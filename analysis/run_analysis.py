#!/usr/bin/env python3
import argparse
import os

import pandas as pd

from .annotate_events import annotate_events
from .integrate_ascat import annotate_local_cnv
from .integrate_rna import annotate_local_expression
from .integrate_sv import annotate_local_sv
from .parse_survirus import merge_survirus_tables, parse_alt_results, parse_main_results, parse_remapped
from .qc_report import write_qc_tables, write_summary_report
from .utils import ensure_columns, read_manifest, setup_logging


def parse_args():
    p = argparse.ArgumentParser(description="Independent downstream integration analysis step.")
    p.add_argument("--manifest", required=True)
    p.add_argument("--gene_gtf", default="")
    p.add_argument("--gene_body_bed", default="")
    p.add_argument("--exon_bed", default="")
    p.add_argument("--intron_bed", default="")
    p.add_argument("--cds_bed", default="")
    p.add_argument("--tss_bed", default="")
    p.add_argument("--promoter_bed", default="")
    p.add_argument("--enhancer_bed", default="")
    p.add_argument("--hpv_bed", default="")
    p.add_argument("--outdir", required=True)

    p.add_argument("--hc_min_split_reads", type=int, default=1)
    p.add_argument("--hc_min_supporting_pairs", type=int, default=5)
    p.add_argument("--hc_require_remap", type=int, default=1)
    p.add_argument("--hc_require_no_alt_extra", type=int, default=1)

    p.add_argument("--cis_distance_bp", type=int, default=100000)
    p.add_argument("--cis_z_cutoff", type=float, default=2.0)
    return p.parse_args()


def write_bed_files(events: pd.DataFrame, outdir: str):
    host_bed = events[["host_chr", "host_pos", "event_uid", "sample_id"]].copy()
    host_bed["start"] = host_bed["host_pos"] - 1
    host_bed["end"] = host_bed["host_pos"]
    host_bed = host_bed[["host_chr", "start", "end", "event_uid", "sample_id"]]
    host_bed.to_csv(os.path.join(outdir, "survirus_host_breakpoints.bed"), sep="\t", index=False, header=False)

    virus_bed = events[["virus_contig", "virus_pos", "event_uid", "sample_id"]].copy()
    virus_bed["start"] = virus_bed["virus_pos"] - 1
    virus_bed["end"] = virus_bed["virus_pos"]
    virus_bed = virus_bed[["virus_contig", "start", "end", "event_uid", "sample_id"]]
    virus_bed.to_csv(os.path.join(outdir, "survirus_virus_breakpoints.bed"), sep="\t", index=False, header=False)

    for w in [50000, 100000, 500000]:
        bed = events[["host_chr", "host_pos", "event_uid", "sample_id"]].copy()
        bed["start"] = (bed["host_pos"] - w).clip(lower=0)
        bed["end"] = bed["host_pos"] + w
        bed = bed[["host_chr", "start", "end", "event_uid", "sample_id"]]
        bed.to_csv(os.path.join(outdir, f"survirus_windows_{int(w/1000)}kb.bed"), sep="\t", index=False, header=False)


def compute_high_conf(df: pd.DataFrame, args) -> pd.Series:
    cond = ((df["split_reads"] >= args.hc_min_split_reads) | (df["supporting_pairs"] >= args.hc_min_supporting_pairs))
    if args.hc_require_remap:
        cond = cond & (df["remap_confirmed"] == 1)
    if args.hc_require_no_alt_extra:
        cond = cond & (df["n_alt_extra"] == 0)
    return cond.astype(int)


def build_sample_summary(events: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    host_regions = ["coding_exon", "noncoding_exon", "intronic", "intergenic", "promoter_proximal", "enhancer_proximal"]
    virus_regions = ["E1", "E2", "E5", "E6", "E7", "L1", "L2", "URR"]

    rows = []
    for _, m in manifest.iterrows():
        sid = m["sample_id"]
        sdf = events[events["sample_id"] == sid]
        if sdf.empty:
            base = {
                "sample_id": sid,
                "subtype": m.get("subtype", ""),
                "has_integration": 0,
                "n_events_raw": 0,
                "n_events_remap_confirmed": 0,
                "n_unique_host_breakpoints": 0,
                "n_unique_virus_breakpoints": 0,
                "major_hpv_type": "NA",
                "sum_supporting_pairs": 0,
                "sum_split_reads": 0,
                "median_supporting_pairs": pd.NA,
                "median_split_reads": pd.NA,
                "mean_host_pbs": pd.NA,
                "mean_coverage": pd.NA,
                "n_events_with_alt_hits": 0,
                "high_conf_event_count": 0,
                "remap_confirmation_rate": pd.NA,
                "n_exonic": 0,
                "n_noncoding_total": 0,
                "dominant_host_region": "NA",
                "dominant_virus_region": "NA",
                "virus_region_combo": "NA",
                "host_region_combo": "NA",
            }
            for hr in host_regions:
                base[f"n_{hr}"] = 0
            for vr in virus_regions:
                base[f"n_virus_{vr}"] = 0
            rows.append(base)
            continue

        host_vc = sdf["host_region_class"].value_counts()
        virus_vc = sdf["virus_region"].value_counts()
        base = {
            "sample_id": sid,
            "subtype": m.get("subtype", ""),
            "has_integration": 1,
            "n_events_raw": len(sdf),
            "n_events_remap_confirmed": int(sdf["remap_confirmed"].sum()),
            "n_unique_host_breakpoints": sdf[["host_chr", "host_pos"]].drop_duplicates().shape[0],
            "n_unique_virus_breakpoints": sdf[["virus_contig", "virus_pos"]].drop_duplicates().shape[0],
            "major_hpv_type": sdf["virus_type"].mode().iloc[0] if not sdf["virus_type"].mode().empty else "NA",
            "sum_supporting_pairs": int(sdf["supporting_pairs"].sum()),
            "sum_split_reads": int(sdf["split_reads"].sum()),
            "median_supporting_pairs": sdf["supporting_pairs"].median(),
            "median_split_reads": sdf["split_reads"].median(),
            "mean_host_pbs": sdf["host_pbs"].mean(),
            "mean_coverage": sdf["coverage"].mean(),
            "n_events_with_alt_hits": int((sdf["n_alt_extra"] > 0).sum()),
            "high_conf_event_count": int(sdf["high_conf_event"].sum()),
            "remap_confirmation_rate": float(sdf["remap_confirmed"].mean()) if len(sdf) else pd.NA,
            "n_exonic": int((sdf["exonic_flag"] == 1).sum()),
            "n_noncoding_total": int((sdf["noncoding_flag"] == 1).sum()),
            "dominant_host_region": host_vc.index[0] if not host_vc.empty else "NA",
            "dominant_virus_region": virus_vc.index[0] if not virus_vc.empty else "NA",
            "virus_region_combo": ";".join(sorted([x for x in sdf["virus_region"].dropna().astype(str).unique() if x and x != "NA"])) or "NA",
            "host_region_combo": ";".join(sorted([x for x in sdf["host_region_class"].dropna().astype(str).unique() if x and x != "NA"])) or "NA",
        }
        for hr in host_regions:
            base[f"n_{hr}"] = int((sdf["host_region_class"] == hr).sum())
        for vr in virus_regions:
            base[f"n_virus_{vr}"] = int((sdf["virus_region"] == vr).sum())
        rows.append(base)
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args.outdir)

    manifest = read_manifest(args.manifest)
    all_events = []
    multiomics = []

    for _, row in manifest.iterrows():
        sid = row["sample_id"]
        subtype = row.get("subtype", "")

        main_df = parse_main_results(row["survirus_results_t1"], sid, subtype)
        if main_df.empty:
            continue

        remap_df = parse_remapped(row.get("survirus_results_remapped_t1", ""))
        alt_df = parse_alt_results(row.get("survirus_results_alternative", ""))
        merged = merge_survirus_tables(main_df, remap_df, alt_df)
        merged["high_conf_event"] = compute_high_conf(merged, args)

        annotated = annotate_events(
            merged,
            gene_gtf=args.gene_gtf,
            gene_body_bed=args.gene_body_bed,
            exon_bed=args.exon_bed,
            intron_bed=args.intron_bed,
            cds_bed=args.cds_bed,
            tss_bed=args.tss_bed,
            promoter_bed=args.promoter_bed,
            enhancer_bed=args.enhancer_bed,
            hpv_bed=args.hpv_bed,
        )
        all_events.append(annotated)

        cnv = annotate_local_cnv(annotated, row)
        sv = annotate_local_sv(cnv, row.get("manta_sv", ""), row.get("tiddit_sv", ""))
        mo = annotate_local_expression(sv, row.get("rna_expr", ""), args.cis_distance_bp, args.cis_z_cutoff)
        multiomics.append(mo)

    if all_events:
        events = pd.concat(all_events, ignore_index=True)
    else:
        events = pd.DataFrame(columns=["sample_id", "subtype", "event_uid", "event_id", "host_chr", "host_strand", "host_pos", "virus_contig", "virus_type", "virus_strand", "virus_pos", "supporting_pairs", "split_reads", "host_pbs", "coverage", "remap_confirmed", "n_alt_total", "n_alt_extra", "multi_mapping_risk", "high_conf_event"])

    events = ensure_columns(events, [
        "host_gene", "nearest_transcript", "nearest_gene", "distance_to_tss", "host_region_class", "host_region_detail",
        "promoter_flag", "enhancer_flag", "coding_flag", "exonic_flag", "intronic_flag", "noncoding_flag",
        "virus_region", "virus_gene_segment"
    ])

    sample_summary = build_sample_summary(events, manifest)
    sample_summary.to_csv(os.path.join(args.outdir, "sample_level_integration_summary.tsv"), sep="\t", index=False)

    event_cols = [
        "sample_id", "subtype", "event_uid", "event_id", "host_chr", "host_strand", "host_pos", "virus_contig", "virus_type",
        "virus_strand", "virus_pos", "supporting_pairs", "split_reads", "host_pbs", "coverage",
        "remap_confirmed", "n_alt_total", "n_alt_extra", "multi_mapping_risk", "high_conf_event",
        "host_gene", "nearest_transcript", "nearest_gene", "distance_to_tss",
        "host_region_class", "host_region_detail", "promoter_flag", "enhancer_flag",
        "coding_flag", "exonic_flag", "intronic_flag", "noncoding_flag",
        "virus_region", "virus_gene_segment",
    ]
    events[event_cols].to_csv(os.path.join(args.outdir, "integration_event_annotation.tsv"), sep="\t", index=False)

    if multiomics:
        mo_df = pd.concat(multiomics, ignore_index=True)
    else:
        mo_df = events.copy()

    mo_cols = [
        "sample_id", "subtype", "event_uid", "host_chr", "host_pos", "virus_type", "supporting_pairs", "split_reads", "remap_confirmed", "n_alt_extra", "nearest_gene", "distance_to_tss", "host_region_class", "virus_gene_segment",
        "segment_start", "segment_end", "local_total_cn", "local_major_cn", "local_minor_cn", "sample_ploidy", "sample_wgd_status",
        "n_sv_breakpoints_50kb", "n_sv_breakpoints_100kb", "n_sv_breakpoints_500kb", "has_local_del", "has_local_dup", "has_local_inv", "has_local_tra", "local_complex_sv_flag",
        "nearest_gene_expr", "nearest_gene_zscore", "cis_effect_flag",
    ]
    mo_df = ensure_columns(mo_df, mo_cols)
    mo_df[mo_cols].to_csv(os.path.join(args.outdir, "integration_local_multiomics.tsv"), sep="\t", index=False)

    if not events.empty:
        write_bed_files(events, args.outdir)
    else:
        for n in ["survirus_host_breakpoints.bed", "survirus_virus_breakpoints.bed", "survirus_windows_50kb.bed", "survirus_windows_100kb.bed", "survirus_windows_500kb.bed"]:
            open(os.path.join(args.outdir, n), "w").close()

    write_qc_tables(events, sample_summary, args.outdir)
    write_summary_report(events, sample_summary, args.outdir)


if __name__ == "__main__":
    main()
