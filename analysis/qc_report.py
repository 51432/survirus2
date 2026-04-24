import os

import pandas as pd


def write_qc_tables(events: pd.DataFrame, sample_summary: pd.DataFrame, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)

    qc_event_count = events.groupby("sample_id", as_index=False).size().rename(columns={"size": "n_events"})
    qc_event_count.to_csv(os.path.join(outdir, "qc_event_count_by_sample.tsv"), sep="\t", index=False)

    qc_remap = sample_summary[["sample_id", "remap_confirmation_rate"]]
    qc_remap.to_csv(os.path.join(outdir, "qc_remap_confirmation_rate.tsv"), sep="\t", index=False)

    qc_alt = events.groupby("sample_id", as_index=False)["multi_mapping_risk"].mean().rename(columns={"multi_mapping_risk": "multi_mapping_risk_rate"})
    qc_alt.to_csv(os.path.join(outdir, "qc_alt_mapping_risk.tsv"), sep="\t", index=False)

    qc_hpv = events.groupby("virus_type", as_index=False).size().rename(columns={"size": "event_count"})
    qc_hpv.to_csv(os.path.join(outdir, "qc_hpv_type_distribution.tsv"), sep="\t", index=False)

    qc_host = events.groupby("host_region_class", as_index=False).size().rename(columns={"size": "event_count"})
    qc_host.to_csv(os.path.join(outdir, "qc_host_region_distribution.tsv"), sep="\t", index=False)


def write_summary_report(events: pd.DataFrame, sample_summary: pd.DataFrame, outdir: str) -> None:
    n_samples = sample_summary["sample_id"].nunique()
    pos_samples = int((sample_summary["has_integration"].astype(int) == 1).sum())
    total_events = len(events)
    hc_events = int(pd.to_numeric(events["high_conf_event"], errors="coerce").fillna(0).astype(int).sum())

    hpv = events["virus_type"].value_counts(dropna=False)
    host_region = events["host_region_class"].value_counts(dropna=False)
    coding = int((events["noncoding_flag"].astype(str) == "0").sum())
    noncoding = int((events["noncoding_flag"].astype(str) == "1").sum())

    lines = [
        "# Analysis Summary Report",
        "",
        f"- Total samples: {n_samples}",
        f"- Integration-positive samples: {pos_samples}",
        f"- Total integration events: {total_events}",
        f"- High-confidence events: {hc_events}",
        "",
        "## HPV type distribution",
    ]
    for k, v in hpv.items():
        lines.append(f"- {k}: {v}")

    lines += [
        "",
        "## Host functional distribution",
        f"- coding: {coding}",
        f"- noncoding: {noncoding}",
        "",
        "## Host region class distribution",
    ]
    for k, v in host_region.items():
        lines.append(f"- {k}: {v}")

    lines += [
        "",
        "## Remap confirmation rate",
        f"- mean: {pd.to_numeric(sample_summary['remap_confirmation_rate'], errors='coerce').mean():.3f}",
        "",
        "## Multi-mapping risk",
        f"- mean event risk: {pd.to_numeric(events['multi_mapping_risk'], errors='coerce').mean():.3f}",
    ]

    with open(os.path.join(outdir, "summary_report.md"), "w") as fw:
        fw.write("\n".join(lines) + "\n")
