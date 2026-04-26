#!/usr/bin/env python3
import argparse
import csv
import logging
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


EVENT_RE = re.compile(
    r"^ID=(?P<event_id>\S+)\s+"
    r"(?P<host_chr>[^:]+):(?P<host_strand>[+-])(?P<host_pos>\d+)\s+"
    r"(?P<virus_contig>[^:]+):(?P<virus_strand>[+-])(?P<virus_pos>\d+)"
)
ALT_ID_RE = re.compile(r"ID=(?P<event_id>\S+)")


@dataclass
class Inputs:
    sample_id: str
    main_path: Path
    remap_path: Path
    alt_path: Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Annotate SurVirus integration events with host/virus features.")
    p.add_argument("--results_root", required=True, help="Directory containing results/<sample_id>/ folders.")
    p.add_argument("--gene_bed", required=True)
    p.add_argument("--tss_bed", required=True)
    p.add_argument("--promoter_bed", required=True)
    p.add_argument("--enhancer_bed", default="", help="Optional.")
    p.add_argument("--hpv_bed", required=True)
    p.add_argument("--sample_manifest", default="", help="Optional TSV with at least sample_id column.")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def setup_logging(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "survirus_annotation.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )


def discover_inputs(results_root: Path, sample_manifest: Optional[Path]) -> List[Inputs]:
    samples: List[str] = []
    if sample_manifest and sample_manifest.exists():
        mdf = pd.read_csv(sample_manifest, sep="\t", dtype=str).fillna("")
        if "sample_id" not in mdf.columns:
            raise ValueError("sample_manifest must include sample_id column")
        samples = [str(x).strip() for x in mdf["sample_id"].tolist() if str(x).strip()]
    else:
        samples = sorted([p.name for p in results_root.iterdir() if p.is_dir()])

    inputs: List[Inputs] = []
    for sid in samples:
        sdir = results_root / sid
        main_path = sdir / "results.t1.txt"
        remap_path = sdir / "results.remapped.t1.txt"
        alt_path = sdir / "results.alternative.txt"
        if not main_path.exists():
            logging.warning("Skip sample %s: missing %s", sid, main_path)
            continue
        inputs.append(Inputs(sample_id=sid, main_path=main_path, remap_path=remap_path, alt_path=alt_path))
    return inputs


def virus_type_from_contig(virus_contig: str) -> str:
    core = virus_contig.split("|")[0]
    m = re.match(r"([A-Za-z]+\d+)", core)
    return m.group(1) if m else core


def _parse_line(line: str) -> Optional[Dict[str, str]]:
    m = EVENT_RE.match(line.strip())
    if not m:
        return None
    d = m.groupdict()
    suffix = line[m.end():].strip()
    for token in suffix.split():
        if "=" not in token:
            continue
        k, v = token.split("=", 1)
        d[k.lower()] = v
    return d


def parse_main(path: Path, sample_id: str) -> pd.DataFrame:
    rows = []
    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            d = _parse_line(line)
            if not d:
                logging.warning("Unparseable line in %s: %s", path, line)
                continue
            rows.append({
                "sample_id": sample_id,
                "event_id": str(d.get("event_id", "")),
                "event_uid": f"{sample_id}|{d.get('event_id', '')}",
                "host_chr": d.get("host_chr", ""),
                "host_strand": d.get("host_strand", "."),
                "host_pos": d.get("host_pos", "0"),
                "virus_contig": d.get("virus_contig", ""),
                "virus_type": virus_type_from_contig(d.get("virus_contig", "")),
                "virus_strand": d.get("virus_strand", "."),
                "virus_pos": d.get("virus_pos", "0"),
                "supporting_pairs": d.get("supporting_pairs", 0),
                "split_reads": d.get("split_reads", 0),
                "host_pbs": d.get("host_pbs", ""),
                "coverage": d.get("coverage", ""),
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for c in ["host_pos", "virus_pos", "supporting_pairs", "split_reads"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    for c in ["host_pbs", "coverage"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def parse_remap(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["event_id", "remap_confirmed", "remap_supporting_pairs", "remap_split_reads"])
    df = parse_main(path, sample_id="")
    if df.empty:
        return pd.DataFrame(columns=["event_id", "remap_confirmed", "remap_supporting_pairs", "remap_split_reads"])
    out = df[["event_id", "supporting_pairs", "split_reads"]].rename(
        columns={"supporting_pairs": "remap_supporting_pairs", "split_reads": "remap_split_reads"}
    )
    out["event_id"] = out["event_id"].astype(str)
    out["remap_confirmed"] = 1
    return out


def parse_alt(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["event_id", "n_alt_total", "n_alt_extra", "multi_mapping_risk"])
    counts: Dict[str, int] = {}
    with path.open() as f:
        for raw in f:
            m = ALT_ID_RE.search(raw)
            if not m:
                continue
            eid = m.group("event_id")
            counts[eid] = counts.get(eid, 0) + 1
    rows = []
    for eid, n in counts.items():
        extra = max(n - 1, 0)
        rows.append({"event_id": eid, "n_alt_total": n, "n_alt_extra": extra, "multi_mapping_risk": int(extra > 0)})
    return pd.DataFrame(rows)


def merge_sample(main_df: pd.DataFrame, remap_df: pd.DataFrame, alt_df: pd.DataFrame) -> pd.DataFrame:
    df = main_df.copy()
    df["event_id"] = df["event_id"].astype(str)
    if not remap_df.empty:
        df = df.merge(remap_df, on="event_id", how="left")
    if not alt_df.empty:
        df = df.merge(alt_df, on="event_id", how="left")

    defaults = {
        "remap_confirmed": 0,
        "remap_supporting_pairs": 0,
        "remap_split_reads": 0,
        "n_alt_total": 0,
        "n_alt_extra": 0,
        "multi_mapping_risk": 0,
    }
    for c, v in defaults.items():
        if c not in df.columns:
            df[c] = v
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(v).astype(int)

    df["high_conf_event"] = (
        ((df["split_reads"] >= 1) | (df["supporting_pairs"] >= 5))
        & (df["remap_confirmed"] == 1)
        & (df["n_alt_extra"] == 0)
    ).astype(int)
    return df


def run_cmd(args: List[str], out_path: Path) -> None:
    with out_path.open("w") as fout:
        subprocess.run(args, check=True, stdout=fout, text=True)


def annotate_with_bedtools(df: pd.DataFrame, outdir: Path, gene_bed: Path, tss_bed: Path, promoter_bed: Path, enhancer_bed: Optional[Path], hpv_bed: Path) -> pd.DataFrame:
    host_tmp = outdir / "_host_points.noheader.bed"
    virus_tmp = outdir / "_virus_points.noheader.bed"

    host_df = df[["host_chr", "host_pos", "event_uid", "host_strand"]].copy()
    host_df["start"] = (host_df["host_pos"] - 1).clip(lower=0)
    host_df["end"] = host_df["host_pos"]
    host_df = host_df[["host_chr", "start", "end", "event_uid", "host_strand"]]
    host_df.to_csv(host_tmp, sep="\t", index=False, header=False)

    virus_df = df[["virus_contig", "virus_pos", "event_uid", "virus_strand"]].copy()
    virus_df["start"] = (virus_df["virus_pos"] - 1).clip(lower=0)
    virus_df["end"] = virus_df["virus_pos"]
    virus_df = virus_df[["virus_contig", "start", "end", "event_uid", "virus_strand"]]
    virus_df.to_csv(virus_tmp, sep="\t", index=False, header=False)

    # Host BED outputs with headers for user-facing intermediate tables.
    host_df.rename(columns={"start": "host_start", "end": "host_end"}).to_csv(outdir / "survirus_host_breakpoints.bed", sep="\t", index=False)
    virus_df.rename(columns={"virus_contig": "virus_chr", "start": "virus_start", "end": "virus_end"}).to_csv(outdir / "survirus_virus_breakpoints.bed", sep="\t", index=False)

    gene_i = outdir / "_host_gene.intersect.tsv"
    tss_c = outdir / "_host_tss.closest.tsv"
    prom_c = outdir / "_host_promoter.count.tsv"
    enh_c = outdir / "_host_enhancer.count.tsv"
    hpv_i = outdir / "_virus_hpv.intersect.tsv"

    run_cmd(["bedtools", "intersect", "-a", str(host_tmp), "-b", str(gene_bed), "-wa", "-wb"], gene_i)
    run_cmd(["bedtools", "closest", "-a", str(host_tmp), "-b", str(tss_bed), "-d", "-t", "first"], tss_c)
    run_cmd(["bedtools", "intersect", "-a", str(host_tmp), "-b", str(promoter_bed), "-c"], prom_c)
    if enhancer_bed and enhancer_bed.exists():
        run_cmd(["bedtools", "intersect", "-a", str(host_tmp), "-b", str(enhancer_bed), "-c"], enh_c)
    else:
        with enh_c.open("w") as f:
            for row in host_df.itertuples(index=False):
                f.write("\t".join(map(str, row)) + "\t0\n")
    run_cmd(["bedtools", "intersect", "-a", str(virus_tmp), "-b", str(hpv_bed), "-wa", "-wb"], hpv_i)

    gene_hits: Dict[str, List[List[str]]] = {}
    with gene_i.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            uid = parts[3]
            gene_hits.setdefault(uid, []).append(parts[5:])

    nearest_gene: Dict[str, str] = {}
    dist_tss: Dict[str, int] = {}
    with tss_c.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            uid = parts[3]
            dist_tss[uid] = int(float(parts[-1]))
            nearest_gene[uid] = parts[8] if len(parts) > 8 else (parts[7] if len(parts) > 7 else "NA")

    prom_flag: Dict[str, int] = {}
    with prom_c.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            prom_flag[parts[3]] = int(parts[-1])

    enh_flag: Dict[str, int] = {}
    with enh_c.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            enh_flag[parts[3]] = int(parts[-1])

    virus_map: Dict[str, Tuple[str, str]] = {}
    with hpv_i.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            uid = parts[3]
            if uid in virus_map:
                continue
            region = parts[8] if len(parts) > 8 else "NA"
            seg = parts[9] if len(parts) > 9 else region
            virus_map[uid] = (region, seg)

    host_gene_list = []
    nearest_gene_list = []
    dist_list = []
    cls_list = []
    detail_list = []
    prom_list = []
    enh_list = []
    noncoding_list = []
    vregion_list = []
    vseg_list = []

    for r in df.itertuples(index=False):
        uid = r.event_uid
        hits = gene_hits.get(uid, [])
        hit_text = " ".join(" ".join(h).lower() for h in hits)

        host_gene = "NA"
        if hits:
            first = hits[0]
            host_gene = first[3] if len(first) > 3 else first[0]

        pflag = int(prom_flag.get(uid, 0) > 0)
        eflag = int(enh_flag.get(uid, 0) > 0)
        is_exonic = any("exon" in " ".join(h).lower() or "cds" in " ".join(h).lower() for h in hits)
        is_coding_exon = any("cds" in " ".join(h).lower() for h in hits)
        is_intronic = bool(hits) and not is_exonic

        if pflag:
            cls = "promoter_proximal"
            detail = "upstream_2kb"
        elif is_exonic:
            cls = "exonic"
            if "5utr" in hit_text or "five_prime_utr" in hit_text:
                detail = "5UTR"
            elif "3utr" in hit_text or "three_prime_utr" in hit_text:
                detail = "3UTR"
            elif "cds" in hit_text:
                detail = "CDS"
            else:
                detail = "exon"
        elif is_intronic:
            cls = "intronic"
            detail = "intron"
        elif eflag:
            cls = "enhancer_proximal"
            detail = "enhancer"
        else:
            cls = "intergenic"
            detail = "distal_intergenic"

        noncoding = 0 if is_coding_exon else 1
        v_region, v_seg = virus_map.get(uid, ("NA", "NA"))

        host_gene_list.append(host_gene)
        nearest_gene_list.append(nearest_gene.get(uid, "NA"))
        dist_list.append(dist_tss.get(uid, -1))
        cls_list.append(cls)
        detail_list.append(detail)
        prom_list.append(pflag)
        enh_list.append(eflag)
        noncoding_list.append(noncoding)
        vregion_list.append(v_region)
        vseg_list.append(v_seg)

    out = df.copy()
    out["host_gene"] = host_gene_list
    out["nearest_gene"] = nearest_gene_list
    out["distance_to_tss"] = dist_list
    out["host_region_class"] = cls_list
    out["host_region_detail"] = detail_list
    out["promoter_flag"] = prom_list
    out["enhancer_flag"] = enh_list
    out["noncoding_flag"] = noncoding_list
    out["virus_region"] = vregion_list
    out["virus_gene_segment"] = vseg_list
    return out


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    setup_logging(outdir)

    inputs = discover_inputs(Path(args.results_root), Path(args.sample_manifest) if args.sample_manifest else None)
    if not inputs:
        raise RuntimeError("No valid sample results found.")

    tables = []
    for inp in inputs:
        logging.info("Processing sample %s", inp.sample_id)
        main_df = parse_main(inp.main_path, inp.sample_id)
        if main_df.empty:
            logging.warning("Sample %s has empty main table.", inp.sample_id)
            continue
        remap_df = parse_remap(inp.remap_path)
        alt_df = parse_alt(inp.alt_path)
        merged = merge_sample(main_df, remap_df, alt_df)
        tables.append(merged)

    if not tables:
        raise RuntimeError("No events parsed from all samples.")

    events = pd.concat(tables, ignore_index=True)
    annotated = annotate_with_bedtools(
        events,
        outdir=outdir,
        gene_bed=Path(args.gene_bed),
        tss_bed=Path(args.tss_bed),
        promoter_bed=Path(args.promoter_bed),
        enhancer_bed=Path(args.enhancer_bed) if args.enhancer_bed else None,
        hpv_bed=Path(args.hpv_bed),
    )

    cols = [
        "sample_id", "event_id", "event_uid", "host_chr", "host_pos", "host_strand",
        "virus_contig", "virus_type", "virus_pos", "virus_strand", "supporting_pairs", "split_reads",
        "host_pbs", "coverage", "remap_confirmed", "n_alt_total", "n_alt_extra", "multi_mapping_risk",
        "high_conf_event", "host_gene", "nearest_gene", "distance_to_tss", "host_region_class",
        "host_region_detail", "promoter_flag", "enhancer_flag", "noncoding_flag", "virus_region",
        "virus_gene_segment",
    ]
    annotated[cols].to_csv(outdir / "integration_event_annotation.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    logging.info("Wrote %s", outdir / "integration_event_annotation.tsv")


if __name__ == "__main__":
    main()
