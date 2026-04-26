# SurVirus downstream event annotation

该模块不修改 SurVirus caller 本体，仅基于 `results/<sample_id>/` 现有输出新增事件级注释。

## 输入

每个样本目录需要：

- `results.t1.txt`
- `results.remapped.t1.txt`（可缺失）
- `results.alternative.txt`（可缺失）

以及外部注释 BED：

- `--gene_bed`
- `--tss_bed`
- `--promoter_bed`
- `--enhancer_bed`（可选）
- `--hpv_bed`

## 运行

```bash
python3 -m analysis.run_analysis \
  --manifest test_manifest.tsv \
  --gene_body_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/gene_body.sorted.bed \
  --exon_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/exons.sorted.bed \
  --intron_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/introns.sorted.bed \
  --cds_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/cds.sorted.bed \
  --tss_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/tss.sorted.bed \
  --promoter_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/promoter_2kb.sorted.bed \
  --enhancer_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/enhancers.cCRE_ELS.hg38.sorted.bed \
  --hpv_bed /data/person/wup/public/liusy_files/reference_genomes/virus/annotation/gff/hpv_annotation_simple.survirus.bed \
  --outdir analysis_test_out_full

```

`sample_manifest` 可选（仅需 `sample_id` 列）；未提供时自动扫描 `results_root` 下所有子目录。

## 输出

- `integration_event_annotation.tsv`（事件级整合注释主表，含 `nearest_transcript` / `nearest_gene` / `distance_to_tss`）
- `survirus_host_breakpoints.bed`
- `survirus_virus_breakpoints.bed`
- `survirus_annotation.log`

### integration_event_annotation.tsv 主要列说明

| 列名 | 含义 |
|---|---|
| `sample_id` | 样本 ID。 |
| `event_id` | SurVirus 在样本内的事件编号（来自 `ID=`）。 |
| `event_uid` | 全局唯一事件 ID，格式 `sample_id|event_id`。 |
| `host_chr` / `host_pos` / `host_strand` | 宿主断点染色体、坐标与方向。 |
| `virus_contig` / `virus_pos` / `virus_strand` | 病毒断点 contig、坐标与方向。 |
| `virus_type` | 从 `virus_contig` 提取的病毒类型（如 `HPV39REF` → `HPV39`）。 |
| `supporting_pairs` | 原始结果中的配对 reads 支持数。 |
| `split_reads` | 原始结果中的 split reads 支持数。 |
| `host_pbs` | 宿主断点 PBS 分值（原始输出字段）。 |
| `coverage` | 断点覆盖度（原始输出字段）。 |
| `remap_confirmed` | 是否在 `results.remapped.t1.txt` 中被再次确认（1/0）。 |
| `n_alt_total` | `results.alternative.txt` 中该事件出现的总条数。 |
| `n_alt_extra` | 额外替代命中数，`max(n_alt_total - 1, 0)`。 |
| `multi_mapping_risk` | 多重比对风险标记；当 `n_alt_extra > 0` 时为 1。 |
| `high_conf_event` | 高置信事件标记（当前规则：`split_reads>=1` 或 `supporting_pairs>=5`，且 `remap_confirmed=1` 且 `n_alt_extra=0`）。 |
| `host_gene` | 若断点落在基因体范围内，填写该基因（否则 `NA`）。 |
| `nearest_transcript` | 与断点最近 TSS 对应 transcript（`bedtools closest`）。 |
| `nearest_gene` | 与断点最近 TSS 对应基因名（`bedtools closest`）。 |
| `distance_to_tss` | 到最近 TSS 的距离（bp，`bedtools closest -d` 输出最后一列）。 |
| `promoter_flag` | 是否与 promoter BED 重叠（1/0）。 |
| `enhancer_flag` | 是否与 enhancer BED 重叠（1/0，未提供 enhancer 时默认 0）。 |
| `host_region_class` | 宿主区域粗分类：`promoter_proximal` / `exonic` / `intronic` / `enhancer_proximal` / `intergenic`。 |
| `host_region_detail` | 更细分类细节，如 `upstream_2kb`、`CDS`、`5UTR`、`intron`、`distal_intergenic`。 |
| `noncoding_flag` | 非编码标记（当前实现中，仅编码外显子 `CDS` 记为 0，其余为 1）。 |
| `virus_region` | 病毒断点在 `--hpv_bed` 上的区域注释。 |
| `virus_gene_segment` | 病毒断点对应的基因/片段注释（若 BED 有对应列则取该列）。 |

## 依赖

- `python >= 3.8`
- `pandas`
- `bedtools`（用于 intersect / closest 注释）

## ASCAT purity/ploidy 回填到sample_subtype_metadata.tsv
### 批量提取 ASCAT purity / ploidy / WGD / GI
#你可以直接新建脚本：

```bash
vim extract_ascat_metrics.py
```
### 写入：
```bash
#!/usr/bin/env python3
import re
from pathlib import Path
import pandas as pd

ascat_root = Path("test/result/variant_calling/ascat")
rows = []

for metrics_file in ascat_root.glob("*/*.metrics.txt"):
    pair_id = metrics_file.parent.name

    # 从 TSDX008_vs_NSDX008 提取 tumor sample = TSDX008
    if "_vs_" in pair_id:
        sample_id = pair_id.split("_vs_")[0]
    else:
        sample_id = pair_id

    try:
        df = pd.read_csv(metrics_file, sep=r"\s+", engine="python")
    except Exception as e:
        print(f"[WARN] failed to read {metrics_file}: {e}")
        continue

    if df.shape[0] < 1:
        continue

    r = df.iloc[0].to_dict()

    rows.append({
        "sample_id": sample_id,
        "ascat_pair_id": pair_id,
        "tumor_purity": r.get("purity", pd.NA),
        "ploidy": r.get("ploidy", pd.NA),
        "goodness_of_fit": r.get("goodness_of_fit", pd.NA),
        "WGD": r.get("WGD", pd.NA),
        "GI": r.get("GI", pd.NA),
        "LOH": r.get("LOH", pd.NA),
        "n_segs": r.get("n_segs", pd.NA),
        "homdel_segs": r.get("homdel_segs", pd.NA),
        "homdel_fraction": r.get("homdel_fraction", pd.NA),
        "mode_minA": r.get("mode_minA", pd.NA),
        "mode_majA": r.get("mode_majA", pd.NA),
    })

out = pd.DataFrame(rows)
out = out.sort_values("sample_id")
out.to_csv("ascat_sample_metrics.tsv", sep="\t", index=False)

print(f"Wrote ascat_sample_metrics.tsv with {out.shape[0]} samples")
```

### 运行：
```bash
python extract_ascat_metrics.py
```
得到：ascat_sample_metrics.tsv
## 把 ASCAT 结果合并进你的 metadata
假设你的 metadata 是：sample_subtype_metadata.tsv

然后运行：
```bash
python - <<'PY'
import pandas as pd

meta = pd.read_csv("sample_subtype_metadata.tsv", sep="\t")
ascat = pd.read_csv("ascat_sample_metrics.tsv", sep="\t")

# 如果 metadata 里已有空的 tumor_purity/ploidy，先删掉
meta = meta.drop(columns=["tumor_purity", "ploidy"], errors="ignore")

out = meta.merge(ascat, on="sample_id", how="left")

out.to_csv("sample_subtype_metadata.with_ascat.tsv", sep="\t", index=False)
print(out[["sample_id","subtype","tumor_purity","ploidy","WGD","GI"]].head())
PY
```




