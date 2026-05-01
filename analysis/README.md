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

## 样本输入格式如下：test_manifest.tsv

| sample_id | survirus_results_t1 | survirus_results_remapped_t1 | survirus_results_alternative | subtype |
| --- | --- | --- | --- | --- |
| TSDX001 | /data/person/wup/liusy/wgs/results/integration/TSDX001/results.t1.txt | /data/person/wup/liusy/wgs/results/integration/TSDX001/results.remapped.t1.txt | /data/person/wup/liusy/wgs/results/integration/TSDX001/results.alternative.txt | NP |
| TSDX002 | /data/person/wup/liusy/wgs/results/integration/TSDX002/results.t1.txt | /data/person/wup/liusy/wgs/results/integration/TSDX002/results.remapped.t1.txt | /data/person/wup/liusy/wgs/results/integration/TSDX002/results.alternative.txt | ES |

## 运行

```bash
cd suevirus
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

nohup python3 -m analysis.run_analysis --manifest test_manifest.tsv --gene_body_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/gene_body.sorted.bed --exon_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/exons.sorted.bed --intron_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/introns.sorted.bed --cds_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/cds.sorted.bed --tss_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/tss.sorted.bed --promoter_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/promoter_2kb.sorted.bed --enhancer_bed /data/person/wup/public/liusy_files/reference_genomes/hg38/reference/host_annotation_beds/enhancers.cCRE_ELS.hg38.sorted.bed --hpv_bed /data/person/wup/public/liusy_files/reference_genomes/virus/annotation/gff/hpv_annotation_simple.survirus.bed --outdir analysis_test_out_full > run_analysis.log 2>&1 &

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
### 得到：ascat_sample_metrics.tsv

| sample_id | ascat_pair_id              | tumor_purity | ploidy | goodness_of_fit | WGD | GI    | LOH   | n_segs | homdel_segs | homdel_largest | homdel_size | homdel_fraction | mode_minA | mode_majA | tumour_mapd         | normal_mapd | n_het_SNP |
|-----------|----------------------------|--------------|--------|----------------|-----|-------|-------|--------|-------------|----------------|-------------|-----------------|-----------|-----------|---------------------|-------------|-----------|
| HP_tumor  | HP_tumor_vs_HP_normal      | 0.68         | 1.9718 | 99.8811        | 0   | 0.1108 | 0.0858 | 39     |             | 2136712     | 177319         | 0.0001      | 1               | 1         | 5kb=0.0293 | 500bp=0             | 836821      |
| TSDX001   | TSDX001_vs_NSDX001         | 0.8          | 1.9414 | 98.5363        | 0   | 0.2715 | 0.1698 | 41     | 2           | 232677         | 264852      | 0.0001          | 1         | 1          | 100kb=0.1475        | 1kb=0       | 832171    |
| TSDX002   | TSDX002_vs_NSDX002         | 0.58         | 2.0219 | 97.1364        | 0   | 0.0802 | 0.0107 | 68     | 2           | 535503         | 864658      | 0.0003          | 1         | 1          | 5kb=0.4724          | 1kb=0       | 825467    |
| TSDX003   | TSDX003_vs_NSDX003         | 0.95         | 1.9122 | 99.818         | 0   | 0.1142 | 0.1037 | 46     | 0           | 0              | 0.0         | 1               | 1         | 1          | 100kb=0.2484        | 1kb=0       | 833788    |

## 把 ASCAT 结果合并进你的 metadata
### 假设你的 metadata 是：sample_subtype_metadata.tsv
|sample_id|subtype|has_wgs|has_wes|has_rna|tumor_purity|ploidy|age_group|tumorncat|lymph_node_involvement|lymph_node_ratio|status|os|dfs|stage|ihc_syn|ihc_iga|ihc_cd56|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|TSDX001|NP|1|0|1||1|mixed|0|0|0|NA|NA|NA|IB1|negative|negative|positive|
|TSDX002|ES|1|0|1||1|pure|1|0|0|NA|NA|NA|IB2|positive|positive|positive|
|TSDX003|ES|1|0|1||1|pure|0|0|0|NA|1133.7|NA|IIA|positive|positive|negative|
|TSDX004|ES|1|0|1||1|pure|0|0|0|NA|071.33333333|71.33333333|IB1|positive|positive|positive|

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

### 最后得到：sample_subtype_metadata.with_ascat.tsv
| sample_id | subtype | has_wgs | has_wes | has_rna | age_group | tumor | ncat | lymph_node_involvement | lymph_node_ratio | status | os | dfs | stage | ihc_syn | ihc_iga | ihc_cd56 | ascat_pair_id | tumor_purity | ploidy | goodness_of_fit | WGD | GI | LOH | n_segs | homdel_segs | homdel_largest | homdel_size | homdel_fraction | mode_minA | mode_majA | tumour_mapd | normal_mapd | n_het_SNP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TSDX001 | NP | 1 | 0 | 1 | 1 | mixed | 0 | 0.0 | 0.0 | | | IB1 | negative | negative | positive | TSDX001_vs_NSDX001 | 0.8 | 1.9414 | 98.5363 | 0 | 0.2715 | 0.1698 | 41.0 | 2.0 | 232677.0 | 264852.0 | 0.0001 | 1.0 | 1.0 | 100kb=0.1475 | 1kb=0 | 832171.0 |
| TSDX002 | ES | 1 | 0 | 1 | 1 | pure | 1 | 0.0 | 0.0 | | | IB2 | positive | positive | positive | TSDX002_vs_NSDX002 | 0.58 | 2.0219 | 97.1364 | 0 | 0.0802 | 0.0107 | 68.0 | 2.0 | 535503.0 | 864658.0 | 0.0003 | 1.0 | 1.0 | 5kb=0.4724 | 1kb=0 | 825467.0 |



# 病毒片段保留模式(E6/E7保留型、全长型、、)
先准备sample_major_hpv_type.tsv，
从你的 sample_level_integration_summary.tsv 提取：
```bash
awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="sample_id")s=i;if($i=="major_hpv_type")h=i}; print "sample_id","major_hpv_type"; next} {print $s,$h}' \
/data/person/wup/liusy/wgs/scripts/surviurs/analysis_test_out_full/sample_level_integration_summary.tsv \
> sample_major_hpv_type.tsv
```

```bash
cd /data/person/wup/liusy/wgs/scripts/survirus
python3 viral_retention_from_virus_side.py \
  --results_root /data/person/wup/liusy/wgs/results/integration \
  --hpv_bed /data/person/wup/public/liusy_files/reference_genomes/virus/annotation/gff/hpv_annotation_simple.survirus.bed \
  --samples TSDX001 TSDX002 TSDX003 TSDX004 TSDX005 TSDX006 TSDX007 TSDX008 TSDX009 TSDX010 TSDX011 TSDX014 TSDX017 TSDX018 TSDX019 TSDX025 \
  --outdir hpv_retention_virus_side_out \
  --top_only

python3 viral_retention_from_virus_side.py \
  --results_root /data/person/wup/liusy/wgs/results/integration \
  --hpv_bed /data/person/wup/public/liusy_files/reference_genomes/virus/annotation/gff/hpv_annotation_simple.survirus.bed \
  --sample_hpv_type_map sample_major_hpv_type.tsv \
  --samples TSDX001 TSDX002 TSDX003 TSDX004 TSDX005 TSDX006 TSDX007 TSDX008 TSDX009 TSDX010 TSDX011 TSDX014 TSDX017 TSDX018 TSDX019 TSDX025 \
  --outdir hpv_retention_virus_side_out

```
#### 最后得到文件包括

hpv_idxstats.virus_side.tsv ： 每个样本每个 contig 的 mapped reads、用于判断主 HPV 型别。
hpv_retention_pattern.virus_side.tsv：每行是一个样本的一个 HPV 区域
hpv_region_depth.virus_side.tsv：每行是一个样本的一个 HPV 型别

**画figure 3a的图**
你现在至少需要这 5 个输入：

sample_subtype_metadata.with_ascat.tsv
**提供 sample_id, subtype, WGD 等**

sample_level_integration_summary.tsv
**提供 high_conf_event_count, remap_confirmation_rate, major_hpv_type 等**

integration_event_annotation.tsv
**提供 host_region_class, virus_gene_segment 等事件级断点信息**

hpv_region_depth.virus_side.tsv
**viral_retention_from_virus_side.py 的输出，提供每个病毒区域的相对覆盖度**

hpv_retention_pattern.virus_side.tsv
**viral_retention_from_virus_side.py 的输出，提供 E6_retained/E7_retained 等标签**

再运行
```bash
cd /data/person/wup/liusy/wgs/scripts/figures/3a
python3 prepare_figure3A_tables.py \
  --metadata /data/person/wup/liusy/sarek/sample_subtype_metadata.with_ascat.tsv \
  --sample_summary  /data/person/wup/liusy/wgs/scripts/survirus/analysis_test_out_full/sample_level_integration_summary.tsv \
  --event_annotation  /data/person/wup/liusy/wgs/scripts/survirus/analysis_test_out_full/integration_event_annotation.tsv \
  --region_depth  /data/person/wup/liusy/wgs/scripts/survirus/hpv_retention_virus_side_out/hpv_region_depth.virus_side.tsv \
  --retention_pattern  /data/person/wup/liusy/wgs/scripts/survirus/hpv_retention_virus_side_out/hpv_retention_pattern.virus_side.tsv \
  --outdir ./ \
  --drop_na_subtype

```

得到文件
figure3A_plot_table.tsv
figure3A_host_region_prop.tsv
figure3A_virus_breakpoint_prop.tsv
figure3A_viral_retention_depth.tsv
figure3A_viral_retention_flags.tsv

运行画图脚本

```bash
conda activate rplot
export PATH="$CONDA_PREFIX/bin:$(echo "$PATH" | tr ':' '\n' | grep -v '/envs/gatk/bin' | paste -sd: -)"
hash -r
which R
which Rscript
```
```bash
Rscript plot_figure3A_complexheatmap.R
```
输出：Figure3A_HPV_feature_overview.pdf

### 画figure 3B：三分型的 HPV 特征全景图
再运行
```bash
cd /data/person/wup/liusy/wgs/scripts/figures/3b
Rscript plot_Figure3B_driver_oncoprint.R
```

### 画figure 3S2-3个体的突变情况
再运行
```bash
cd /data/person/wup/liusy/wgs/scripts/figures/3b
#画驱动突变基因的oncoplot
Rscript plot_Figure3B_driver_oncoprint.R
#画两个个体病例卡片式表格：左边放突变型样本 TSDX025，右边放 POLE 超突变样本
Rscript plot_FigureS3B_individual_driver_tables.R
#画 POLE 超突变样本的突变summary figure
Rscript plot_FigureS3B_individual_driver_tables.R

```















