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
python run_survirus_annotation.py \
  --results_root /path/to/results \
  --gene_bed genes.bed \
  --tss_bed tss.bed \
  --promoter_bed promoter_2kb.bed \
  --enhancer_bed enhancers.bed \
  --hpv_bed hpv_annot.bed \
  --sample_manifest sample_manifest.tsv \
  --outdir downstream_annotation
```

`sample_manifest` 可选（仅需 `sample_id` 列）；未提供时自动扫描 `results_root` 下所有子目录。

## 输出

- `integration_event_annotation.tsv`（事件级整合注释主表，含 `nearest_transcript` / `nearest_gene` / `distance_to_tss`）
- `survirus_host_breakpoints.bed`
- `survirus_virus_breakpoints.bed`
- `survirus_annotation.log`

## 依赖

- `python >= 3.8`
- `pandas`
- `bedtools`（用于 intersect / closest 注释）
