# step: analysis

这是一个独立于 SurVirus caller 的下游分析步骤（固定步骤名：`analysis`）。

## 功能

- 从 `manifest TSV` 自动读取每个样本的 SurVirus 结果及可选多组学输入。
- 自动输出：
  - `sample_level_integration_summary.tsv`
  - `integration_event_annotation.tsv`
  - `integration_local_multiomics.tsv`
  - BED 中间文件（host/virus 断点与 50/100/500kb 窗口）
  - QC 表与 `summary_report.md`
- 支持缺失列/缺失文件自动跳过，不因单样本缺失而失败。

## 输入 manifest

必填列：
- `sample_id`
- `survirus_results_t1`

选填列：
- `subtype`
- `survirus_results_remapped_t1`
- `survirus_results_alternative`
- `ascat_segments`
- `manta_sv`
- `tiddit_sv`
- `rna_expr`
- `sample_ploidy`
- `sample_wgd_status`

## 运行方式

```bash
python -m analysis.run_analysis \
  --manifest sample_manifest.tsv \
  --gene_gtf genes.gtf \
  --tss_bed genes_tss.bed \
  --promoter_bed promoters_2kb.bed \
  --enhancer_bed enhancers.bed \
  --hpv_bed hpv_annotation.bed \
  --outdir analysis_out
```

## 可调参数（示例）

- high-confidence 规则参数：
  - `--hc_min_split_reads`（默认 1）
  - `--hc_min_supporting_pairs`（默认 5）
  - `--hc_require_remap`（默认 1）
  - `--hc_require_no_alt_extra`（默认 1）
- cis effect 规则参数：
  - `--cis_distance_bp`（默认 100000）
  - `--cis_z_cutoff`（默认 2.0）

## 最小测试样例

见 `analysis/testdata/`，包含：
- 最小 manifest
- SurVirus t1/remapped/alternative
- 简化基因注释与 HPV 注释
- RNA 表达矩阵

示例：

```bash
python -m analysis.run_analysis \
  --manifest analysis/testdata/manifest.tsv \
  --gene_gtf analysis/testdata/genes.gtf \
  --tss_bed analysis/testdata/tss.bed \
  --promoter_bed analysis/testdata/promoter.bed \
  --hpv_bed analysis/testdata/hpv.bed \
  --outdir analysis/testdata/out
```
