# SurVirus（HPC/Slurm 批量版）

本仓库在尽量保留原始 SurVirus 核心流程的基础上，新增了适合 HPC/Slurm 的“外层封装”运行方式：
- 使用 `samples.tsv` 管理样本输入
- 支持单样本调试运行
- 支持 Slurm job array 样本级并行
- 保留原始 `surveyor.py` + C++ 可执行程序的主流程与算法

> 说明：原始 SurVirus 为 Python2 程序，请确保运行环境可用 `python2`。

---

## 1. 这个版本改造了什么

### 保留不变
- `surveyor.py` 主流程和所有核心算法逻辑不改。
- 仍然通过 SurVirus 原生流程完成筛选、重比对、断点分析与过滤。

### 新增内容
- `run_survirus_pipeline.py`：统一入口（读取 `samples.tsv`、参数检查、按样本调用 SurVirus）
- `run_survirus_array.slurm`：Slurm array 提交脚本
- `submit_survirus_array.sh`：按实时可用资源自动调优并提交 array（共享集群推荐）
- `run_survirus_single.sh`：单样本快速调试脚本
- `examples/samples.tsv`：示例输入文件

---

## 2. 环境准备

## 2.1 编译 C++ 程序

```bash
./build_libs.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

## 2.2 依赖软件

- Python2（用于 `surveyor.py`）
- Python3（用于新封装脚本 `run_survirus_pipeline.py`）
- NumPy / PyFaidx / PySam
- bwa-mem2（必须；本改造版统一使用 bwa-mem2）
- samtools
- sdust

---

## 3. Reference（默认路径）

本封装默认使用以下参考路径（可通过参数覆盖）：

- host
  - `/data/person/wup/public/liusy_files/reference_genomes/hg38_no_alt/reference/hg38_no_alt.fa`
- virus
  - `/data/person/wup/public/liusy_files/reference_genomes/virus/reference/virus.fa`
- host_virus
  - `/data/person/wup/public/liusy_files/reference_genomes/host_virus/reference/host_virus.fa`

请确保三套 reference 都已建立索引（`bwa-mem2 index` + `samtools faidx`）。

---

## 4. `samples.tsv` 格式说明

`samples.tsv` 为 **TAB 分隔**，必须有表头，至少包含 3 列：

- `sample_id`
- `input_R1`
- `input_R2`

示例：

```tsv
sample_id	input_R1	input_R2
sampleA	/path/to/sampleA.R1.fastp.gz	/path/to/sampleA.R2.fastp.gz
sampleB	/path/to/sampleB.R1.fastp.gz	/path/to/sampleB.R2.fastp.gz
```

### 生成 samples.tsv 的bash脚本
```bash
#!/bin/bash

# 设置数据目录路径
DATA_DIR="/data/person/wup/public/liusy_files/sccc/preprocessed_bam/wgs/fastp"
OUTPUT_FILE="samples.tsv"

# 写入表头
echo -e "sample_id\tinput_R1\tinput_R2" > "$OUTPUT_FILE"

# 提取所有唯一的样本ID
declare -A samples

# 遍历所有 .R1.fastp.gz 文件
for r1_file in "$DATA_DIR"/*.R1.fastp.gz; do
    # 获取文件名（不含路径）
    basename=$(basename "$r1_file")
    # 提取样本ID（去掉 .R1.fastp.gz）
    sample_id=${basename%.R1.fastp.gz}
    
    # 构建对应的R2文件路径
    r2_file="$DATA_DIR/${sample_id}.R2.fastp.gz"
    
    # 检查R2文件是否存在
    if [[ -f "$r2_file" ]]; then
        # 输出到TSV文件
        echo -e "${sample_id}\t${r1_file}\t${r2_file}" >> "$OUTPUT_FILE"
    else
        echo "警告: 找不到 $sample_id 的R2文件" >&2
    fi
done

echo "已生成 $OUTPUT_FILE，包含 $(($(wc -l < "$OUTPUT_FILE") - 1)) 个样本"
```

程序会做以下检查并给出清晰报错：
- `samples.tsv` 是否存在
- 表头是否包含必须列
- 空值
- 重复 `sample_id`
- FASTQ 文件是否存在
- Slurm task id 是否越界
- 输出目录若已存在会自动删除并覆盖（同一样本重跑时直接覆盖旧结果）

---

## 5. 单样本运行（本地调试）

### 方法 A：统一入口直接指定 sample_id

```bash
python3 run_survirus_pipeline.py \
  --samples samples.tsv \
  --outdir results \
  --sample-id sampleA \
  --threads 6 \
  --python2 python2 \
  --surveyor surveyor.py \
  --bwa /path/to/bwa-mem2 \
  --samtools samtools \
  --dust /path/to/sdust
```

### 方法 B：使用便捷脚本

```bash
./run_survirus_single.sh samples.tsv sampleA results 6
```

---

## 6. Slurm array 并行运行

### 6.1 计算样本数量（跳过表头）

```bash
N=$(( $(wc -l < samples.tsv) - 1 ))
```

### 6.2 提交任务

```bash
sbatch --array=1-${N}%8 run_survirus_array.slurm samples.tsv results
# WGS 数据请追加 --wgs
sbatch --array=1-${N}%8 run_survirus_array.slurm samples.tsv results --wgs
```

说明：
- `SLURM_ARRAY_TASK_ID` 为 1-based，对应 `samples.tsv` 数据行（不含表头）
- `%8` 表示最大并发 8 个任务，请根据集群情况调整
- 每个任务只处理 1 个样本
- WES/捕获测序可不加 `--wgs`；WGS 数据建议加 `--wgs` 以匹配 SurVirus README 的推荐参数
- 可通过环境变量指定 bwa-mem2 路径：`BWA_MEM2=/path/to/bwa-mem2 sbatch ...`

### 6.3 推荐：按当前可用资源自动调优提交（共享集群友好）

新增脚本：`submit_survirus_array.sh`。  
它会实时读取 `cpu1/cpu2` 分区当前空闲 CPU，按你设置的“最多使用 1/2 或 1/3”自动计算 `--array` 并发度。

```bash
# 默认使用 1/3 空闲 CPU（适合共享环境）
./submit_survirus_array.sh --samples samples.tsv --outdir results

# 使用 1/2 空闲 CPU（更激进）
./submit_survirus_array.sh --samples samples.tsv --outdir results --fraction 2

# WGS 数据（自动调优脚本也支持 --wgs）
./submit_survirus_array.sh --samples samples.tsv --outdir results --wgs

# 自定义每任务线程和内存
./submit_survirus_array.sh \
  --samples samples.tsv \
  --outdir results \
  --threads-per-task 8 \
  --mem-gb 48 \
  --fraction 3
```

该脚本会输出：样本数、idle CPU、safe CPU、最终 array 并发和完整 sbatch 命令，便于你审计资源占用。

---

## 7. 输出结果与中间文件（按用途整理）

下面把 SurVirus 在每个样本目录（如 `results/sampleA/`）里产生的**最终结果**和**中间文件**集中列出来，并说明它们可用于什么后续分析。

### 7.1 顶层目录

```text
results/
  logs/
    <sample_id>.survirus.log
  <sample_id>/
    ...（该样本全部结果 + 中间文件）
```

另有 Slurm 标准输出（通常在你提交目录下）：

```text
logs/
  slurm-<jobid>_<arrayid>.out
  slurm-<jobid>_<arrayid>.err
```

---

### 7.2 每个样本目录中的关键“最终结果文件”

以下文件位于 `results/<sample_id>/`：

- `results.t1.txt`  
  - **含义**：主结果（默认过滤后保留的整合事件）。  
  - **可用于**：下游统计（每样本事件数）、队列汇总、做 host/virus 断点注释。

- `results.remapped.t1.txt`  
  - **含义**：基于 remapped 断点结果的过滤输出版本。  
  - **可用于**：与 `results.t1.txt` 交叉验证，评估重比对前后事件稳定性。

- `results.discarded.txt`  
  - **含义**：被过滤掉的候选事件。  
  - **可用于**：质控审计（为什么被丢弃）、阈值敏感性分析（调 `minClipSize` / `maxSCDist` 后复查）。

- `results.alternative.txt`  
  - **含义**：候选断点序列在 host 基因组上的替代比对位置（多重匹配线索）。  
  - **可用于**：评估断点唯一性、识别潜在重复区域假阳性。

---

### 7.3 每个样本目录中的中间文件（建议保留）

- `config.txt`  
  - **含义**：本次运行参数快照（threads、clip 阈值、read_len 等）。  
  - **可用于**：复现分析、跨批次参数一致性检查。

- `contig_map`  
  - **含义**：参考序列 contig 名称与内部 ID 映射。  
  - **可用于**：解释二进制/内部索引文件中的 contig 编号。

- `log.txt`  
  - **含义**：`remapper` 阶段标准错误输出。  
  - **可用于**：排查 remap 异常、低支持度候选的诊断。

- `results.txt`、`results.remapped.txt`  
  - **含义**：过滤前的原始候选结果。  
  - **可用于**：自定义过滤策略（自行替换 `filter` 规则）或方法学评估。

- `host_bp_seqs.fa`、`virus_bp_seqs.fa`  
  - **含义**：每个候选断点的宿主/病毒侧共识序列。  
  - **可用于**：做 BLAST、motif、微同源（micro-homology）分析。

- `host_bp_seqs.masked.bed`、`virus_bp_seqs.masked.bed`  
  - **含义**：`sdust` 低复杂度区域掩码结果。  
  - **可用于**：低复杂度过滤、评估序列复杂度对断点置信度的影响。

- `host_bp_seqs.bam`  
  - **含义**：断点序列回贴到 host 的比对结果（用于 alternative 位点提取）。  
  - **可用于**：进一步分析多重定位/重复序列背景。

- `host-side.cs.bam`、`virus-side.cs.bam`  
  - **含义**：跨 workspace 合并后的 host/virus 侧 reads 比对 BAM（坐标排序）。  
  - **可用于**：IGV 可视化、支持 reads 深挖。

- `qnames`、`cigars`、`scores.bin`  
  - **含义**：region-read 关联索引与打分中间数据。  
  - **可用于**：算法调试与性能剖析（一般用户可不直接使用）。

- `readsx/`  
  - **内容**：`<event_id>.bam` + `<event_id>.bam.bai`（每个候选事件一个 BAM）。  
  - **可用于**：按事件做精细人工审阅（最推荐用于 IGV 逐事件验证）。

---

### 7.4 `bam_0/`（FASTQ 模式）中的中间文件清单

当使用 `--fq`（本封装默认）时，样本目录下会有 `bam_0/`。常见文件：

- `stats.txt`：该 workspace 的插入片段上限统计。  
- `retained-pairs_1.fq`、`retained-pairs_2.fq`：保留下来的候选 read pairs。  
- `retained-pairs.remapped.bam`、`retained-pairs.remapped.cs.bam`：候选 pairs 回贴 host+virus 的 BAM。  
- `virus-clips.fa`、`host-clips.fa`：软剪切片段序列（FASTA）。  
- `virus-clips.bam`、`virus-clips.cs.bam`：virus clips 回贴 host 的 BAM。  
- `host-clips.bam`、`host-clips.cs.bam`：host clips 回贴 virus 的 BAM。  
- `virus-anchors.bam`、`virus-anchors.cs.bam`：virus 锚定 reads。  
- `host-anchors.bam`、`host-anchors.cs.bam`：host 锚定 reads。  
- `virus-side.fq`、`host-side.fq`：分类后的 virus/host 侧 reads。  
- `virus-side.bam`、`virus-side.cs.bam`：virus-side 比对结果。  
- `host-side.bam`、`host-side.cs.bam`：host-side 比对结果。  

> 若输入是 BAM/CRAM 而非 FASTQ，会产生 `bam_0/`, `bam_1/` ... 多个 workspace，结构类似。

---

### 7.5 推荐的后续分析路径

1. **先看最终结果**：`results.t1.txt` 与 `results.remapped.t1.txt`。  
2. **做质控与误差排查**：结合 `results.discarded.txt`、`results.alternative.txt`、`*.masked.bed`。  
3. **做可视化验证**：优先使用 `readsx/<event_id>.bam`（按事件最清晰），再结合 `host-side.cs.bam` / `virus-side.cs.bam`。  
4. **做序列层分析**：使用 `host_bp_seqs.fa`、`virus_bp_seqs.fa`（BLAST / motif / micro-homology）。  
5. **做可复现记录**：归档 `config.txt`、`contig_map`、`logs/<sample_id>.survirus.log`、`log.txt`。

---

## 8. 常用参数

`run_survirus_pipeline.py` 关键参数：

- `--samples`：samples.tsv 路径（必填）
- `--outdir`：结果根目录（必填）
- `--sample-id`：指定单样本运行
- `--slurm-array`：使用 `SLURM_ARRAY_TASK_ID` 选择样本
- `--task-id`：手动指定样本行号（本地调试）
- `--threads`：每样本线程数
- `--bwa`：bwa-mem2 路径（默认 `bwa-mem2`）
- `--dry-run`：只打印命令不执行
- `--wgs`：将 `--wgs` 透传给 `surveyor.py`（WGS 数据建议开启；WES/捕获通常不需要）

`submit_survirus_array.sh` 关键参数：
- `--fraction`：最多使用当前空闲 CPU 的 `1/N`（推荐 2 或 3）
- `--threads-per-task`：每个样本任务申请 CPU 核数
- `--mem-gb`：每个样本任务内存
- `--max-parallel`：手动上限（防止并发过高）
- `--wgs`：提交时自动把 `--wgs` 传给 `run_survirus_array.slurm`

---

## 9. 在服务器上下载并解压本 pipeline

如果你是从压缩包分发：

```bash
# 1) 下载（示例）
wget <your_release_url>/survirus2.tar.gz

# 2) 解压
tar -xzf survirus2.tar.gz
cd survirus2

# 3) 编译
./build_libs.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make

# 4) 试运行（dry-run）
python3 run_survirus_pipeline.py --samples examples/samples.tsv --outdir results --sample-id sampleA --dry-run
```

如果你从 Git 仓库拉取：

```bash
git clone <your_repo_url> survirus2
cd survirus2
./build_libs.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

---

## 10. 设计原则

本改造遵循“外层封装 + 最小侵入式修改”：
- 不大规模改动核心算法
- 将批量输入、并行调度、报错友好性放到外层脚本中实现
- 方便初学者快速上手和维护

---

## 11. 常见报错与处理

### 报错：`Classic bwa not found: bwa ...`

这是旧版本脚本的报错（当时还依赖经典 `bwa`）。当前版本已改为直接使用 `bwa-mem2`，并且：
- 若你误传了 `--bwa bwa`，但系统里有 `bwa-mem2`，会自动切换；
- 若系统找不到 `bwa-mem2`，请显式指定路径：

```bash
export BWA_MEM2=/data/person/wup/liusy/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2
sbatch --array=1-${N}%8 run_survirus_array.slurm samples.tsv results
```

### 报错：`extract_clips: No such file or directory` / `reads_categorizer: No such file or directory`

说明 SurVirus C++ 可执行程序未编译完整（或 `surveyor.py` 路径不在编译目录内）。请在仓库根目录重新编译：

```bash
./build_libs.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

### 报错：`ERROR! Unable to open ... .bwt.2bit.64`

说明 reference 没有建立 bwa-mem2 索引。请对 host / virus / host_virus 都执行：

```bash
bwa-mem2 index /path/to/reference.fa
samtools faidx /path/to/reference.fa
```

### 报错：`dust: command not found`

本流程默认使用 `sdust`（不是 `dust`）。请指定 sdust 路径，例如：

```bash
export DUST_EXEC=/data/person/wup/liusy/software/sdust-master/sdust
sbatch --array=1-${N}%8 run_survirus_array.slurm samples.tsv results
```

说明：当前版本 `surveyor.py` 在检测到 `dust` 不存在且 `sdust` 存在时，会自动切换到 `sdust`。
