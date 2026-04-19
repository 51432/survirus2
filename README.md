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

程序会做以下检查并给出清晰报错：
- `samples.tsv` 是否存在
- 表头是否包含必须列
- 空值
- 重复 `sample_id`
- FASTQ 文件是否存在
- Slurm task id 是否越界
- 输出目录是否已存在（默认不覆盖，需 `--force`）

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
  --dust dust
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
```

说明：
- `SLURM_ARRAY_TASK_ID` 为 1-based，对应 `samples.tsv` 数据行（不含表头）
- `%8` 表示最大并发 8 个任务，请根据集群情况调整
- 每个任务只处理 1 个样本
- 可通过环境变量指定 bwa-mem2 路径：`BWA_MEM2=/path/to/bwa-mem2 sbatch ...`

### 6.3 推荐：按当前可用资源自动调优提交（共享集群友好）

新增脚本：`submit_survirus_array.sh`。  
它会实时读取 `cpu1/cpu2` 分区当前空闲 CPU，按你设置的“最多使用 1/2 或 1/3”自动计算 `--array` 并发度。

```bash
# 默认使用 1/3 空闲 CPU（适合共享环境）
./submit_survirus_array.sh --samples samples.tsv --outdir results

# 使用 1/2 空闲 CPU（更激进）
./submit_survirus_array.sh --samples samples.tsv --outdir results --fraction 2

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

## 7. 输出目录结构

默认输出目录结构如下：

```text
results/
  logs/
    sampleA.survirus.log
    sampleB.survirus.log
  sampleA/
    ...（SurVirus 原始输出）
  sampleB/
    ...（SurVirus 原始输出）
```

同时 Slurm 也会生成：

```text
logs/
  slurm-<jobid>_<arrayid>.out
  slurm-<jobid>_<arrayid>.err
```

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
- `--force`：允许复用已存在样本输出目录
- `--dry-run`：只打印命令不执行

`submit_survirus_array.sh` 关键参数：
- `--fraction`：最多使用当前空闲 CPU 的 `1/N`（推荐 2 或 3）
- `--threads-per-task`：每个样本任务申请 CPU 核数
- `--mem-gb`：每个样本任务内存
- `--max-parallel`：手动上限（防止并发过高）

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
