## 中文版本 (Chinese Version)

**Filename:** `README_CN.md`

# Genome Gap Filler Pipeline

## 目录
1. [简介](#简介)
2. [功能](#功能)
3. [先决条件](#先决条件)
4. [安装](#安装)
    - [使用 Conda 环境](#使用-conda-环境)
5. [使用方法](#使用方法)
    - [运行完整管道](#运行完整管道)
    - [运行单独步骤](#运行单独步骤)
6. [管道步骤概述](#管道步骤概述)
7. [故障排除](#故障排除)
8. [贡献](#贡献)
9. [许可证](#许可证)
10. [联系方式](#联系方式)

---

### 简介

**Genome Gap Filler Pipeline** 是一个全面的工具，旨在识别并填补基因组组装中的缺口。通过利用各种生物信息学工具和高效的数据处理技术，该管道自动化了检测缺口、提取相关读取、生成一致性序列以及将其整合到参考基因组中的过程。

### 功能

- **自动化管道：** 无缝集成从缺口检测到基因组组装的多个步骤。
- **灵活配置：** 允许指定每个缺口的最少读取数，并调整线程数以优化性能。
- **并行处理：** 利用多线程加速计算密集型任务。
- **清晰日志记录：** 提供信息丰富的进度更新和错误消息，以监控管道的执行。
- **简易安装：** 使用 Conda 进行环境管理，确保所有依赖项正确安装。

### 先决条件

在运行管道之前，请确保已安装以下工具和包：

- **生物信息学工具：**
  - [seqkit](https://bioinf.shenwei.me/seqkit/)
  - [hifiasm](https://github.com/chhylp123/hifiasm)
  - [minimap2](https://github.com/lh3/minimap2)

- **Python 包：**
  - `pandas`
  - `tqdm`
  - `pyfaidx`

### 安装

#### 使用 Conda 环境

按照以下步骤设置所需环境：

1. **安装 Conda：**

   如果尚未安装 Conda，请下载并安装 [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 或 [Anaconda](https://www.anaconda.com/products/distribution)。

2. **创建 Conda 环境：**

   您可以使用 YAML 配置文件或直接通过命令行创建环境。

   - **方法 1：使用 YAML 配置文件**

     创建一个名为 `genome_gap_filler_env.yaml` 的文件，内容如下：

     name: genome_gap_filler
     channels:
       - bioconda
       - conda-forge
       - defaults
     dependencies:
       - python=3.9
       - pandas
       - tqdm
       - pyfaidx
       - seqkit
       - hifiasm
       - minimap2

     然后，创建环境：

     conda env create -f genome_gap_filler_env.yaml

   - **方法 2：直接使用命令行**

     conda create -n genome_gap_filler \
       -c bioconda -c conda-forge -c defaults \
       python=3.9 \
       pandas \
       tqdm \
       pyfaidx \
       seqkit \
       hifiasm \
       minimap2

3. **激活环境：**

   conda activate genome_gap_filler

4. **验证安装：**

   - **Python 包：**

     python -c "import pandas, tqdm, pyfaidx; print('Python packages are installed correctly.')"

   - **外部工具：**

     seqkit --version
     hifiasm --version
     minimap2 --version

### 使用方法

#### 运行完整管道

要执行从缺口过滤到基因组更新的整个管道，请使用 `full_pipeline` 子命令：

./genome_gap_filler.py full_pipeline \
  --input_paf input.paf \
  --gaps_bed gaps.bed \
  --fastq_files read1.fastq read2.fastq \
  --consensus_output_dir consensus_seqs \
  --final_output updated_reference.fasta \
  --threads 8 \
  --min_reads 2

**参数说明：**

- `--input_paf`: 初始 PAF 文件的路径。
- `--gaps_bed`: 包含缺口信息的 BED 文件路径。
- `--fastq_files`: 输入 FASTQ 文件列表（以空格分隔）。
- `--consensus_output_dir`: 用于存储一致性 FASTA 文件的目录。
- `--final_output`: 最终包含填补缺口的新参考基因组 FASTA 文件路径。
- `--threads`: 并行处理的线程数。
- `--min_reads`: 每个缺口所需的最少读取数（默认值为 1）。

#### 运行单独步骤

您也可以使用相应的子命令单独运行每个步骤。

- **步骤 0：过滤跨越缺口的 PAF 读取**

  ./genome_gap_filler.py step0_filter_paf \
    --input_paf input.paf \
    --gaps_bed gaps.bed \
    --output_paf filtered.paf

- **步骤 1：按每个缺口拆分 PAF 文件**

  ./genome_gap_filler.py step1_split_paf \
    --filtered_paf filtered.paf \
    --gaps_bed gaps.bed \
    --output_dir split_paf \
    --min_reads 2

- **步骤 2：从 FASTQ 文件中提取并拆分读取**

  ./genome_gap_filler.py step2_extract_split_fastq \
    --paf_dir split_paf \
    --fastq_files read1.fastq read2.fastq \
    --output_dir gap_fastq \
    --temp_reads_list temp_reads_list.txt \
    --extracted_fastq extracted_reads.fastq

- **步骤 3：为每个缺口生成一致性序列**

  ./genome_gap_filler.py step3_generate_consensus \
    --input_dir gap_fastq \
    --output_dir consensus_seqs \
    --threads 8

- **步骤 4：替换参考基因组中的缺口**

  ./genome_gap_filler.py step4_replace_gaps \
    --reference reference.fasta \
    --gaps gaps.bed \
    --consensus_dir consensus_seqs \
    --output updated_reference.fasta \
    --threads 8

### 管道步骤概述

1. **步骤 0：过滤跨越缺口的 PAF 读取**
   - 从 PAF 文件中过滤跨越指定缺口的读取。

2. **步骤 1：按每个缺口拆分 PAF 文件**
   - 根据 BED 文件将过滤后的 PAF 文件拆分为缺口特定的 PAF 文件。
   - 仅为具有最少读取数的缺口创建 PAF 文件。

3. **步骤 2：从 FASTQ 文件中提取并拆分读取**
   - 从提供的 FASTQ 文件中提取在拆分后的 PAF 文件中列出的读取。
   - 将提取的读取拆分为缺口特定的 FASTQ 文件。

4. **步骤 3：为每个缺口生成一致性序列**
   - 使用 `hifiasm` 组装每个缺口的提取读取。
   - 从组装中生成一致性序列。

5. **步骤 4：替换参考基因组中的缺口**
   - 使用 `minimap2` 将一致性序列比对到参考基因组。
   - 将参考基因组中的缺口区域替换为一致性序列。

### 故障排除

- **`gap_fastq` 目录为空：**
  - 确保步骤 1 成功生成了符合最少读取数要求的 PAF 文件。
  - 验证这些 PAF 文件中是否包含有效的读取条目。

- **NameError 问题：**
  - 确保您使用的是最新版本的 `genome_gap_filler.py` 脚本。
  - 在运行子命令时提供所有必需的参数。

- **工具执行错误：**
  - 验证所有外部工具（`seqkit`、`hifiasm`、`minimap2`）是否正确安装并在 Conda 环境中可访问。
  - 检查工具版本以确保兼容性。

- **权限问题：**
  - 确保您具有对所有输入和输出目录的读写权限。

### 贡献

欢迎贡献！请按照以下步骤进行：

1. Fork 本仓库。
2. 为您的功能或修复创建一个新分支。
3. 提交更改，附上清晰的提交信息。
4. 打开一个 pull request，详细说明您的更改。

### 许可证

本项目采用 [MIT 许可证](LICENSE)。

### 联系方式

如有任何问题或需要支持，请联系 [yicongyang@genetics.ac.cn](mailto:yicongyang@genetics.ac.cn)。

---

