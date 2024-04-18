# UniGapClose
基因组Gap Closing工具：使用比对质量最高的区域来消除基因组中的间隙（gaps）。本工具适用于高质量基因组、contigs，以及三代测序数据。

## 使用说明

### 安装依赖

在运行之前，请确保你已经安装了以下依赖项：

- Python 3.7+
- Biopython 1.78+
- minimap2

你可以使用以下命令来安装相关依赖：

```bash
pip install biopython
conda install -c bioconda minimap2
```

### 下载代码
克隆这个仓库到你的本地计算机：
```bash
git clone https://github.com/CongyangY/UniGapClose.git
cd UniGapClose
```

### 软件运行 
#### one-step


#### step by step
1. 基因组比对
```shell
unimap --cs -t 120 -cs asm5 target.chr.fa query.fa > your_alignment.paf
```

2. 生成基因组的gap数据
```bash
python ./GapFinder.py your_genome.fa raw_gaps.bed
```
- your_genome.fa: 基因组DNA序列
- raw_gaps.bed: 统计gaps位置信息的输出文件

3. 找到最佳的跨越gap是的paf records
```bash
python ./FindBestPAFRecord.py your_alignment.paf raw_gaps.bed your_alignment.best.paf
```
- your_alignment.paf: 比对的PAF格式文件
- raw_gaps.bed: 第一步得到的gaps位置信息文件
- your_alignment.best.paf: 算法生成的最佳的PAF records 文件

4. 替换gap区域的序列
```shell
python ./ReplaceGaps.py your_alignment.best.paf target.chr.fa query.fa chr.replaceGaps.fa chr.replaceGaps.gap.bed
```

5. scaffold vs contig (self)