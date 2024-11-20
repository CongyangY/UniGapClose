**Filename:** `README_EN.md`

# Genome Gap Filler Pipeline

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Prerequisites](#prerequisites)
4. [Installation](#installation)
    - [Using Conda Environment](#using-conda-environment)
5. [Usage](#usage)
    - [Running the Full Pipeline](#running-the-full-pipeline)
    - [Running Individual Steps](#running-individual-steps)
6. [Pipeline Steps Overview](#pipeline-steps-overview)
7. [Troubleshooting](#troubleshooting)
8. [Contributing](#contributing)
9. [License](#license)
10. [Contact](#contact)

---

### Introduction

The **Genome Gap Filler Pipeline** is a comprehensive tool designed to identify and fill gaps in genome assemblies. By leveraging various bioinformatics tools and efficient data processing techniques, this pipeline automates the process of detecting gaps, extracting relevant reads, generating consensus sequences, and integrating them into the reference genome.

### Features

- **Automated Pipeline:** Seamlessly integrates multiple steps from gap detection to genome assembly.
- **Flexible Configuration:** Allows specifying minimum reads per gap and adjusting thread counts for performance optimization.
- **Parallel Processing:** Utilizes multi-threading to speed up computationally intensive tasks.
- **Clear Logging:** Provides informative progress updates and error messages to monitor the pipeline's execution.
- **Easy Installation:** Utilizes Conda for environment management, ensuring all dependencies are correctly installed.

### Prerequisites

Before running the pipeline, ensure that the following tools and packages are installed:

- **Bioinformatics Tools:**
  - [seqkit](https://bioinf.shenwei.me/seqkit/)
  - [hifiasm](https://github.com/chhylp123/hifiasm)
  - [minimap2](https://github.com/lh3/minimap2)

- **Python Packages:**
  - `pandas`
  - `tqdm`
  - `pyfaidx`

### Installation

#### Using Conda Environment

To set up the necessary environment, follow these steps:

1. **Install Conda:**

   If you haven't installed Conda yet, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

2. **Create the Conda Environment:**

   You can create the environment using a YAML configuration file or directly via the command line.

   - **Method 1: Using a YAML Configuration File**

     Create a file named `genome_gap_filler_env.yaml` with the following content:

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

     Then, create the environment:

     conda env create -f genome_gap_filler_env.yaml

   - **Method 2: Using the Command Line Directly**

     conda create -n genome_gap_filler \
       -c bioconda -c conda-forge -c defaults \
       python=3.9 \
       pandas \
       tqdm \
       pyfaidx \
       seqkit \
       hifiasm \
       minimap2

3. **Activate the Environment:**

   conda activate genome_gap_filler

4. **Verify Installations:**

   - **Python Packages:**

     python -c "import pandas, tqdm, pyfaidx; print('Python packages are installed correctly.')"

   - **External Tools:**

     seqkit --version
     hifiasm --version
     minimap2 --version

### Usage

#### Running the Full Pipeline

To execute the entire pipeline from gap filtering to genome updating, use the `full_pipeline` subcommand:

./genome_gap_filler.py full_pipeline \
  --input_paf input.paf \
  --gaps_bed gaps.bed \
  --fastq_files read1.fastq read2.fastq \
  --consensus_output_dir consensus_seqs \
  --final_output updated_reference.fasta \
  --threads 8 \
  --min_reads 2

**Parameter Descriptions:**

- `--input_paf`: Path to the initial PAF file.
- `--gaps_bed`: Path to the BED file containing gap information.
- `--fastq_files`: List of input FASTQ files (space-separated).
- `--consensus_output_dir`: Directory to store consensus FASTA files.
- `--final_output`: Path for the final reference genome with gaps replaced.
- `--threads`: Number of threads for parallel processing.
- `--min_reads`: Minimum number of reads required per gap (default is 1).

#### Running Individual Steps

You can also run each step separately using the corresponding subcommands.

- **Step 0: Filter PAF Reads that Span Gaps**

  ./genome_gap_filler.py step0_filter_paf \
    --input_paf input.paf \
    --gaps_bed gaps.bed \
    --output_paf filtered.paf

- **Step 1: Split PAF by Each Gap**

  ./genome_gap_filler.py step1_split_paf \
    --filtered_paf filtered.paf \
    --gaps_bed gaps.bed \
    --output_dir split_paf \
    --min_reads 2

- **Step 2: Extract and Split Reads from FASTQ Files**

  ./genome_gap_filler.py step2_extract_split_fastq \
    --paf_dir split_paf \
    --fastq_files read1.fastq read2.fastq \
    --output_dir gap_fastq \
    --temp_reads_list temp_reads_list.txt \
    --extracted_fastq extracted_reads.fastq

- **Step 3: Generate Consensus Sequences for Each Gap**

  ./genome_gap_filler.py step3_generate_consensus \
    --input_dir gap_fastq \
    --output_dir consensus_seqs \
    --threads 8

- **Step 4: Replace Gaps in Reference Genome**

  ./genome_gap_filler.py step4_replace_gaps \
    --reference reference.fasta \
    --gaps gaps.bed \
    --consensus_dir consensus_seqs \
    --output updated_reference.fasta \
    --threads 8

### Pipeline Steps Overview

1. **Step 0: Filter PAF Reads that Span Gaps**
   - Filters reads from the PAF file that span the specified gaps.

2. **Step 1: Split PAF by Each Gap**
   - Splits the filtered PAF file into gap-specific PAF files based on the BED file.
   - Only creates PAF files for gaps with a minimum number of supporting reads.

3. **Step 2: Extract and Split Reads from FASTQ Files**
   - Extracts reads listed in the split PAF files from the provided FASTQ files.
   - Splits the extracted reads into gap-specific FASTQ files.

4. **Step 3: Generate Consensus Sequences for Each Gap**
   - Assembles the extracted reads for each gap using `hifiasm`.
   - Generates consensus sequences from the assemblies.

5. **Step 4: Replace Gaps in Reference Genome**
   - Aligns consensus sequences to the reference genome using `minimap2`.
   - Replaces the gap regions in the reference genome with the consensus sequences.

### Troubleshooting

- **Empty `gap_fastq` Directory:**
  - Ensure that Step 1 successfully generated PAF files with the required minimum number of reads.
  - Verify that the PAF files contain valid read entries.

- **NameError Issues:**
  - Make sure you are using the latest version of the `genome_gap_filler.py` script.
  - Ensure all required parameters are provided when running subcommands.

- **Tool Execution Errors:**
  - Verify that all external tools (`seqkit`, `hifiasm`, `minimap2`) are correctly installed and accessible in the Conda environment.
  - Check the versions of the tools to ensure compatibility.

- **Permission Issues:**
  - Ensure you have the necessary read/write permissions for all input and output directories.

### Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Commit your changes with clear commit messages.
4. Open a pull request detailing your changes.

### License

This project is licensed under the [MIT License](LICENSE).



### Contact

For any questions or support, please contact [yicongyang@genetics.ac.cn](mailto:yicongyang@genetics.ac.cn).

---
