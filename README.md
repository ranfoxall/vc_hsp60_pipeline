# vc_hsp60_pipeline README

## Overview

This pipeline processes Hsp60 amplicon sequencing data using DADA2, producing ASVs, counts tables, FASTA files, and a phyloseq object ready for downstream analysis. It can also merge multiple runs or plates into a single dataset.

---

## Step 1: Clone the repository

```bash
# Clone the pipeline repository
git clone https://github.com/yourusername/vc_hsp60_pipeline.git
cd vc_hsp60_pipeline
```

This will create a directory called `vc_hsp60_pipeline` containing the `scripts/` folder and any test data.

---

## Step 2: Install Miniconda (if not already installed)

Download Miniconda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html) and follow instructions for your OS.  

---

## Step 3: Create the conda environment

```bash
# From within the repository directory
conda env create -f environment.yml
conda activate vc_hsp60_pipeline
```

This installs R, required R packages, and system dependencies.  

> **Note:** QIIME2 is **not included** in this environment. See Step 4.

---

## Step 4: Install QIIME2 separately

QIIME2 requires its own conda environment. Follow official instructions:  
[https://docs.qiime2.org/2023.11/install/](https://docs.qiime2.org/2023.11/install/)

After installation, make sure the `qiime` command is available in your PATH when you run the pipeline.

---

## Step 5: Prepare your reads folder

Organize your raw FASTQ files in a folder. Example structure:

```
raw_reads/
  Sample1_S1_L001_R1_001.fastq.gz
  Sample1_S1_L001_R2_001.fastq.gz
  Sample2_S1_L001_R1_001.fastq.gz
  Sample2_S1_L001_R2_001.fastq.gz
  ...
```

> The filenames must follow this format for the script to correctly pair forward/reverse reads.

---

## Step 6: Run the main Hsp60 pipeline

```bash
Rscript scripts/run_vc_hsp60_pipeline.R --reads_path /full/path/to/raw_reads --output_prefix Hsp60_OvN
```

**Outputs:**

File | Description
---- | -----------
Hsp60_OvN_SeqTab.csv | Sequence table filtered by target length
Hsp60_OvN_Counts_ASV_num.tsv | Counts table with numbered ASVs (ASV_1, ASV_2, …)
Hsp60_OvN_ASVs.fa | FASTA of ASVs (for taxonomy assignment)
Hsp60_OvN_phyloseq.rds | Phyloseq object with numbered ASVs and sample metadata

---

## Step 7: (Optional) Merge Multiple Plates / Runs

If you have multiple sequencing runs or plates:

```bash
Rscript scripts/merge_vc_hsp60_ASVs.R merged_run \
    tests/test_plate1_Counts_seqASV_b.tsv tests/test_plate1_ASVs.fa \
    tests/test_plate2_Counts_seqASV_b.tsv tests/test_plate2_ASVs.fa
```

**Outputs:**  

File | Description
---- | -----------
merged_run_merged_ASVs_b.tsv | Counts table with ASV sequences as row names
merged_run_merged_ASVs_b.fa  | Merged ASV FASTA
merged_run_merged_ASVs_num.tsv | Counts table with numbered ASVs
merged_run_merged_phyloseq.rds | Phyloseq object with sample metadata

> Notes: Only sequence-based ASVs (_b) are mergeable; do not merge numbered ASV tables.

---

## Step 8: Use phyloseq object in R

```R
library(phyloseq)
ps <- readRDS("Hsp60_OvN_phyloseq.rds")
otu_table(ps)
sample_data(ps)
```

This object contains:
- `otu_table()` : ASV counts (numbered ASVs)
- `sample_data()` : SampleID and TotalReads
- `tax_table()` : can be filled after taxonomy assignment

---

## Directory structure summary

```
vc_hsp60_pipeline/
  scripts/
    run_vc_hsp60_pipeline.R
    merge_vc_hsp60_ASVs.R
  tests/
    example FASTQ files (optional)
  environment.yml
  README.md
```
