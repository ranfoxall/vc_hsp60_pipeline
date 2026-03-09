# vc_hsp60_pipeline

A Vibrio-centric hsp60 amplicon sequencing pipeline built on DADA2 and QIIME2.

This pipeline processes paired-end Illumina reads generated with the *Vibrio*-centric assay described in King et al. 2019 (*Front. Microbiol.* 10:2907). It produces ASV tables, FASTA files, and a phyloseq object ready for downstream analysis. Taxonomy classification uses the cpn60 QIIME2 classifier and is run as a separate step, allowing DADA2 and QIIME2 to run in their own conda environments.

> **Note:** This pipeline is designed specifically for the *Vibrio*-centric hsp60 assay. It is **not suitable** for universal hsp60 (cpn60) primers.

---

## What you'll need

- A Linux, macOS, or Windows system, or an HPC cluster (tested on UNH Premise with Linux)
- Conda or Miniconda
- Paired-end demultiplexed FASTQ files in standard Illumina format
- The cpn60 QIIME2 classifier (see Step 4)

> **Note:** QIIME2 has limited Windows support. If you are on Windows and want to run the classification step, consider using [WSL (Windows Subsystem for Linux)](https://learn.microsoft.com/en-us/windows/wsl/install) or a Linux HPC cluster.

---

## Step 1: Clone the repository

```bash
git clone https://github.com/yourusername/vc_hsp60_pipeline.git
cd vc_hsp60_pipeline
```

> **UNH Premise:** Run this from your home or project directory, e.g.:
> ```bash
> cd /mnt/home/whistler/yourname
> git clone https://github.com/yourusername/vc_hsp60_pipeline.git
> ```

---

## Step 2: Install Miniconda (if not already installed)

Download and install Miniconda: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

> **UNH Premise:** Conda is already available via the `anaconda/colsa` module — skip this step.

---

## Step 3: Create the conda environment

```bash
conda env create -f environment.yml
conda activate vc_hsp60_pipeline
```

> **UNH Premise:**
> ```bash
> module purge
> module load anaconda/colsa
> conda env create -f environment.yml
> conda activate vc_hsp60_pipeline
> ```

> **Note:** QIIME2 is **not** included in this environment — it runs in its own separate environment. See Step 5.

---

## Step 3b: Install R/Bioconductor packages

Most Bioconductor packages cannot be reliably installed via conda. Run the provided installer **once** after activating the environment:

```bash
conda activate vc_hsp60_pipeline
Rscript install_packages.R
```

> **UNH Premise:**
> ```bash
> module purge
> module load anaconda/colsa
> conda activate vc_hsp60_pipeline
> Rscript install_packages.R
> ```

This installs all packages in the correct order, skips anything already installed, and prints a version summary when done.

<details>
<summary>Troubleshooting package installation errors</summary>

- **`lzma.h: No such file or directory`** (Rhtslib fails):
  ```bash
  conda install -c conda-forge xz
  Rscript install_packages.R
  ```

- **`rhdf5filters` fails with a C compiler conflict** — the `environment.yml` installs this via conda to prevent it, but if it occurs:
  ```bash
  mamba install -c bioconda bioconductor-rhdf5filters bioconductor-rhdf5 bioconductor-biomformat bioconductor-phyloseq
  Rscript install_packages.R
  ```

- **`png` fails on Linux** — install the system library:
  ```bash
  sudo apt install libpng-dev        # Debian/Ubuntu
  sudo dnf install libpng-devel      # RHEL/Rocky/Fedora
  ```
</details>

---

## Step 4: Get QIIME2

QIIME2 is only needed for taxonomy classification (Step 7). It must be installed in its own conda environment — **do not** install it into `vc_hsp60_pipeline`.

> **UNH Premise:** QIIME2 is already available — no installation needed. The SLURM scripts will activate it automatically.

For all other systems, install following the official instructions: [https://docs.qiime2.org](https://docs.qiime2.org). Note the name of your QIIME2 conda environment — you will need it in Step 6.

---

## Step 5: Get the cpn60 classifier

The classifier is needed only for taxonomy classification (Step 7).

> **UNH Premise:** The classifier is already available at:
> ```
> /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
> ```
> The SLURM scripts point to this path automatically — no download needed.

**All other users:** Download from GitHub:
```bash
wget https://github.com/HillLabSask/cpn60-Classifier/releases/download/v11.1/cpn60-q2-feature-classifier-v11.tar.gz
tar -xzf cpn60-q2-feature-classifier-v11.tar.gz
# Classifier will be at: cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
```

> **Important:** The pre-trained classifier was built with QIIME2 v2022.11 (scikit-learn v0.24.1). If your QIIME2 version is different, you may need to retrain it. See the [Troubleshooting](#troubleshooting) section.

---

## Step 6: Prepare your reads

Organize your demultiplexed FASTQ files in a single folder. Files must follow standard Illumina naming:

```
raw_reads/
  Sample1_S1_L001_R1_001.fastq.gz
  Sample1_S1_L001_R2_001.fastq.gz
  Sample2_S2_L001_R1_001.fastq.gz
  Sample2_S2_L001_R2_001.fastq.gz
```

Sample names are extracted by stripping `_S##_L###_R1_001.fastq.gz`, so `1B_1_new_S1_L001_R1_001.fastq.gz` becomes `1B_1_new`.

> Every `_R1_` file must have a matching `_R2_` file.

---

## Step 7: Run the pipeline

There are two ways to run the pipeline — via SLURM (recommended for HPC) or directly on the command line.

### Option A: SLURM (recommended for UNH Premise and other HPC clusters)

Open the SLURM script for your use case and edit the variables at the top:

| Script | Use when |
|--------|----------|
| `run_vchsp60_full.slurm` | Run everything in one job (trimming + DADA2 + classification) |
| `run_vchsp60_dada2.slurm` | Run trimming + DADA2 only, then classify separately |
| `run_vchsp60_classify.slurm` | Run classification only (after DADA2 is done) |

**Variables to set at the top of each script:**

```bash
UNH_PREMISE=1        # set to 1 for UNH Premise (auto-sets environments + classifier path)
                     # set to 0 for other systems and fill in QIIME2_ENV, PIPELINE_ENV, CLASSIFIER

READS_PATH=/path/to/your/raw_reads    # folder containing your .fastq.gz files
OUTPUT_PREFIX=Hsp60_MyRun             # short name for this run — used to name all outputs
WORKDIR=/path/to/your/workdir         # where outputs will be written
SCRIPT_PATH=/path/to/vc_hsp60_pipeline  # where you cloned this repo
```

> **UNH Premise example:**
> ```bash
> UNH_PREMISE=1
> READS_PATH=/mnt/home/whistler/yourname/hsp60_vc/raw_reads
> OUTPUT_PREFIX=Hsp60_MyRun
> WORKDIR=/mnt/home/whistler/yourname/hsp60_vc
> SCRIPT_PATH=/mnt/home/whistler/yourname/vc_hsp60_pipeline
> ```

Then submit:
```bash
sbatch run_vchsp60_full.slurm
# or for separate jobs:
JOB1=$(sbatch --parsable run_vchsp60_dada2.slurm)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 run_vchsp60_classify.slurm)
```

Monitor your jobs:
```bash
squeue -u $USER                          # check job status
tail -f vchsp60_full_<job_id>.log        # stream log output
sacct -u $USER                           # view completed/failed jobs
scancel <job_id>                         # cancel a job
```

---

### Option B: Command line (non-HPC)

**Step 7a — DADA2 pipeline:**
```bash
conda activate vc_hsp60_pipeline

Rscript scripts/vc_hsp60_pipeline.R \
  --reads_path /path/to/raw_reads \
  --output_prefix Hsp60_MyRun
```

**Step 7b — Taxonomy classification:**
```bash
conda activate <your_qiime2_env>

Rscript scripts/vc_hsp60_classify.R \
  --output_prefix Hsp60_MyRun \
  --classifier /path/to/cpn60_classifier_v11.qza
```

> **UNH Premise users** — use this classifier path:
> ```
> --classifier /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
> ```

---

## Outputs

### From `vc_hsp60_pipeline.R`

| File | Description |
|------|-------------|
| `{prefix}_Counts_seqASV_b.tsv` | Counts table with full ASV sequences as row names |
| `{prefix}_Counts_numASV.tsv` | Counts table with numbered ASVs (ASV_1, ASV_2, …) |
| `{prefix}_ASVs.fa` | FASTA of ASVs — input for classification |
| `{prefix}_phyloseq.rds` | Phyloseq object with counts and sample metadata |
| `{prefix}_quality_profiles.pdf` | Quality profiles for representative samples |
| `{prefix}_read_tracking.tsv` | Per-sample read counts at each step: input, filtered, denoised, merged, non-chimeric, % retained |
| `{prefix}_error_model.rds` | Cached error model — reused automatically on reruns |
| `{prefix}_dropped_samples.txt` | List of samples with zero reads after filtering *(if any)* |

### From `vc_hsp60_classify.R`

| File | Description |
|------|-------------|
| `{prefix}_taxonomy.tsv` | Taxonomy assignments from cpn60 classifier |
| `{prefix}_ASVs_counts_taxonomy.tsv` | Counts + taxonomy merged table |
| `{prefix}_phyloseq_taxonomy.rds` | Phyloseq object with taxonomy added |

---

## Restarting a failed job

The SLURM scripts have restart flags so you can skip steps that already completed. Set any flag to `1` to skip that step:

```bash
SKIP_IMPORT=0      # set to 1 to skip QIIME2 import
SKIP_CUTADAPT=0    # set to 1 to skip adapter trimming
SKIP_DADA2=0       # set to 1 to skip DADA2
SKIP_CLASSIFY=0    # set to 1 to skip classification (full script only)
```

When skipping cutadapt, you must also set `TRIMMED_QZA` to your trimmed reads artifact from the previous run:
```bash
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
```

**Common scenario — trimming worked but DADA2 failed:**
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
SKIP_DADA2=0
```

**Common scenario — DADA2 worked but classification failed (`run_vchsp60_full.slurm` only):**
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
SKIP_DADA2=1
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
SKIP_CLASSIFY=0
```

### Reusing a saved error model

Error model learning is the most time-consuming part of DADA2. The pipeline handles this in two ways:

**Automatic caching:** After the first successful run, `{prefix}_error_model.rds` is saved in your working directory. On any rerun with the same `OUTPUT_PREFIX`, it is loaded automatically — no action needed.

**Manual override:** If you want to reuse an error model from a *different* run (e.g. your job failed before the cache was saved, or you want to reuse a model from a previous dataset), set `ERROR_MODEL_RDS` at the top of the SLURM script:
```bash
ERROR_MODEL_RDS=Hsp60_MyRun_error_model.rds
```

Or on the command line:
```bash
Rscript scripts/vc_hsp60_pipeline.R   --reads_path /path/to/reads   --output_prefix Hsp60_MyRun   --error_model_rds Hsp60_MyRun_error_model.rds
```

To force the pipeline to relearn error rates from scratch (e.g. if you change your input reads or truncation lengths), leave `ERROR_MODEL_RDS` blank and delete the auto-cached file:
```bash
rm Hsp60_MyRun_error_model.rds
```

### Checking read loss with the tracking table

After a successful run, `{prefix}_read_tracking.tsv` shows how many reads passed each step per sample:

| sample | reads_in | filtered | denoised | merged | nonchim | pct_retained |
|--------|----------|----------|----------|--------|---------|--------------|
| Sample1 | 150000 | 120000 | 118500 | 112000 | 108000 | 72.0 |

If a sample has very low `pct_retained`, check:
- **Low filtered:** reads may be poor quality or wrong length — check `{prefix}_quality_profiles.pdf`
- **Low merged:** truncation lengths may be too short for forward and reverse reads to overlap
- **Low nonchim:** high chimera rate can indicate primer issues or library prep problems

---

## Merging multiple runs or plates

If you have data from multiple sequencing runs, merge them using the merge script:

```bash
Rscript scripts/merge_vc_hsp60_ASVs.R merged_run \
  results/plate1_Counts_seqASV_b.tsv results/plate1_ASVs.fa \
  results/plate2_Counts_seqASV_b.tsv results/plate2_ASVs.fa
```

> Only sequence-based ASV tables (`_b`) are mergeable — do not merge numbered ASV tables directly.

---

## Using the phyloseq object in R

```R
library(phyloseq)
ps <- readRDS("Hsp60_MyRun_phyloseq.rds")
otu_table(ps)    # ASV counts
sample_data(ps)  # sample metadata
tax_table(ps)    # taxonomy (available after running vc_hsp60_classify.R)
```

---

## Troubleshooting

**Classification fails with a scikit-learn version error:**
The pre-trained classifier was built with scikit-learn 0.24.1 (QIIME2 v2022.11). If your QIIME2 version is different, retrain the classifier:

```bash
conda activate <your_qiime2_env>

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /path/to/cpn60_v11_seqs.fasta \
  --output-path cpn60_v11_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /path/to/cpn60_v11_taxonomy_table.txt \
  --output-path cpn60_v11_taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads cpn60_v11_seqs.qza \
  --i-reference-taxonomy cpn60_v11_taxonomy.qza \
  --o-classifier cpn60_classifier_v11_retrained.qza
```

> **UNH Premise:** The reference files are in the shared directory:
> ```
> /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/
>   cpn60_v11_seqs.fasta
>   cpn60_v11_taxonomy_table.txt
> ```
> This only needs to be done once. Contact [Toni Westbrook](mailto:anthony.westbrook@unh.edu) if you have trouble with module conflicts.

**Samples dropped with zero reads after filtering:**
The pipeline will automatically exclude these samples and write their names to `{prefix}_dropped_samples.txt`. This is usually caused by very low input read counts or poor quality reads. Check the `_quality_profiles.pdf` to assess read quality before running.

**QIIME2 not found in PATH:**
Make sure your QIIME2 conda environment is active before running `vc_hsp60_classify.R`:
```bash
conda activate <your_qiime2_env>
```

---

## Arguments reference

### vc_hsp60_pipeline.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--reads_path` | Yes | — | Path to folder containing raw paired FASTQ files |
| `--output_prefix` | Yes | — | Prefix for all output files |
| `--error_model` | No | `loess` | Error model: `loess` (recommended), `default`, or `compare` (runs both, picks best — slowest) |

### vc_hsp60_classify.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_prefix` | Yes | — | Output prefix used in vc_hsp60_pipeline.R |
| `--classifier` | Yes | — | Path to cpn60 QIIME2 classifier `.qza` |
| `--phyloseq_rds` | No | `{prefix}_phyloseq.rds` | Override phyloseq file path if not in working directory |

---

## Repository structure

```
vc_hsp60_pipeline/
  scripts/
    vc_hsp60_pipeline.R       # DADA2 pipeline
    vc_hsp60_classify.R       # QIIME2 classification
    merge_vc_hsp60_ASVs.R     # Multi-run merge utility
  run_vchsp60_full.slurm      # SLURM: full pipeline
  run_vchsp60_dada2.slurm     # SLURM: trimming + DADA2 only
  run_vchsp60_classify.slurm  # SLURM: classification only
  environment.yml             # conda environment
  install_packages.R          # Bioconductor installer (run once)
  README.md
  CITATION.cff
  LICENSE
```

---

## Citation

If you use this pipeline in your research, please cite it as:

> Foxall, R., Whistler, C. and Jones, S. (2026). *vc_hsp60_pipeline: A Vibrio-centric Hsp60 amplicon sequencing pipeline*. Whistler Lab, University of New Hampshire. https://github.com/yourusername/vc_hsp60_pipeline

A `CITATION.cff` file is included — GitHub will display a **"Cite this repository"** button automatically.

Please also cite the original assay:

> King WL, Siboni N, Kahlke T, Green TJ, Labbate M and Seymour JR (2019). A New High Throughput Sequencing Assay for Characterizing the Diversity of Natural *Vibrio* Communities and Its Application to a Pacific Oyster Mortality Event. *Frontiers in Microbiology* 10:2907. https://doi.org/10.3389/fmicb.2019.02907

### Tool citations

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA and Holmes SP (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581–583. https://doi.org/10.1038/nmeth.3869

> Bolyen E, Rideout JR, Dillon MR et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology* 37:852–857. https://doi.org/10.1038/s41587-019-0209-9

> Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17(1):10–12. https://doi.org/10.14806/ej.17.1.200

> McMurdie PJ and Holmes S (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8(4):e61217. https://doi.org/10.1371/journal.pone.0061217

> R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

---

## Funding and acknowledgements

This pipeline was developed by Randi Foxall in the Whistler Lab at the University of New Hampshire under the supervision of Dr. Cheryl Whistler, with support from Dr. Stephen Jones.

This work was supported by the New Hampshire Agricultural Experiment Station through the NHAES CREATE (Collaborative Research Enhancement and Team Exploration) program, and by NH EPSCoR.

---

## Authors

| Name | Role | Affiliation |
|------|------|-------------|
| Randi Foxall | Developer | Whistler Lab, University of New Hampshire |
| Cheryl Whistler | Principal Investigator & Funder | Department of Molecular, Cellular and Biomedical Sciences, University of New Hampshire |
| Stephen Jones | Co-Investigator & Funder | Department of Natural Resources and the Environment, UNH |

---

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt this material for non-commercial purposes with appropriate credit. Commercial use requires written permission from the authors.
