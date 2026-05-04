# vc_hsp60_pipeline

A Vibrio-centric hsp60 amplicon sequencing pipeline built on DADA2 and QIIME2.

This pipeline processes paired-end Illumina reads generated with the *Vibrio*-centric assay described in King et al. 2019 (*Front. Microbiol.* 10:2907). It produces ASV tables, FASTA files, and a phyloseq object ready for downstream analysis. Taxonomy classification uses a two-step approach matching King et al. 2019: ASVs are first screened against a curated Vibrio hsp60 reference set using BLAST (≥90% identity, ≥90% coverage), then Vibrio ASVs are classified with a Vibrio-specific sklearn classifier. Non-Vibrio ASVs are optionally classified with the universal cpn60 classifier. Two phyloseq objects are produced: one containing only Vibrio ASVs for community analysis, and one containing all ASVs for estimating Vibrio relative abundance within the broader community.

> **Note:** This pipeline is designed specifically for the *Vibrio*-centric hsp60 assay. It is **not suitable** for universal hsp60 (cpn60) primers.

---

## Quick start

| I am a... | Start here |
|-----------|-----------|
| **UNH Premise user** | Skip to [UNH Premise quick start](#unh-premise-quick-start) |
| **New user on another system** | Follow the full setup steps below |
| **Returning user re-running an analysis** | See [Restarting a failed job](#restarting-a-failed-job) |

---

## What you'll need

- A Linux, macOS, or Windows system, or an HPC cluster (tested on UNH Premise)
- Conda or Miniconda
- Paired-end demultiplexed FASTQ files in standard Illumina naming format
- The cpn60 QIIME2 classifier (see [Step 5](#step-5-get-the-cpn60-classifier))

> **Windows users:** QIIME2 has limited Windows support. For the classification step, use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) or a Linux HPC cluster.

---

## UNH Premise quick start

If you are working on UNH Premise, environments and the classifier are already set up. You just need to:

**1. Clone the repo:**
```bash
cd /mnt/home/whistler/yourname
git clone https://github.com/ranfoxall/vc_hsp60_pipeline.git
```

**2. Create the pipeline environment (once only):**
```bash
module purge
module load anaconda/colsa
conda env create -f vc_hsp60_pipeline/environment.yml
conda activate vc_hsp60_pipeline
Rscript vc_hsp60_pipeline/install_packages.R
```

**3. Edit and submit the SLURM script:**
```bash
cd /mnt/home/whistler/yourname/vc_hsp60_pipeline
cp run_vchsp60_full.slurm my_run.slurm
# Edit my_run.slurm — set READS_PATH, OUTPUT_PREFIX, WORKDIR, SCRIPT_PATH
sbatch my_run.slurm
```

The SLURM scripts set `UNH_PREMISE=1` by default, which automatically configures the QIIME2 environment and classifier path. See [Step 7](#step-7-run-the-pipeline) for details on what variables to set.

---

## Full setup (all systems)

### Step 1: Clone the repository

```bash
git clone https://github.com/ranfoxall/vc_hsp60_pipeline.git
cd vc_hsp60_pipeline
```

---

### Step 2: Install Miniconda (if not already installed)

Download and install Miniconda: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

> **UNH Premise:** Conda is already available via the `anaconda/colsa` module — skip this step.

---

### Step 3: Create the pipeline conda environment

This environment contains R, DADA2, and all required packages for the core pipeline.

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

> **Note:** QIIME2 is **not** included in this environment — it runs in its own separate environment (see Step 4).

---

### Step 3b: Install R/Bioconductor packages

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

### Step 4: Get QIIME2

QIIME2 is only needed for taxonomy classification. It must be installed in its own conda environment — **do not** install it into `vc_hsp60_pipeline`.

> **UNH Premise:** QIIME2 (`qiime2-2020.2`) is already available — no installation needed. The SLURM scripts activate it automatically.

**All other systems:** Install following the official instructions at [https://docs.qiime2.org](https://docs.qiime2.org). Note the name of your QIIME2 conda environment — you will need it in Step 7.

---

### Step 5: Get the classification reference files

Classification requires three files: a BLAST reference FASTA, a Vibrio-only sklearn classifier, and the universal cpn60 classifier for non-Vibrio ASVs.

> **UNH Premise:** All files are already available — no download needed. The SLURM scripts point to them automatically:
> ```
> # BLAST reference (King et al. 2019 Vibrio hsp60 repset)
> /mnt/home/whistler/foxall/hsp60_ref/vc_hsp60_ref/repset_final_130219.fas
>
> # Vibrio-only classifier (trained from repset)
> /mnt/home/whistler/foxall/hsp60_ref/vc_hsp60_ref/vc_hsp60_classifier_sklearn022.qza
>
> # Universal cpn60 classifier (for non-Vibrio ASVs)
> /mnt/home/whistler/shared/cpn60-Classifier/cpn60_classifier_v11_sklearn142.qza
> ```

**All other users:**

Download the Vibrio hsp60 reference set from the supplemental data of King et al. 2019 — you need `repset_final_130219.fas` and `Vibrio_taxonomy.txt`. Then train the Vibrio-only classifier:

```bash
conda activate <your_qiime2_env>

# deduplicate taxonomy file (required)
python3 -c "
import sys
seen = set()
for line in open('Vibrio_taxonomy.txt'):
    acc = line.split('\t')[0]
    if acc not in seen:
        seen.add(acc)
        sys.stdout.write(line)
" > Vibrio_taxonomy_dedup.txt

# import reference files
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path repset_final_130219.fas \
  --output-path repset_final_130219.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path Vibrio_taxonomy_dedup.txt \
  --output-path Vibrio_taxonomy.qza

# train classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads repset_final_130219.qza \
  --i-reference-taxonomy Vibrio_taxonomy.qza \
  --o-classifier vc_hsp60_classifier.qza
```

Download the universal cpn60 classifier from GitHub:
```bash
wget https://github.com/HillLabSask/cpn60-Classifier/releases/download/v11.1/cpn60-q2-feature-classifier-v11.tar.gz
tar -xzf cpn60-q2-feature-classifier-v11.tar.gz
```

> **Important:** Both classifiers must match the scikit-learn version in your QIIME2 environment. See [Troubleshooting](#troubleshooting) if you get a version mismatch error.

---

### Step 6: Prepare your reads

Organize your demultiplexed FASTQ files in a single folder. Files must follow standard Illumina naming:

```
raw_reads/
  Sample1_S1_L001_R1_001.fastq.gz
  Sample1_S1_L001_R2_001.fastq.gz
  Sample2_S2_L001_R1_001.fastq.gz
  Sample2_S2_L001_R2_001.fastq.gz
```

Sample names are extracted automatically by stripping `_S##_L###_R1_001.fastq.gz`. For example, `1B_1_new_S1_L001_R1_001.fastq.gz` becomes sample `1B_1_new`.

> Every `_R1_` file must have a matching `_R2_` file or it will be skipped.

---

### Step 7: Run the pipeline

There are two ways to run the pipeline — via SLURM (recommended for HPC clusters) or directly on the command line.

#### Option A: SLURM (recommended for UNH Premise and other HPC clusters)

Choose the SLURM script that fits your situation:

| Script | Use when |
|--------|----------|
| `run_vchsp60_full.slurm` | Run everything in one job: trimming + DADA2 + classification |
| `run_vchsp60_dada2.slurm` | Run trimming + DADA2 only (classify separately when done) |
| `run_vchsp60_classify.slurm` | Run classification only (after DADA2 is already done) |

**Edit the variables at the top of the script before submitting.** The key variables are:

```bash
UNH_PREMISE=1     # set to 1 for UNH Premise (auto-configures environments + classifier)
                  # set to 0 for other systems, then fill in QIIME2_ENV, PIPELINE_ENV, CLASSIFIER below

READS_PATH=/path/to/your/raw_reads   # folder containing your .fastq.gz files
OUTPUT_PREFIX=Hsp60_MyRun            # short name for this run — all output files start with this
WORKDIR=/path/to/your/workdir        # where outputs will be written
SCRIPT_PATH=/path/to/vc_hsp60_pipeline  # where you cloned this repo

POOL=pseudo       # DADA2 pooling: pseudo (recommended), FALSE (fastest), TRUE (most sensitive)
METHOD=vibrio     # classification: vibrio (recommended) or sklearn
CONFIDENCE=0.7    # taxonomy confidence threshold (0–1); lower = more assignments, less certainty
```

> **UNH Premise example:**
> ```bash
> UNH_PREMISE=1
> READS_PATH=/mnt/home/whistler/yourname/hsp60_vc/raw_reads
> OUTPUT_PREFIX=Hsp60_MyRun
> WORKDIR=/mnt/home/whistler/yourname/hsp60_vc
> SCRIPT_PATH=/mnt/home/whistler/yourname/vc_hsp60_pipeline
> ```

**Submit the job:**
```bash
sbatch run_vchsp60_full.slurm
```

Or run DADA2 and classification as separate dependent jobs:
```bash
JOB1=$(sbatch --parsable run_vchsp60_dada2.slurm)
sbatch --dependency=afterok:$JOB1 run_vchsp60_classify.slurm
```

**Monitor your jobs:**
```bash
squeue -u $USER                           # see running jobs
tail -f vchsp60_full_<job_id>.log         # stream log output live
sacct -j <job_id> --format=JobID,State,Elapsed,MaxRSS   # check memory and status
scancel <job_id>                          # cancel a job
```

---

#### Option B: Command line (non-HPC)

**Step 7a — DADA2 pipeline:**
```bash
conda activate vc_hsp60_pipeline

Rscript scripts/vc_hsp60_pipeline.R \
  --reads_path /path/to/raw_reads \
  --output_prefix Hsp60_MyRun \
  --pool pseudo
```

**Step 7b — Taxonomy classification:**
```bash
conda activate vc_hsp60_pipeline

# vibrio method (default — recommended)
Rscript scripts/vc_hsp60_classify.R \
  --output_prefix Hsp60_MyRun \
  --method vibrio \
  --blast_db /path/to/repset_final_130219.fas \
  --vibrio_classifier /path/to/vc_hsp60_classifier.qza \
  --classifier /path/to/cpn60_classifier_v11.qza \
  --confidence 0.7

# sklearn method (universal classifier only)
Rscript scripts/vc_hsp60_classify.R \
  --output_prefix Hsp60_MyRun \
  --method sklearn \
  --classifier /path/to/cpn60_classifier_v11.qza \
  --confidence 0.7
```

**Step 7c — Add taxonomy to phyloseq object:**
```bash
conda activate vc_hsp60_pipeline

Rscript scripts/vc_hsp60_add_taxonomy.R \
  --output_prefix Hsp60_MyRun \
  --method vibrio
```

For the vibrio method this produces two phyloseq objects — `{prefix}_phyloseq_vibrio.rds` and `{prefix}_phyloseq_full.rds`. For the sklearn method it produces `{prefix}_phyloseq_taxonomy.rds`.

---

## Outputs

### From `vc_hsp60_pipeline.R`

| File | Description |
|------|-------------|
| `{prefix}_Counts_seqASV_b.tsv` | ASV counts table with full sequences as row names |
| `{prefix}_Counts_numASV.tsv` | ASV counts table with numbered ASVs (ASV_1, ASV_2, …) |
| `{prefix}_ASVs.fa` | FASTA of ASV sequences — input for classification |
| `{prefix}_phyloseq.rds` | Phyloseq object with counts and sample metadata |
| `{prefix}_quality_profiles.pdf` | Read quality profiles for representative samples |
| `{prefix}_read_tracking.tsv` | Per-sample read counts at each step (see below) |
| `{prefix}_error_model.rds` | Cached error model — reloaded automatically on reruns |
| `{prefix}_dada_objects.rds` | Cached denoising results — reloaded automatically if merging fails |
| `{prefix}_dropped_samples.txt` | Samples with zero reads after filtering *(if any)* |
| `{prefix}_merge_failed_samples.txt` | Samples that failed merging *(if any)* |

### From `vc_hsp60_classify.R`

| File | Description |
|------|-------------|
| `{prefix}_taxonomy_table.csv` | Taxonomy split into rank columns (Kingdom–Species) for phyloseq |
| `{prefix}_taxonomy_confidence.csv` | Per-ASV confidence scores from sklearn classifier |
| `{prefix}_ASVs_counts_taxonomy.tsv` | Counts table merged with taxonomy |
| `{prefix}_blast_screen.csv` | Per-ASV BLAST screen results — Vibrio/non-Vibrio, % identity, query coverage *(vibrio method only)* |
| `{prefix}_vibrio_summary.tsv` | Per-sample Vibrio proportion statistics *(vibrio method only)* |
| `{prefix}_vibrio_species_abundance.tsv` | Per-sample Vibrio species relative abundance, long format *(vibrio method only)* |

### From `vc_hsp60_add_taxonomy.R`

| File | Description |
|------|-------------|
| `{prefix}_phyloseq_vibrio.rds` | Vibrio-only phyloseq object — use for species-level community analysis *(vibrio method only)* |
| `{prefix}_phyloseq_full.rds` | All ASVs with full taxonomy — use for Vibrio relative abundance estimates *(vibrio method only)* |
| `{prefix}_phyloseq_taxonomy.rds` | Phyloseq object with taxonomy and ASV sequences *(sklearn method only)* |

---

## Understanding the Vibrio summary table

When using the vibrio method, `{prefix}_vibrio_summary.tsv` provides per-sample Vibrio proportion statistics:

| Column | Description |
|--------|-------------|
| `sample` | Sample ID |
| `total_reads` | Total reads in the sample |
| `vibrio_reads` | Reads classified as Vibrio by the Vibrio-only classifier |
| `vibrio_universal_reads` | Reads classified as Vibrio/Aliivibrio by the universal classifier |
| `combined_vibrio_reads` | Total Vibrio reads from both classifiers |
| `non_vibrio_reads` | Classified non-Vibrio reads |
| `unassigned_reads` | Reads with no taxonomy assignment |
| `prop_vibrio` | Proportion of phyloseq reads that are Vibrio |
| `prop_vibrio_of_input` | Proportion of raw input reads that are Vibrio *(if read tracking table is present)* |
| `prop_vibrio_of_merged` | Proportion of merged reads that are Vibrio *(if read tracking table is present)* |
| `reads_in` | Raw input reads from pipeline *(if read tracking table is present)* |
| `nonchim` | Merged reads from pipeline *(if read tracking table is present)* |
| `vibrio_asv_count` | Number of Vibrio ASVs detected |
| `vibrio_species_count` | Number of unique Vibrio species detected |

The `{prefix}_vibrio_species_abundance.tsv` file provides per-sample species relative abundance in long format — one row per sample-species combination, with `count`, `relative_abundance` (within Vibrio reads), and `total_vibrio_reads` columns.

---

## Understanding the read tracking table

After a successful DADA2 run, `{prefix}_read_tracking.tsv` shows how many reads passed each step per sample:

| sample | reads_in | filtered | denoised | merged | nonchim | pct_retained |
|--------|----------|----------|----------|--------|---------|--------------|
| Sample1 | 150000 | 120000 | 118500 | 112000 | 108000 | 72.0 |

**What each column means:**
- `reads_in` — raw reads entering the pipeline
- `filtered` — reads passing quality filtering and truncation
- `denoised` — unique sequences after DADA2 denoising
- `merged` — read pairs successfully merged
- `nonchim` — reads remaining after chimera removal
- `pct_retained` — percentage of input reads in the final output

**If `pct_retained` is very low, check:**
- **Low `filtered`:** reads may be poor quality or truncation lengths too short — check `{prefix}_quality_profiles.pdf`
- **Low `merged`:** truncation lengths may be too short for reads to overlap — the Hsp60 amplicon is ~430bp, requiring F+R ≥ 450bp total
- **Low `nonchim`:** high chimera rate can indicate primer sequences remaining on reads

---

## Restarting a failed job

The pipeline saves intermediate results automatically so you don't have to start over if a job fails.

**Restart flags** — set to `1` to skip steps already completed:
```bash
SKIP_IMPORT=0      # skip QIIME2 import
SKIP_CUTADAPT=0    # skip adapter trimming
SKIP_DADA2=0       # skip DADA2
SKIP_CLASSIFY=0    # skip classification (full script only)
```

When skipping cutadapt, set `TRIMMED_QZA` to point to your trimmed reads from the previous run:
```bash
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
```

**Common restart scenarios:**

*Trimming worked but DADA2 failed:*
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
SKIP_DADA2=0
```

*Denoising completed but merging failed:*

No extra flags needed. The pipeline saves denoising results to `{prefix}_dada_objects.rds` automatically. On resubmit it will skip denoising and go straight to merging.
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
SKIP_DADA2=0
```

> To force denoising to rerun (e.g. after changing truncation lengths or pooling method):
> ```bash
> rm Hsp60_MyRun_dada_objects.rds
> ```

*DADA2 worked but classification failed:*
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
SKIP_DADA2=1
TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
SKIP_CLASSIFY=0
```

### Reusing a saved error model

The error model is cached to `{prefix}_error_model.rds` after the first run and reloaded automatically on reruns with the same prefix. To reuse an error model from a different run, set:
```bash
ERROR_MODEL_RDS=Hsp60_MyRun_error_model.rds
```

To force the pipeline to relearn the error model from scratch:
```bash
rm Hsp60_MyRun_error_model.rds
# leave ERROR_MODEL_RDS blank
```

---

## Runtime and performance

DADA2 is the most computationally intensive step. Runtimes below are approximate for UNH Premise:

| Dataset size | Pooling | Approximate runtime |
|---|---|---|
| ~30 samples, mixed read depths | `pseudo` | 1–2 days |
| ~60 samples, high read depths | `pseudo` | 2–3 days |
| ~30 samples, mixed read depths | `FALSE` | 2–6 hours |

### Pooling methods

| Method | Speed | Sensitivity | Recommended for |
|--------|-------|-------------|----------------|
| `pseudo` | Medium | Medium-high | Most datasets — default |
| `FALSE` | Fastest | Lower | High read depth, high-abundance taxa |
| `TRUE` | Slowest | Highest | Small datasets with rare taxa |

### Classification confidence threshold

The `--confidence` argument (default: `0.7`) controls how certain the classifier must be before assigning a taxonomy. The value used is printed at the start of the classification log.

- **Higher values** (e.g. `0.9`) — fewer assignments, more conservative
- **Lower values** (e.g. `0.5`) — more assignments, less certainty
- Unassigned ranks are filled with the rank prefix (e.g. `g__`, `s__`) following standard phyloseq convention

---

## Merging multiple runs or plates

If you have data from multiple sequencing runs or plates, merge them using the merge script:

```bash
Rscript scripts/merge_vc_hsp60_ASVs.R merged_run \
  results/plate1_Counts_seqASV_b.tsv results/plate1_ASVs.fa \
  results/plate2_Counts_seqASV_b.tsv results/plate2_ASVs.fa
```

> Only sequence-based ASV tables (`_Counts_seqASV_b.tsv`) are mergeable — do not merge numbered ASV tables directly.

---

## Using the phyloseq object in R

```R
library(phyloseq)

# Vibrio-only object — for species-level community analysis
ps_vibrio <- readRDS("Hsp60_MyRun_phyloseq_vibrio.rds")

# Full object — for Vibrio relative abundance within broader community
ps_full <- readRDS("Hsp60_MyRun_phyloseq_full.rds")

otu_table(ps_vibrio)    # ASV counts
sample_data(ps_vibrio)  # sample metadata
tax_table(ps_vibrio)    # taxonomy (Kingdom through Species)
refseq(ps_vibrio)       # ASV sequences
```

---

## Troubleshooting

**Classification fails with a scikit-learn version error:**

The classifier must match the scikit-learn version in your QIIME2 environment. Retrain using the pre-imported reference files:

```bash
conda activate <your_qiime2_env>

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /path/to/cpn60_v11_seqs.qza \
  --i-reference-taxonomy /path/to/cpn60_v11_taxonomy.qza \
  --o-classifier cpn60_classifier_v11_retrained.qza
```

> **UNH Premise:** Reference files are in the shared directory:
> ```
> /mnt/home/whistler/shared/cpn60-Classifier/
>   cpn60_v11_seqs.qza        # reference sequences (pre-imported)
>   cpn60_v11_taxonomy.qza    # taxonomy (pre-imported)
> ```

**Very low read retention in the `nonchim` column:**

If reads pass filtering and merging but are almost entirely removed as chimeras, check whether gene-specific primer sequences are still present on the reads. Primer sequences remaining on reads will cause them to appear chimeric. Use `cutadapt` with `--front-f` / `--front-r` to trim 5' primers before running DADA2.

**Samples dropped with zero reads after filtering:**

The pipeline automatically excludes these samples and writes their names to `{prefix}_dropped_samples.txt`. This is usually caused by very low input read counts, poor read quality, or truncation lengths that are too long. Check `{prefix}_quality_profiles.pdf`.

**Denoising completed but merging fails with "Non-corresponding objects" error:**

This means the cached `{prefix}_dada_objects.rds` was saved with different filtered reads than what exists now (e.g. you changed truncation lengths). Delete the cache and rerun:
```bash
rm Hsp60_MyRun_dada_objects.rds
```

**QIIME2 not found in PATH:**
```bash
# Activate your QIIME2 environment before running classify:
conda activate <your_qiime2_env>
# UNH Premise:
module purge && module load anaconda/colsa && source activate qiime2-2020.2
```

---

## Arguments reference

### vc_hsp60_pipeline.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--reads_path` | Yes | — | Path to folder containing raw paired FASTQ files |
| `--output_prefix` | Yes | — | Prefix for all output files |
| `--error_model` | No | `loess` | Error model: `loess` (recommended), `default`, or `compare` |
| `--error_model_rds` | No | — | Path to a saved error model `.rds` to skip error learning |
| `--pool` | No | `pseudo` | DADA2 pooling: `pseudo`, `FALSE`, or `TRUE` |

### vc_hsp60_classify.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_prefix` | Yes | — | Must match prefix used in `vc_hsp60_pipeline.R` |
| `--method` | No | `vibrio` | Classification method: `vibrio` (recommended) or `sklearn` |
| `--blast_db` | vibrio | — | Path to Vibrio hsp60 reference FASTA (`repset_final_130219.fas`) |
| `--vibrio_classifier` | vibrio | — | Path to Vibrio-only sklearn classifier `.qza` |
| `--classifier` | No | — | Path to universal cpn60 classifier `.qza` (used for non-Vibrio ASVs in vibrio method) |
| `--qiime2_env` | No | `qiime2-2020.2` | Name of QIIME2 conda environment |
| `--confidence` | No | `0.7` | Taxonomy confidence threshold (0–1) |
| `--phyloseq_rds` | No | `{prefix}_phyloseq.rds` | Override phyloseq file path |

### vc_hsp60_add_taxonomy.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_prefix` | Yes | — | Must match prefix used in previous steps |
| `--method` | No | `vibrio` | Must match method used in `vc_hsp60_classify.R` |
| `--phyloseq_rds` | No | `{prefix}_phyloseq.rds` | Override phyloseq file path |

---

## Repository structure

```
vc_hsp60_pipeline/
  scripts/
    vc_hsp60_pipeline.R        # DADA2 pipeline
    vc_hsp60_classify.R        # QIIME2 taxonomy classification
    vc_hsp60_add_taxonomy.R    # Merges taxonomy + sequences into phyloseq RDS
    merge_vc_hsp60_ASVs.R      # Multi-run merge utility
    install_packages.R         # Bioconductor installer (run once)
  run_vchsp60_full.slurm       # SLURM: full pipeline
  run_vchsp60_dada2.slurm      # SLURM: trimming + DADA2 only
  run_vchsp60_classify.slurm   # SLURM: classification only
  environment.yml              # conda environment
  README.md
  CITATION.cff
  LICENSE
```

---

## Citation

If you use this pipeline in your research, please cite:

> Foxall, R., Whistler, C. and Jones, S. (2026). *vc_hsp60_pipeline: A Vibrio-centric Hsp60 amplicon sequencing pipeline*. Whistler Lab, University of New Hampshire. https://github.com/ranfoxall/vc_hsp60_pipeline

A `CITATION.cff` file is included — GitHub will display a **"Cite this repository"** button automatically.

Please also cite the original assay:

> King WL, Siboni N, Kahlke T, Green TJ, Labbate M and Seymour JR (2019). A New High Throughput Sequencing Assay for Characterizing the Diversity of Natural *Vibrio* Communities and Its Application to a Pacific Oyster Mortality Event. *Frontiers in Microbiology* 10:2907. https://doi.org/10.3389/fmicb.2019.02907

### Tool citations

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA and Holmes SP (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581–583. https://doi.org/10.1038/nmeth.3869

> Bolyen E, Rideout JR, Dillon MR et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology* 37:852–857. https://doi.org/10.1038/s41587-019-0209-9

> Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17(1):10–12. https://doi.org/10.14806/ej.17.1.200

> McMurdie PJ and Holmes S (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8(4):e61217. https://doi.org/10.1371/journal.pone.0061217

> Camacho C, Coulouris G, Avagyan V et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. https://doi.org/10.1186/1471-2105-10-421

> R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

---

## Funding and acknowledgements

This pipeline was developed by Randi Foxall in the Whistler Lab at the University of New Hampshire under the supervision of Dr. Cheryl Whistler, with support from Dr. Stephen Jones.

Funding provided by the New Hampshire Agricultural Experiment Station through the NHAES CREATE program, and by NH EPSCoR.

---

## Authors

| Name | Role | Affiliation |
|------|------|-------------|
| Randi Foxall | Developer | Whistler Lab, University of New Hampshire |
| Cheryl Whistler | Principal Investigator & Funder | Dept. of Molecular, Cellular and Biomedical Sciences, UNH |
| Stephen Jones | Co-Investigator & Funder | Dept. of Natural Resources and the Environment, UNH |

---

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt this material for non-commercial purposes with appropriate credit. Commercial use requires written permission from the authors.