#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=512000
#SBATCH -p shared
#SBATCH --job-name="vchsp60_full"
#SBATCH --output=vchsp60_full_%j.log
#SBATCH --error=vchsp60_full_%j.err
#SBATCH --exclude=node117,node118

# =============================================================================
# run_vchsp60_full.slurm
# Vibrio-centric hsp60 amplicon pipeline — full run (trimming + DADA2 + classify)
#
# BEFORE SUBMITTING: edit the variables below, then run:
#   sbatch run_vchsp60_full.slurm
# =============================================================================

UNH_PREMISE=1

if [ $UNH_PREMISE -eq 1 ]; then
  QIIME2_ENV=qiime2-2020.2
  PIPELINE_ENV=vc_hsp60_pipeline
  CLASSIFIER=/mnt/home/whistler/shared/cpn60-Classifier/cpn60_classifier_v11_sklearn142.qza
else
  QIIME2_ENV=qiime2-2024.10.1                    # name of your QIIME2 conda environment
  PIPELINE_ENV=vc_hsp60_pipeline                 # name of your pipeline conda environment
  CLASSIFIER=/path/to/cpn60_classifier_v11.qza   # path to downloaded classifier .qza
fi


# ---------------- EDIT THESE ----------------

# Path to your raw FASTQ folder
READS_PATH=/path/to/your/raw_reads

# Short name for this run — all output files start with this
OUTPUT_PREFIX=Hsp60_MyRun

# Working directory — where outputs will be written
WORKDIR=/path/to/your/workdir

# Path to the cloned vc_hsp60_pipeline repo (contains the scripts/ folder)
SCRIPT_PATH=/path/to/vc_hsp60_pipeline

# DADA2 error model: loess (recommended), default, or compare
ERROR_MODEL=loess

# DADA2 pooling: pseudo (recommended), FALSE (fastest), TRUE (most sensitive)
POOL=pseudo

# Classification method: vibrio (default) or sklearn
# vibrio = BLAST screen + Vibrio sklearn + universal sklearn (King et al. 2019)
# sklearn = universal cpn60 sklearn classifier only
METHOD=vibrio

# For vibrio method — BLAST reference and Vibrio-only classifier
# UNH Premise: BLAST_DB=/mnt/home/whistler/foxall/hsp60_ref/vc_hsp60_ref/repset_final_130219.fas
# UNH Premise: VIBRIO_CLASSIFIER=/mnt/home/whistler/foxall/hsp60_ref/vc_hsp60_ref/vc_hsp60_classifier_sklearn022.qza
BLAST_DB=/path/to/repset_final_130219.fas
VIBRIO_CLASSIFIER=/path/to/vc_hsp60_classifier.qza

# Universal cpn60 classifier — used for non-Vibrio ASVs (vibrio method)
# or as the sole classifier (sklearn method)
CLASSIFIER=/mnt/home/whistler/shared/cpn60-Classifier/cpn60_classifier_v11_sklearn142.qza
CONFIDENCE=0.7

# Bases to trim from 5' end of reads to remove gene-specific primers
# Set to 0 if primers are not present on your reads (e.g. _orig samples)
# For Vib-hspF1/R400 primers: TRIM_LEFT_F=0, TRIM_LEFT_R=0
TRIM_LEFT_F=0
TRIM_LEFT_R=0

# Path to a previously saved error model .rds — leave blank to learn from scratch
ERROR_MODEL_RDS=

# [restart flags] — set to 1 to skip steps already completed
SKIP_IMPORT=0
SKIP_CUTADAPT=0
SKIP_DADA2=0
SKIP_CLASSIFY=0

# Required when SKIP_CUTADAPT=0
TRIMMED_QZA=

# Note: if denoising completed but merging failed, the pipeline will
# automatically reuse {OUTPUT_PREFIX}_dada_objects.rds if it exists.

# =============================================================================
# do not edit below this line
# =============================================================================

cd $WORKDIR

# ---------------- IMPORT ----------------
if [ $SKIP_IMPORT -eq 0 ] && [ $SKIP_CUTADAPT -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  srun qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $READS_PATH \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path ${OUTPUT_PREFIX}_demux-paired-reads.qza

  if [ $? -ne 0 ]; then echo "ERROR: QIIME2 import failed."; exit 1; fi
fi

# ---------------- CUTADAPT ----------------
if [ $SKIP_CUTADAPT -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  srun qiime cutadapt trim-paired \
    --i-demultiplexed-sequences ${OUTPUT_PREFIX}_demux-paired-reads.qza \
    --p-cores 24 \
    --p-adapter-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --p-adapter-r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --p-match-read-wildcards \
    --o-trimmed-sequences ${OUTPUT_PREFIX}_trimmed-reads.qza

  if [ $? -ne 0 ]; then echo "ERROR: Illumina adapter trimming failed."; exit 1; fi

  srun qiime demux summarize \
    --i-data ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --o-visualization ${OUTPUT_PREFIX}_trimmed-reads.qzv

  rm -rf ${OUTPUT_PREFIX}_Exported_Reads

  srun qiime tools extract \
    --input-path ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --output-path ${OUTPUT_PREFIX}_Exported_Reads

  if [ $? -ne 0 ]; then echo "ERROR: FASTQ extraction failed."; exit 1; fi

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi

else
  if [ -z "$TRIMMED_QZA" ]; then
    echo "ERROR: SKIP_CUTADAPT=0 but TRIMMED_QZA is not set."
    exit 1
  fi

  echo "Skipping cutadapt. Re-extracting FASTQs from: $TRIMMED_QZA"
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  rm -rf ${OUTPUT_PREFIX}_Exported_Reads

  srun qiime tools extract \
    --input-path $TRIMMED_QZA \
    --output-path ${OUTPUT_PREFIX}_Exported_Reads

  if [ $? -ne 0 ]; then echo "ERROR: FASTQ extraction failed."; exit 1; fi

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi
fi

# ---------------- DADA2 ----------------
if [ $SKIP_DADA2 -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $PIPELINE_ENV

  srun Rscript $SCRIPT_PATH/scripts/vc_hsp60_pipeline.R \
    --reads_path $TRIMMED_PATH \
    --output_prefix $OUTPUT_PREFIX \
    --error_model $ERROR_MODEL \
    --pool $POOL \
    --trim_left_f $TRIM_LEFT_F \
    --trim_left_r $TRIM_LEFT_R \
    ${ERROR_MODEL_RDS:+--error_model_rds $ERROR_MODEL_RDS}

  if [ $? -ne 0 ]; then echo "ERROR: DADA2 pipeline failed."; exit 1; fi
else
  echo "Skipping DADA2."
fi

# ---------------- CLASSIFY ----------------
if [ $SKIP_CLASSIFY -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $PIPELINE_ENV

  srun Rscript $SCRIPT_PATH/scripts/vc_hsp60_classify.R \
    --output_prefix $OUTPUT_PREFIX \
    --method $METHOD \
    --blast_db $BLAST_DB \
    --vibrio_classifier $VIBRIO_CLASSIFIER \
    --classifier $CLASSIFIER \
    --qiime2_env $QIIME2_ENV \
    --confidence $CONFIDENCE

  if [ $? -ne 0 ]; then echo "ERROR: Classification failed."; exit 1; fi

  module purge
  module load anaconda/colsa
  source activate $PIPELINE_ENV

  srun Rscript $SCRIPT_PATH/scripts/vc_hsp60_add_taxonomy.R \
    --output_prefix $OUTPUT_PREFIX

  if [ $? -ne 0 ]; then echo "ERROR: Adding taxonomy to phyloseq failed."; exit 1; fi
fi
