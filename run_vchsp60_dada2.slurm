#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=512000
#SBATCH -p shared
#SBATCH --job-name="vchsp60_dada2"
#SBATCH --output=vchsp60_dada2_%j.log
#SBATCH --error=vchsp60_dada2_%j.err
#SBATCH --exclude=node117,node118

# =============================================================================
# run_vchsp60_dada2.slurm
# Vibrio-centric hsp60 amplicon pipeline — trimming + DADA2 only
# Run this first, then submit run_vchsp60_classify.slurm when it completes.
#
# BEFORE SUBMITTING: edit the variables below, then run:
#   sbatch run_vchsp60_dada2.slurm
# =============================================================================

# [system] --------------------------------------------------------------------
# Set to 1 if running on UNH Premise — environments set automatically.
# Set to 0 for all other systems and fill in the variables in the else block.
UNH_PREMISE=1

if [ $UNH_PREMISE -eq 1 ]; then
  QIIME2_ENV=qiime2-2020.2
  PIPELINE_ENV=vc_hsp60_pipeline
else
  QIIME2_ENV=qiime2-2024.10.1                    # name of your QIIME2 conda environment
  PIPELINE_ENV=vc_hsp60_pipeline                 # name of your DADA2 conda environment
fi

# [paths] — edit these before submitting --------------------------------------

# Path to your raw FASTQ folder (must contain .fastq.gz files directly)
READS_PATH=/path/to/your/raw_reads

# Short name for this run — all output files will start with this
OUTPUT_PREFIX=Hsp60_MyRun

# Working directory — where outputs will be written
# UNH Premise example: /mnt/home/whistler/yourname/hsp60_vc
WORKDIR=/path/to/your/workdir

# Path to the cloned vc_hsp60_pipeline repo (contains the scripts/ folder)
# UNH Premise example: /mnt/home/whistler/yourname/vc_hsp60_pipeline
SCRIPT_PATH=/path/to/vc_hsp60_pipeline

# DADA2 error model: loess (recommended), default, or compare (runs both, picks best — slower)
ERROR_MODEL=loess

# path to a previously saved error model .rds — leave blank to learn from scratch
# e.g. ERROR_MODEL_RDS=Hsp60_MyRun_error_model.rds
ERROR_MODEL_RDS=

# [restart flags] — set to 1 to skip steps already completed in a previous run
SKIP_IMPORT=0
SKIP_CUTADAPT=0
SKIP_DADA2=0

# Required when SKIP_CUTADAPT=1: path to trimmed reads .qza from previous run
# e.g. TRIMMED_QZA=Hsp60_MyRun_trimmed-reads.qza
TRIMMED_QZA=

# =============================================================================
# do not edit below this line
# =============================================================================

cd $WORKDIR

# [step 1] import -------------------------------------------------------------
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

# [step 2] trim Illumina adapters ---------------------------------------------
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
    --o-trimmed-sequences ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --verbose

  if [ $? -ne 0 ]; then echo "ERROR: Illumina adapter trimming failed."; exit 1; fi

  srun qiime demux summarize \
    --i-data ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --o-visualization ${OUTPUT_PREFIX}_trimmed-reads.qzv

  rm -rf ${OUTPUT_PREFIX}_Exported_Reads

  srun qiime tools extract \
    --input-path ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --output-path ${OUTPUT_PREFIX}_Exported_Reads

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  echo "Trimmed reads exported to: $TRIMMED_PATH"
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi

else
  if [ -z "$TRIMMED_QZA" ]; then
    echo "ERROR: SKIP_CUTADAPT=1 but TRIMMED_QZA is not set."
    echo "Set TRIMMED_QZA to your trimmed reads .qza (e.g. Hsp60_MyRun_trimmed-reads.qza)"
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

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  echo "Trimmed reads exported to: $TRIMMED_PATH"
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi
fi

# [step 3] DADA2 --------------------------------------------------------------
if [ $SKIP_DADA2 -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $PIPELINE_ENV

  srun Rscript $SCRIPT_PATH/scripts/vc_hsp60_pipeline.R \
    --reads_path $TRIMMED_PATH \
    --output_prefix $OUTPUT_PREFIX \
    --error_model $ERROR_MODEL \
    ${ERROR_MODEL_RDS:+--error_model_rds $ERROR_MODEL_RDS}

  if [ $? -ne 0 ]; then echo "ERROR: DADA2 pipeline failed."; exit 1; fi
else
  echo "Skipping DADA2."
fi

echo "Done. Submit run_vchsp60_classify.slurm to run taxonomy classification."
