#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=512000
#SBATCH -p shared
#SBATCH --job-name="vchsp60_dada2"
#SBATCH --output=vchsp60_dada2_%j.log
#SBATCH --error=vchsp60_dada2_%j.err
#SBATCH --exclude=node117,node118

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

UNH_PREMISE=1

if [ $UNH_PREMISE -eq 1 ]; then
  QIIME2_ENV=qiime2-2020.2
  PIPELINE_ENV=vc_hsp60_pipeline
else
  QIIME2_ENV=qiime2-2024.10.1                    # name of your QIIME2 conda environment
  PIPELINE_ENV=vc_hsp60_pipeline                 # name of your pipeline conda environment
fi

# ---------------- EDIT THESE ----------------

READS_PATH=/path/to/your/raw_reads
OUTPUT_PREFIX=Hsp60_MyRun
WORKDIR=/path/to/your/workdir
SCRIPT_PATH=/path/to/vc_hsp60_pipeline

ERROR_MODEL=loess
POOL=pseudo

# Bases to trim from 5' end of reads to remove gene-specific primers
# Set to 0 if primers are not present on your reads
# For Vib-hspF1/R400 primers: TRIM_LEFT_F=23, TRIM_LEFT_R=24
TRIM_LEFT_F=0
TRIM_LEFT_R=0

ERROR_MODEL_RDS=

SKIP_IMPORT=0
SKIP_CUTADAPT=0
SKIP_DADA2=0

TRIMMED_QZA=

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

  if [ $? -ne 0 ]; then echo "ERROR: FASTQ extraction failed."; exit 1; fi

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi

else
  if [ -z "$TRIMMED_QZA" ]; then
    echo "ERROR: SKIP_CUTADAPT=1 but TRIMMED_QZA is not set."
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

echo "Done. Submit run_vchsp60_classify.slurm to run taxonomy classification."