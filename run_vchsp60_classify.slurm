#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=32000
#SBATCH -p shared
#SBATCH --job-name="vchsp60_classify"
#SBATCH --output=vchsp60_classify_%j.log
#SBATCH --error=vchsp60_classify_%j.err
#SBATCH --exclude=node117,node118

# =============================================================================
# run_vchsp60_classify.slurm
# Vibrio-centric hsp60 amplicon pipeline — taxonomy classification only
# Run this after run_vchsp60_dada2.slurm completes successfully.
#
# BEFORE SUBMITTING: edit the variables below, then run:
#   sbatch run_vchsp60_classify.slurm
# =============================================================================

# [system] --------------------------------------------------------------------
# Set to 1 if running on UNH Premise — environments and classifier set automatically.
# Set to 0 for all other systems and fill in the variables in the else block.
UNH_PREMISE=1

if [ $UNH_PREMISE -eq 1 ]; then
  QIIME2_ENV=qiime2-2020.2
  CLASSIFIER=/mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
else
  QIIME2_ENV=qiime2-2024.10.1                    # name of your QIIME2 conda environment
  CLASSIFIER=/path/to/cpn60_classifier_v11.qza   # path to downloaded classifier .qza
fi

# [paths] — edit these before submitting --------------------------------------

# Must match the OUTPUT_PREFIX used in run_vchsp60_dada2.slurm
OUTPUT_PREFIX=Hsp60_MyRun

# Working directory — same as used in run_vchsp60_dada2.slurm
# UNH Premise example: /mnt/home/whistler/yourname/hsp60_vc
WORKDIR=/path/to/your/workdir

# Path to the cloned vc_hsp60_pipeline repo (contains the scripts/ folder)
# UNH Premise example: /mnt/home/whistler/yourname/vc_hsp60_pipeline
SCRIPT_PATH=/path/to/vc_hsp60_pipeline

# =============================================================================
# do not edit below this line
# =============================================================================

cd $WORKDIR

module purge
module load anaconda/colsa
source activate $QIIME2_ENV

srun Rscript $SCRIPT_PATH/scripts/vc_hsp60_classify.R \
  --output_prefix $OUTPUT_PREFIX \
  --classifier $CLASSIFIER

if [ $? -ne 0 ]; then
  echo "ERROR: Classification failed."
  exit 1
fi

echo "Done."
