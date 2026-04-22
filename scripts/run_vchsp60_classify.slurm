#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=32000
#SBATCH -p shared
#SBATCH --job-name="vchsp60_classify"
#SBATCH --output=vchsp60_classify_%j.log
#SBATCH --error=vchsp60_classify_%j.err
#SBATCH --exclude=node117,node118

# run_vchsp60_classify.slurm
# Vibrio-centric hsp60 amplicon pipeline — taxonomy classification only
# Run this after run_vchsp60_dada2.slurm completes successfully.
#
# BEFORE SUBMITTING: edit the variables below, then run:
#   sbatch run_vchsp60_classify.slurm

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

# Must match the OUTPUT_PREFIX used in run_vchsp60_dada2.slurm
OUTPUT_PREFIX=Hsp60_MyRun

# Working directory — same as used in run_vchsp60_dada2.slurm
WORKDIR=/path/to/your/workdir

# Path to the cloned vc_hsp60_pipeline repo (contains the scripts/ folder)
SCRIPT_PATH=/path/to/vc_hsp60_pipeline

# Classification method: vibrio (default) or sklearn
METHOD=vibrio

# For vibrio method — BLAST reference and Vibrio-only classifier
BLAST_DB=/path/to/repset_final_130219.fas
VIBRIO_CLASSIFIER=/path/to/vc_hsp60_classifier_sklearn022.qza

# Taxonomy confidence threshold (0-1, default 0.7)
CONFIDENCE=0.7

# =============================================================================
# do not edit below this line
# =============================================================================

cd $WORKDIR

# [step 1] classify -----------------------------------------------------------
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

if [ $? -ne 0 ]; then
  echo "ERROR: Classification failed."
  exit 1
fi

# [step 2] merge taxonomy + refseq into phyloseq RDS --------------------------
module purge
module load anaconda/colsa
source activate $PIPELINE_ENV

srun Rscript $SCRIPT_PATH/scripts/vc_hsp60_add_taxonomy.R \
  --output_prefix $OUTPUT_PREFIX

if [ $? -ne 0 ]; then
  echo "ERROR: Adding taxonomy to phyloseq failed."
  exit 1
fi

echo "Done."
