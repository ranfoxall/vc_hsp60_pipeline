# vc_hsp60_classify.R
# QIIME2 taxonomy classification for Vibrio-centric hsp60 ASVs
#
# Run AFTER vc_hsp60_pipeline.R with your QIIME2 environment active.
# On UNH Premise, use the SLURM script (run_vchsp60_classify.slurm).
# Manually:
#   conda activate <your_qiime2_env>
#   Rscript vc_hsp60_classify.R \
#     --output_prefix Hsp60_MyRun \
#     --classifier    /path/to/cpn60_classifier_v11.qza
#
# UNH Premise classifier path:
#   /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Funder:     Dr. Cheryl Whistler, Dr. Stephen Jones
# Lab:        Whistler Lab, Dept. of Molecular, Cellular and Biomedical Sciences, UNH
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
})

# [1] arguments -------------------------------------------------------
# parse and validate command line inputs
option_list <- list(
  make_option(c("--output_prefix"), type="character",
              help="Output prefix from vc_hsp60_pipeline.R — input files derived from this"),
  make_option(c("--classifier"),    type="character",
              help="Path to cpn60 QIIME2 classifier .qza"),
  make_option(c("--phyloseq_rds"),  type="character", default=NULL,
              help="(Optional) Path to phyloseq .rds — derived from prefix if not provided")
)
opt <- parse_args(OptionParser(option_list=option_list))

for (arg in c("output_prefix", "classifier")) {
  if (is.null(opt[[arg]])) stop("Please provide --", arg)
}

output_prefix    <- opt$output_prefix
cpn60_classifier <- opt$classifier
phyloseq_rds     <- if (!is.null(opt$phyloseq_rds)) opt$phyloseq_rds else paste0(output_prefix, "_phyloseq.rds")

asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
counts_table   <- paste0(output_prefix, "_Counts_numASV.tsv")

for (f in c(asv_fasta_path, counts_table, cpn60_classifier)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

cat("output_prefix: ", output_prefix, "\n")
cat("classifier:    ", cpn60_classifier, "\n")

# [2] check QIIME2 ----------------------------------------------------
# verify QIIME2 is available in PATH before attempting classification
# fail early with a clear message if QIIME2 is not active
qiime_check <- system2("which", args="qiime", stdout=TRUE, stderr=TRUE)
if (length(qiime_check) == 0 || !grepl("qiime", qiime_check[1])) {
  stop(
    "QIIME2 not found in PATH.\n",
    "Activate your QIIME2 conda environment before running:\n",
    "  conda activate <your_qiime2_env>\n",
    "UNH Premise users: conda activate qiime2-2020.2"
  )
}
cat("QIIME2 found at:", qiime_check[1], "\n")

# [3] classify --------------------------------------------------------
# import ASV FASTA to QIIME2, run cpn60 sklearn classifier, export taxonomy TSV
queries_qza         <- paste0(output_prefix, "_ASVs.qza")
taxonomy_qza        <- paste0(output_prefix, "_taxonomy.qza")
taxonomy_tsv        <- paste0(output_prefix, "_taxonomy.tsv")
taxonomy_export_dir <- paste0(output_prefix, "_taxonomy_export")

cat("Importing ASVs to QIIME2...\n")
ret <- system2("qiime", args=c(
  "tools", "import",
  "--type", "FeatureData[Sequence]",
  "--input-path",  asv_fasta_path,
  "--output-path", queries_qza
))
if (ret != 0) stop("QIIME2 import failed. Check that your QIIME2 environment is activated.")

cat("Classifying ASVs...\n")
ret <- system2("qiime", args=c(
  "feature-classifier", "classify-sklearn",
  "--i-classifier",     cpn60_classifier,
  "--i-reads",          queries_qza,
  "--o-classification", taxonomy_qza
))
if (ret != 0) stop(
  "Classification failed.\n",
  "Common cause: scikit-learn version mismatch.\n",
  "cpn60_classifier_v11.qza was trained with scikit-learn 0.24.1 (QIIME2 v2022.11).\n",
  "If your QIIME2 is newer, retrain the classifier — see README for instructions."
)

cat("Exporting taxonomy...\n")
if (!dir.exists(taxonomy_export_dir)) dir.create(taxonomy_export_dir)
system2("qiime", args=c(
  "tools", "export",
  "--input-path",  taxonomy_qza,
  "--output-path", taxonomy_export_dir
))
file.rename(file.path(taxonomy_export_dir, "taxonomy.tsv"), taxonomy_tsv)
cat("Taxonomy written to:", taxonomy_tsv, "\n")

# [4] merge counts + taxonomy -----------------------------------------
# join taxonomy assignments to numbered ASV count table by ASV ID
cat("Merging counts and taxonomy...\n")

counts_num <- read.table(counts_table, sep="\t", header=TRUE,
                         stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
taxonomy   <- read.table(taxonomy_tsv, sep="\t", header=TRUE,
                         stringsAsFactors=FALSE, row.names=1)
taxonomy   <- taxonomy[rownames(counts_num), , drop=FALSE]

counts_tax      <- cbind(counts_num, taxonomy)
counts_tax_file <- paste0(output_prefix, "_ASVs_counts_taxonomy.tsv")
write.table(counts_tax, counts_tax_file, sep="\t", quote=FALSE, col.names=NA)
cat("Counts + taxonomy written to:", counts_tax_file, "\n")

# [5] update phyloseq -------------------------------------------------
# add taxonomy to the phyloseq object produced by vc_hsp60_pipeline.R
if (!file.exists(phyloseq_rds)) {
  warning("Phyloseq RDS not found: ", phyloseq_rds, " -- skipping.")
} else {
  cat("Adding taxonomy to phyloseq object...\n")
  ps          <- readRDS(phyloseq_rds)
  tax_mat     <- as.matrix(taxonomy)
  ps          <- merge_phyloseq(ps, tax_table(tax_mat))
  ps_tax_file <- paste0(output_prefix, "_phyloseq_taxonomy.rds")
  saveRDS(ps, ps_tax_file)
  cat("Phyloseq with taxonomy written to:", ps_tax_file, "\n")
}

# [6] summary ---------------------------------------------------------
# print output file paths
cat("\nClassification complete. Outputs:\n")
cat("  ", taxonomy_tsv, "\n")
cat("  ", counts_tax_file, "\n")
if (file.exists(paste0(output_prefix, "_phyloseq_taxonomy.rds")))
  cat("  ", paste0(output_prefix, "_phyloseq_taxonomy.rds"), "\n")
