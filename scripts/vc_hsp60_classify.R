# vc_hsp60_classify.R
# QIIME2 taxonomy classification for Vibrio-centric hsp60 ASVs
#
# Run AFTER vc_hsp60_pipeline.R with your QIIME2 environment active.
# On UNH Premise, use the SLURM script (run_vchsp60_classify.slurm).
# Manually:
#   conda activate qiime2-2020.2
#   Rscript vc_hsp60_classify.R \
#     --output_prefix Hsp60_MyRun \
#     --classifier    /path/to/cpn60_classifier_v11.qza \
#     --confidence    0.7
#
# UNH Premise classifier path:
#   /mnt/home/whistler/shared/cpn60-Classifier/cpn60_classifier_v11_sklearn142.qza
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Funder:     Dr. Cheryl Whistler, Dr. Stephen Jones
# Lab:        Whistler Lab, Dept. of Molecular, Cellular and Biomedical Sciences, UNH
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0

# [1] arguments -------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

parse_args_simple <- function(args) {
  result <- list(output_prefix=NULL, classifier=NULL, phyloseq_rds=NULL, confidence="0.7")
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--output_prefix" && i < length(args)) {
      result$output_prefix <- args[i+1]; i <- i + 2
    } else if (args[i] == "--classifier" && i < length(args)) {
      result$classifier <- args[i+1]; i <- i + 2
    } else if (args[i] == "--phyloseq_rds" && i < length(args)) {
      result$phyloseq_rds <- args[i+1]; i <- i + 2
    } else if (args[i] == "--confidence" && i < length(args)) {
      result$confidence <- args[i+1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  result
}

opt <- parse_args_simple(args)

for (arg in c("output_prefix", "classifier")) {
  if (is.null(opt[[arg]])) stop("Please provide --", arg)
}

output_prefix    <- opt$output_prefix
cpn60_classifier <- opt$classifier
phyloseq_rds     <- if (!is.null(opt$phyloseq_rds)) opt$phyloseq_rds else paste0(output_prefix, "_phyloseq.rds")
confidence       <- opt$confidence

cat("output_prefix: ", output_prefix, "\n")
cat("classifier:    ", cpn60_classifier, "\n")
cat("confidence:    ", confidence, "\n")

# [2] check QIIME2 ----------------------------------------------------
qiime_check <- system2("which", args="qiime", stdout=TRUE, stderr=TRUE)
if (length(qiime_check) == 0 || !grepl("qiime", qiime_check[1])) {
  stop(
    "QIIME2 not found in PATH.\n",
    "Activate your QIIME2 conda environment before running:\n",
    "  conda activate qiime2-2020.2"
  )
}
cat("QIIME2 found at:", qiime_check[1], "\n")

# [3] classify --------------------------------------------------------
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
  "--p-confidence",     confidence,
  "--o-classification", taxonomy_qza
))
if (ret != 0) stop(
  "Classification failed.\n",
  "Common cause: scikit-learn version mismatch.\n",
  "On UNH Premise, use cpn60_classifier_v11_sklearn142.qza.\n",
  "If using a different system, retrain the classifier to match your QIIME2's sklearn version — see README."
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

# [4] build taxonomy table --------------------------------------------
cat("Building taxonomy table...\n")

taxonomy_raw <- read.table(taxonomy_tsv, sep="\t", header=TRUE,
                            stringsAsFactors=FALSE, row.names=1)

tax_ranks   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_prefix <- c("k__",     "p__",    "c__",   "o__",   "f__",    "g__",   "s__")

tax_split <- strsplit(taxonomy_raw$Taxon, ";\\s*")
tax_mat <- t(sapply(tax_split, function(x) {
  length(x) <- length(tax_ranks)
  x
}))
rownames(tax_mat) <- rownames(taxonomy_raw)
colnames(tax_mat) <- tax_ranks

for (i in seq_along(tax_ranks)) {
  empty <- is.na(tax_mat[, i]) | tax_mat[, i] == "" | tax_mat[, i] == "-" |
           tax_mat[, i] == rank_prefix[i]
  tax_mat[empty, i] <- rank_prefix[i]
}

taxonomy_table      <- as.data.frame(tax_mat, stringsAsFactors=FALSE)
taxonomy_table_file <- paste0(output_prefix, "_taxonomy_table.csv")
write.csv(taxonomy_table, taxonomy_table_file, quote=FALSE)
cat("Taxonomy table written to:", taxonomy_table_file, "\n")

confidence_file <- paste0(output_prefix, "_taxonomy_confidence.csv")
write.csv(data.frame(Confidence=taxonomy_raw$Confidence, row.names=rownames(taxonomy_raw)),
          confidence_file, quote=FALSE)
cat("Confidence scores written to:", confidence_file, "\n")

counts_num   <- read.table(counts_table, sep="\t", header=TRUE,
                           stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
taxonomy_raw <- taxonomy_raw[rownames(counts_num), , drop=FALSE]
counts_tax      <- cbind(counts_num, taxonomy_raw)
counts_tax_file <- paste0(output_prefix, "_ASVs_counts_taxonomy.tsv")
write.table(counts_tax, counts_tax_file, sep="\t", quote=FALSE, col.names=NA)
cat("Counts + taxonomy written to:", counts_tax_file, "\n")

# [5] summary ---------------------------------------------------------
cat("\nClassification complete. Outputs:\n")
cat("  ", taxonomy_table_file, "\n")
cat("  ", confidence_file, "\n")
cat("  ", counts_tax_file, "\n")
cat("\nNext step: vc_hsp60_add_taxonomy.R will merge taxonomy and refseq into phyloseq RDS.\n")