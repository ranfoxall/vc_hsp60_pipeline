#!/usr/bin/env Rscript

# merge_vc_hsp60_ASVs.R
# Merges ASV count tables from multiple hsp60 sequencing runs.
#
# Merging is done by exact ASV sequence — identical sequences across runs
# are treated as the same ASV and their counts are summed. New sequential
# ASV numbers are assigned to the merged set. Outputs are named to match
# the conventions expected by vc_hsp60_classify.R so classification can
# be run directly on the merged dataset.
#
# Usage:
#   Rscript merge_vc_hsp60_ASVs.R <output_prefix> \
#     <run1_Counts_seqASV_b.tsv> <run1_ASVs.fa> \
#     <run2_Counts_seqASV_b.tsv> <run2_ASVs.fa> \
#     [<run3_Counts_seqASV_b.tsv> <run3_ASVs.fa> ...]
#
# Example:
#   Rscript merge_vc_hsp60_ASVs.R vcHsp60_merged \
#     vcHsp60_old_redo_Counts_seqASV_b.tsv vcHsp60_old_redo_ASVs.fa \
#     vcHsp60_MyRun4_Counts_seqASV_b.tsv vcHsp60_MyRun4_ASVs.fa
#
# Outputs:
#   {prefix}_Counts_seqASV_b.tsv  — merged counts, ASV sequences as row names
#   {prefix}_Counts_numASV.tsv    — merged counts, ASV_1/ASV_2/... as row names
#   {prefix}_ASVs.fa              — merged FASTA with ASV_# headers
#   {prefix}_phyloseq.rds         — phyloseq object (counts + sample metadata)
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Funder:     Dr. Cheryl Whistler, Dr. Stephen Jones
# Lab:        Whistler Lab, Dept. of Molecular, Cellular and Biomedical Sciences, UNH
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3 || length(args) %% 2 != 1) {
  stop(
    "Usage: Rscript merge_vc_hsp60_ASVs.R <output_prefix> ",
    "<counts1.tsv> <fasta1.fa> [<counts2.tsv> <fasta2.fa> ...]\n",
    "Counts files must be sequence-based ASV tables (*_Counts_seqASV_b.tsv)"
  )
}

prefix    <- args[1]
file_args <- args[-1]

counts_files <- file_args[seq(1, length(file_args), by=2)]
fasta_files  <- file_args[seq(2, length(file_args), by=2)]

cat("Merging", length(counts_files), "runs into prefix:", prefix, "\n")
for (i in seq_along(counts_files)) {
  cat("  Run", i, "—", counts_files[i], "/", fasta_files[i], "\n")
}

# [1] read and validate input files ------------------------------------
all_counts <- list()

for (i in seq_along(counts_files)) {
  if (!file.exists(counts_files[i])) stop("Counts file not found: ", counts_files[i])
  if (!file.exists(fasta_files[i]))  stop("FASTA file not found: ", fasta_files[i])

  counts <- read.table(counts_files[i], header=TRUE, sep="\t",
                       row.names=1, check.names=FALSE)
  cat("Run", i, "—", nrow(counts), "ASVs,", ncol(counts), "samples\n")

  # validate that row names look like sequences (not ASV_1 style)
  if (any(grepl("^ASV_", rownames(counts)))) {
    stop(
      "Counts file ", counts_files[i], " appears to use numbered ASV IDs.\n",
      "Please use the sequence-based counts table (*_Counts_seqASV_b.tsv)."
    )
  }

  all_counts[[i]] <- counts
}

# [2] merge by sequence ------------------------------------------------
# convert to long format with sequence as key, then full_join
all_counts_named <- lapply(all_counts, function(df) {
  tibble::rownames_to_column(as.data.frame(df), "ASV_sequence")
})

cat("Merging by ASV sequence...\n")
merged_counts <- Reduce(function(x, y) {
  full_join(x, y, by="ASV_sequence")
}, all_counts_named)

merged_counts[is.na(merged_counts)] <- 0

# extract sequences before setting rownames
merged_seqs                <- merged_counts$ASV_sequence
rownames(merged_counts)    <- merged_seqs
merged_counts$ASV_sequence <- NULL

cat("Merged table:", nrow(merged_counts), "unique ASVs,",
    ncol(merged_counts), "total samples\n")

# check for duplicate sample names across runs
if (any(duplicated(colnames(merged_counts)))) {
  dupes <- colnames(merged_counts)[duplicated(colnames(merged_counts))]
  stop(
    "Duplicate sample names found across runs: ",
    paste(dupes, collapse=", "), "\n",
    "Please rename samples in your input files to be unique."
  )
}

# [3] write sequence-based outputs ------------------------------------
seq_counts_file <- paste0(prefix, "_Counts_seqASV_b.tsv")
write.table(merged_counts, seq_counts_file, sep="\t",
            quote=FALSE, col.names=NA, row.names=TRUE)
cat("Sequence-based counts written to:", seq_counts_file, "\n")

# [4] assign new sequential ASV numbers --------------------------------
asv_ids              <- paste0("ASV_", seq_len(nrow(merged_counts)))
rownames(merged_counts) <- asv_ids

# [5] write numbered outputs -------------------------------------------
num_counts_file <- paste0(prefix, "_Counts_numASV.tsv")
write.table(merged_counts, num_counts_file, sep="\t",
            quote=FALSE, col.names=NA, row.names=TRUE)
cat("Numbered counts written to:", num_counts_file, "\n")

# [6] write FASTA ------------------------------------------------------
asv_fasta_path <- paste0(prefix, "_ASVs.fa")
fasta_lines    <- unlist(lapply(seq_along(asv_ids), function(i) {
  c(paste0(">", asv_ids[i]), merged_seqs[i])
}))
writeLines(fasta_lines, asv_fasta_path)
cat("FASTA written to:", asv_fasta_path, "\n")

# [7] phyloseq object --------------------------------------------------
if (requireNamespace("phyloseq", quietly=TRUE)) {
  library(phyloseq)
  counts_matrix <- as.matrix(merged_counts)
  mode(counts_matrix) <- "numeric"

  sample_metadata <- data.frame(
    SampleID   = colnames(counts_matrix),
    TotalReads = colSums(counts_matrix),
    row.names  = colnames(counts_matrix),
    stringsAsFactors = FALSE
  )

  ps <- phyloseq(
    otu_table(counts_matrix, taxa_are_rows=TRUE),
    sample_data(sample_metadata)
  )

  ps_file <- paste0(prefix, "_phyloseq.rds")
  saveRDS(ps, ps_file)
  cat("Phyloseq RDS written to:", ps_file, "\n")
} else {
  cat("phyloseq not installed — skipping RDS creation.\n")
}

# [8] summary ----------------------------------------------------------
cat("\nMerge complete.\n")
cat("Total unique ASVs:", nrow(merged_counts), "\n")
cat("Total samples:    ", ncol(merged_counts), "\n")
cat("\nOutputs:\n")
cat("  ", seq_counts_file, "  (use for future merges)\n")
cat("  ", num_counts_file, "\n")
cat("  ", asv_fasta_path, "\n")
if (exists("ps_file")) cat("  ", ps_file, "\n")
cat("\nNext step: run vc_hsp60_classify.R with --output_prefix", prefix, "\n")