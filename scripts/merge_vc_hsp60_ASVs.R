#!/usr/bin/env Rscript

# ----------------------------
# vc_hsp60 ASV Merge Script (fixed)
# ----------------------------
# Usage:
# Rscript merge_vc_hsp60_ASVs.R <output_prefix> <counts_and_fasta_files...>
#
# Example:
# Rscript merge_vc_hsp60_ASVs.R merged_test \
#   test_plate1_Counts_seqASV_b.tsv test_plate1_ASVs.fa \
#   test_plate2_Counts_seqASV_b.tsv test_plate2_ASVs.fa
# ----------------------------

library(dplyr)
library(readr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3 | length(args) %% 2 != 1){
  stop("Usage: Rscript merge_vc_hsp60_ASVs.R <output_prefix> <counts_and_fasta_files...>")
}

prefix <- args[1]
file_args <- args[-1]

# ----------------------------
# Separate counts and fasta paths
# ----------------------------
counts_files <- file_args[seq(1, length(file_args), by = 2)]
fasta_files  <- file_args[seq(2, length(file_args), by = 2)]

all_counts <- list()
all_seqs   <- list()

# ----------------------------
# Read each plate/run
# ----------------------------
for(i in seq_along(counts_files)){
  counts <- read.table(counts_files[i], header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  all_counts[[i]] <- counts

  seqs <- readLines(fasta_files[i])
  # Extract sequences from fasta (every 2nd line)
  seq_lines <- seqs[seq(2, length(seqs), by=2)]
  all_seqs[[i]] <- seq_lines
}

# ----------------------------
# Merge counts by ASV sequence (clean headers)
# ----------------------------
all_counts_named <- lapply(seq_along(all_counts), function(i) {
  df <- as.data.frame(all_counts[[i]]) %>% tibble::rownames_to_column("ASV_sequence")
  # If columns are generic V1,V2..., rename to Sample1, Sample2...
  if(all(grepl("^V", colnames(df)[-1]))) {
    colnames(df)[-1] <- paste0("Sample", seq_len(ncol(df) - 1))
  }
  df
})

merged_counts <- Reduce(function(x, y) {
  full_join(x, y, by = "ASV_sequence")
}, all_counts_named)

merged_counts[is.na(merged_counts)] <- 0

# ASV sequences as rownames
rownames(merged_counts) <- merged_counts$ASV_sequence
merged_seqs <- rownames(merged_counts)
merged_counts$ASV_sequence <- NULL

# ----------------------------
# Write mergeable _b dataset (sequence as ID)
# ----------------------------
b_counts_file <- paste0(prefix, "_merged_ASVs_b.tsv")
b_fasta_file  <- paste0(prefix, "_merged_ASVs_b.fa")

# Counts table: rows = sequences, columns = sample names
write.table(merged_counts, b_counts_file, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)

# FASTA: >sequence as header, sequence as line
fasta_lines_b <- unlist(lapply(merged_seqs, function(seq) c(paste0(">", seq), seq)))
writeLines(fasta_lines_b, b_fasta_file)

cat("Mergeable _b dataset written:\n", b_counts_file, "\n", b_fasta_file, "\n")

# ----------------------------
# Write numbered ASV dataset (_num)
# ----------------------------
num_counts_file <- paste0(prefix, "_merged_ASVs_num.tsv")
num_fasta_file  <- paste0(prefix, "_merged_ASVs_num.fa")

# Create numbered ASVs
num_ids <- paste0("ASV_", seq_len(nrow(merged_counts)))
rownames(merged_counts) <- num_ids

# Counts table with ASV_# rownames
write.table(merged_counts, num_counts_file, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)

# FASTA: >ASV_# as header, sequence as line
fasta_lines_num <- unlist(lapply(seq_along(merged_seqs), function(i) {
  c(paste0(">", num_ids[i]), merged_seqs[i])
}))
writeLines(fasta_lines_num, num_fasta_file)

cat("Numbered ASV dataset written:\n", num_counts_file, "\n", num_fasta_file, "\n")

# ----------------------------
# Optional: create phyloseq RDS
# ----------------------------
if(requireNamespace("phyloseq", quietly = TRUE)){
  library(phyloseq)

  counts_matrix <- as.matrix(merged_counts)
  mode(counts_matrix) <- "numeric"

  # Simple sample metadata
  sample_metadata <- data.frame(SampleID = colnames(counts_matrix),
                                stringsAsFactors = FALSE)
  rownames(sample_metadata) <- sample_metadata$SampleID

  ps <- phyloseq(otu_table(counts_matrix, taxa_are_rows = TRUE),
                 sample_data(sample_metadata))

  phyloseq_file <- paste0(prefix, "_merged_phyloseq.rds")
  saveRDS(ps, phyloseq_file)
  cat("Phyloseq RDS written:", phyloseq_file, "\n")
} else {
  cat("phyloseq not installed; skipping RDS creation.\n")
}
