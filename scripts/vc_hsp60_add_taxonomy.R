# vc_hsp60_add_taxonomy.R
# Merges taxonomy table and refseq into the phyloseq RDS produced by vc_hsp60_pipeline.R.
# Run from the vc_hsp60_pipeline conda environment (where phyloseq is installed),
# after vc_hsp60_classify.R has produced the taxonomy_table.csv.
#
# Called automatically by run_vchsp60_full.slurm and run_vchsp60_classify.slurm.
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Funder:     Dr. Cheryl Whistler, Dr. Stephen Jones
# Lab:        Whistler Lab, Dept. of Molecular, Cellular and Biomedical Sciences, UNH
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0

library(phyloseq)
args <- commandArgs(trailingOnly=TRUE)

parse_args_simple <- function(args) {
  result <- list(output_prefix=NULL, phyloseq_rds=NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--output_prefix" && i < length(args)) {
      result$output_prefix <- args[i+1]; i <- i + 2
    } else if (args[i] == "--phyloseq_rds" && i < length(args)) {
      result$phyloseq_rds <- args[i+1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  result
}

opt <- parse_args_simple(args)
if (is.null(opt$output_prefix)) stop("Please provide --output_prefix")

prefix   <- opt$output_prefix
ps_file  <- if (!is.null(opt$phyloseq_rds)) opt$phyloseq_rds else paste0(prefix, "_phyloseq.rds")
tax_file <- paste0(prefix, "_taxonomy_table.csv")
fa_file  <- paste0(prefix, "_ASVs.fa")
out_file <- paste0(prefix, "_phyloseq_taxonomy.rds")

if (!file.exists(ps_file))  stop("Phyloseq RDS not found: ", ps_file)
if (!file.exists(tax_file)) stop("Taxonomy table not found: ", tax_file)
if (!file.exists(fa_file))  stop("FASTA file not found: ", fa_file)

cat("Loading phyloseq object:", ps_file, "\n")
ps <- readRDS(ps_file)

cat("Loading taxonomy table:", tax_file, "\n")
tax_mat <- as.matrix(read.csv(tax_file, row.names=1))

cat("Loading ASV sequences:", fa_file, "\n")
fa <- Biostrings::readDNAStringSet(fa_file)

cat("Merging taxonomy and refseq into phyloseq object...\n")
ps <- merge_phyloseq(ps, tax_table(tax_mat), refseq(fa))
saveRDS(ps, out_file)
cat("Phyloseq with taxonomy and refseq written to:", out_file, "\n")
