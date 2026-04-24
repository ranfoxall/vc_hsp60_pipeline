# vc_hsp60_add_taxonomy.R
# Merges taxonomy and refseq into phyloseq RDS files produced by vc_hsp60_pipeline.R.
# Run from the vc_hsp60_pipeline conda environment after vc_hsp60_classify.R completes.
#
# For the vibrio method, produces two phyloseq objects:
#   {prefix}_phyloseq_vibrio.rds — Vibrio ASVs only (BLAST screen + Vibrio classifier)
#   {prefix}_phyloseq_full.rds   — all ASVs with full taxonomy
#
# For the sklearn method, produces one phyloseq object:
#   {prefix}_phyloseq_taxonomy.rds — all ASVs with taxonomy
#
# Called automatically by run_vchsp60_full.slurm and run_vchsp60_classify.slurm.
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Funder:     Dr. Cheryl Whistler, Dr. Stephen Jones
# Lab:        Whistler Lab, Dept. of Molecular, Cellular and Biomedical Sciences, UNH
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0

suppressPackageStartupMessages(library(phyloseq))

args <- commandArgs(trailingOnly=TRUE)

parse_args_simple <- function(args) {
  result <- list(output_prefix=NULL, phyloseq_rds=NULL, method="vibrio")
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--output_prefix" && i < length(args)) {
      result$output_prefix <- args[i+1]; i <- i + 2
    } else if (args[i] == "--phyloseq_rds" && i < length(args)) {
      result$phyloseq_rds <- args[i+1]; i <- i + 2
    } else if (args[i] == "--method" && i < length(args)) {
      result$method <- args[i+1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  result
}

opt    <- parse_args_simple(args)
if (is.null(opt$output_prefix)) stop("Please provide --output_prefix")

prefix   <- opt$output_prefix
method   <- opt$method
ps_file  <- if (!is.null(opt$phyloseq_rds)) opt$phyloseq_rds else paste0(prefix, "_phyloseq.rds")
tax_file <- paste0(prefix, "_taxonomy_table.csv")
fa_file  <- paste0(prefix, "_ASVs.fa")

for (f in c(ps_file, tax_file, fa_file)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

cat("Loading phyloseq object:", ps_file, "\n")
ps <- readRDS(ps_file)

cat("Loading taxonomy table:", tax_file, "\n")
tax_mat <- as.matrix(read.csv(tax_file, row.names=1))

cat("Loading ASV sequences:", fa_file, "\n")
fa <- Biostrings::readDNAStringSet(fa_file)

# full phyloseq object — all ASVs with taxonomy
cat("Building full phyloseq object...\n")
ps_full <- merge_phyloseq(ps, tax_table(tax_mat), refseq(fa))

if (method == "vibrio") {
  # full object
  full_file <- paste0(prefix, "_phyloseq_full.rds")
  saveRDS(ps_full, full_file)
  cat("Full phyloseq object written to:", full_file, "\n")

  # Vibrio-only object — subset to ASVs with Genus containing "vibrio" (case insensitive)
  # this captures both Vibrio and Aliivibrio if present
  cat("Building Vibrio-only phyloseq object...\n")
  tax_df       <- as.data.frame(tax_table(ps_full))
  vibrio_asvs  <- rownames(tax_df)[grepl("vibrio", tax_df$Genus, ignore.case=TRUE)]

  if (length(vibrio_asvs) > 0) {
    ps_vibrio    <- prune_taxa(vibrio_asvs, ps_full)
    # remove samples with zero Vibrio reads
    ps_vibrio    <- prune_samples(sample_sums(ps_vibrio) > 0, ps_vibrio)
    vibrio_file  <- paste0(prefix, "_phyloseq_vibrio.rds")
    saveRDS(ps_vibrio, vibrio_file)
    cat("Vibrio-only phyloseq object written to:", vibrio_file, "\n")
    cat("  Vibrio ASVs:   ", ntaxa(ps_vibrio), "\n")
    cat("  Samples:       ", nsamples(ps_vibrio), "\n")
  } else {
    cat("WARNING: No Vibrio ASVs found in taxonomy table — Vibrio-only RDS not produced.\n")
  }

} else {
  # sklearn method — single output
  out_file <- paste0(prefix, "_phyloseq_taxonomy.rds")
  saveRDS(ps_full, out_file)
  cat("Phyloseq with taxonomy written to:", out_file, "\n")
}

cat("\nDone.\n")
