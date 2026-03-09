# ==========================================
# vc_hsp60_pipeline: End-to-end Hsp60 microbiome pipeline
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Funder:     Dr. Cheryl Whistler, Dr. Stephen Jones
# Lab:        Whistler Lab, Dept. of Molecular, Cellular and Biomedical Sciences, UNH
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0
# ==========================================

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(dplyr)
  library(magrittr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
  library(phyloseq)
})

# [1] arguments -------------------------------------------------------
# parse and validate command line inputs
option_list <- list(
  make_option(c("--reads_path"), type="character", help="Path to raw FASTQ reads"),
  make_option(c("--output_prefix"), type="character", help="Prefix for output files"),
  make_option(c("--error_model"), type="character", default="loess",
              help="Error model: 'loess' (default), 'default', or 'compare'"),
  make_option(c("--error_model_rds"), type="character", default=NULL,
              help="(Optional) Path to a previously saved error model .rds — skips error learning")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$reads_path) | is.null(opt$output_prefix)) {
  stop("Please provide --reads_path and --output_prefix")
}

if (!opt$error_model %in% c("loess", "default", "compare")) {
  stop("--error_model must be one of: loess, default, compare")
}

reads_path       <- opt$reads_path
output_prefix    <- opt$output_prefix
error_model      <- opt$error_model
error_model_rds  <- opt$error_model_rds

cat("Reading FASTQ files from:", reads_path, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Error model:", error_model, "\n")
if (!is.null(error_model_rds)) cat("Error model RDS:", error_model_rds, "\n")

# [2] read files -------------------------------------------------------
# find paired FASTQ files and extract sample names
fnFs <- sort(list.files(reads_path, pattern="_R1_001\\.fastq\\.gz$", full.names=TRUE))
fnRs <- sub("_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz", fnFs)

paired <- file.exists(fnRs)
fnFs <- fnFs[paired]
fnRs <- fnRs[paired]

stopifnot(length(fnFs) == length(fnRs))

# Robust sample name extraction: strip _S##_L###_R1_001.fastq.gz
sample_names <- sub("_S[0-9]+_L[0-9]+_R1_001\\.fastq\\.gz$", "", basename(fnFs))
names(fnFs) <- sample_names
names(fnRs) <- sample_names

cat("Found", length(fnFs), "samples:\n")
print(sample_names)

# [3] quality profiles and truncation lengths --------------------------
# plot quality profiles for representative samples and auto-calculate
# truncation lengths based on median Q30 position

# Choose a few representative FASTQ files (forward and reverse)
repFs <- fnFs[1:min(3, length(fnFs))]
repRs <- fnRs[1:min(3, length(fnRs))]

# Generate quality profiles
pdf(file = paste0(output_prefix, "_quality_profiles.pdf"))
plotQualityProfile(repFs)
plotQualityProfile(repRs)
dev.off()

# Function to calculate truncation at median Q30 threshold
calc_trunc_len <- function(fastq_files, min_len=100, max_len=251, q_threshold=30) {
  qual_list <- lapply(fastq_files, function(f) {
    fastq <- ShortRead::readFastq(f)
    as(quality(fastq), "matrix")
  })
  qual_mat <- do.call(rbind, qual_list)
  med_qual <- apply(qual_mat, 2, median)
  trunc <- max(which(med_qual >= q_threshold))
  trunc <- max(min_len, min(trunc, max_len))
  return(trunc)
}

truncLenF <- calc_trunc_len(repFs, min_len=200, max_len=251)
truncLenR <- calc_trunc_len(repRs, min_len=200, max_len=251)

cat("Run-level truncation lengths:\n")
cat("Forward:", truncLenF, "\n")
cat("Reverse:", truncLenR, "\n")

# [4] filter and trim --------------------------------------------------
# remove low-quality reads, trim to truncation lengths, drop empty samples
filt_path <- paste0(output_prefix, "_filtered")
if (!dir.exists(filt_path)) {
  dir.create(filt_path)
  cat("Created filtered directory at:", filt_path, "\n")
} else {
  cat("Filtered directory already exists at:", filt_path, "\n")
}

filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

cat("Filtering and trimming reads...\n")

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen  = c(truncLenF, truncLenR),
  maxN      = 0,
  maxEE     = c(2, 2),
  truncQ    = 2,
  rm.phix   = TRUE,
  compress  = TRUE,
  multithread = TRUE
)
cat("Read counts after filtering:\n")
print(out)

# Drop samples with zero reads after filtering — these will crash DADA2
kept    <- out[, "reads.out"] > 0 & file.exists(filtFs) & file.size(filtFs) > 0
dropped <- sample_names[!kept]
if (length(dropped) > 0) {
  cat("WARNING: The following samples had zero reads after filtering and will be excluded:\n")
  cat(paste(" ", dropped, collapse="\n"), "\n")
  write(dropped, paste0(output_prefix, "_dropped_samples.txt"))
  cat("Dropped sample names written to:", paste0(output_prefix, "_dropped_samples.txt"), "\n")
}
filtFs       <- filtFs[kept]
filtRs       <- filtRs[kept]
sample_names <- sample_names[kept]

if (length(filtFs) == 0) stop("No reads survived filtering. Check truncation lengths and input quality.")
cat(length(filtFs), "samples remain after filtering.\n")

# [5] error model ------------------------------------------------------
# learn DADA2 error rates from the data; cached to .rds after first run
# to save time on reruns — delete the .rds file to force relearning

# LOESS error function (defined here for use by any model option)
loessErrfun_mod <- function(trans) {
  qq  <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=16, ncol=length(qq))
  rownames(est) <- paste0(rep(c("A","C","G","T"), each=4), "2",
                          rep(c("A","C","G","T"), 4))
  colnames(est) <- colnames(trans)
  for (nti in c("A","C","G","T")) {
    for (ntj in c("A","C","G","T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot  <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        df   <- data.frame(q=qq, rlogp=log10((errs+1)/(tot+1)))
        mod.lo <- loess(rlogp ~ q, df, span=0.95, degree=1)
        pred   <- predict(mod.lo, qq)
        pred[is.na(pred)] <- 0
        est[paste0(nti,"2",ntj),] <- 10^pred
      }
    }
  }
  est[est > 0.25] <- 0.25
  est[est < 1e-7] <- 1e-7
  return(est)
}

# Model fit metric (used only for 'compare' mode)
error_fit_metric <- function(errObj) {
  trans      <- errObj$trans
  trans_norm <- sweep(trans, 2, colSums(trans), "/")
  err_mat    <- errObj$err_out
  err_mat    <- err_mat[rownames(trans_norm), colnames(trans_norm)]
  return(sum((err_mat - trans_norm)^2))
}

err_cache_file <- paste0(output_prefix, "_error_model.rds")

# load error model from --error_model_rds if provided, then auto-cache, then learn fresh
if (!is.null(error_model_rds)) {
  if (!file.exists(error_model_rds)) stop("Error model RDS not found: ", error_model_rds)
  cat("Loading error model from --error_model_rds:", error_model_rds, "\n")
  err_cache <- readRDS(error_model_rds)
  errF <- err_cache$errF
  errR <- err_cache$errR
  cat("Error model loaded.\n\n")

} else if (file.exists(err_cache_file)) {
  cat("Loading saved error model from:", err_cache_file, "\n")
  err_cache <- readRDS(err_cache_file)
  errF <- err_cache$errF
  errR <- err_cache$errR
  cat("Error model loaded.\n\n")

} else {

  if (error_model == "loess") {
    # --- LOESS only (default, recommended for Hsp60 amplicon data) ---
    cat("Learning error rates using LOESS model...\n")
    errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e9,
                        errorEstimationFunction=loessErrfun_mod)
    errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e9,
                        errorEstimationFunction=loessErrfun_mod)
    cat("LOESS error model complete.\n\n")

  } else if (error_model == "default") {
    # --- DADA2 default model only ---
    cat("Learning error rates using DADA2 default model...\n")
    errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e9)
    errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e9)
    cat("Default error model complete.\n\n")

  } else if (error_model == "compare") {
    # --- Run both models and select best fit ---
    cat("Learning error rates using both default and LOESS models (this will take longer)...\n")

    cat("  Running default model...\n")
    errF_default <- learnErrors(filtFs, multithread=TRUE, nbases=1e9)
    errR_default <- learnErrors(filtRs, multithread=TRUE, nbases=1e9)

    cat("  Running LOESS model...\n")
    errF_loess <- learnErrors(filtFs, multithread=TRUE, nbases=1e9,
                              errorEstimationFunction=loessErrfun_mod)
    errR_loess <- learnErrors(filtRs, multithread=TRUE, nbases=1e9,
                              errorEstimationFunction=loessErrfun_mod)

    fit_default_F <- error_fit_metric(errF_default)
    fit_loess_F   <- error_fit_metric(errF_loess)
    fit_default_R <- error_fit_metric(errR_default)
    fit_loess_R   <- error_fit_metric(errR_loess)

    cat("Forward fit scores:  Default =", fit_default_F, " | LOESS =", fit_loess_F, "\n")
    cat("Reverse fit scores:  Default =", fit_default_R, " | LOESS =", fit_loess_R, "\n")

    errF <- if (fit_loess_F < fit_default_F) {
      cat("Selected LOESS for Forward.\n"); errF_loess
    } else {
      cat("Selected DEFAULT for Forward.\n"); errF_default
    }
    errR <- if (fit_loess_R < fit_default_R) {
      cat("Selected LOESS for Reverse.\n"); errR_loess
    } else {
      cat("Selected DEFAULT for Reverse.\n"); errR_default
    }
    cat("Error model selection complete.\n\n")
  }

  # save error model for reuse on reruns
  cat("Saving error model to:", err_cache_file, "\n")
  saveRDS(list(errF=errF, errR=errR), err_cache_file)
}

# [6] dereplicate and denoise ------------------------------------------
# collapse identical reads, then run DADA2 to infer true sequences
cat("\n--- Step 6: Dereplicating reads ---\n")
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample_names
names(derepRs) <- sample_names
cat("Dereplication complete.\n")

cat("\n--- Step 6b: DADA2 denoising (this is the slow step) ---\n")
cat("Denoising forward reads...\n")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
cat("Denoising reverse reads...\n")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
cat("Denoising complete.\n")

cat("\n--- Step 6c: Merging paired reads ---\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
cat("Merging complete.\n")

cat("\n--- Step 6d: Making sequence table ---\n")
seqtab  <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")

# [7] length filter ----------------------------------------------------
# keep only sequences within the expected Hsp60 amplicon range (390-475 bp)
cat("\n--- Step 7: Filtering by amplicon length (390-475 bp) ---\n")
target_range <- 390:475
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% target_range]
cat("Sequences in target length range (390-475 bp):", ncol(seqtab), "\n")

# [8] chimera removal --------------------------------------------------
# remove chimeric sequences using consensus method
cat("\n--- Step 8: Removing chimeras ---\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
cat("Sequences after chimera removal:", ncol(seqtab_nochim), "\n")

# [8b] read tracking ---------------------------------------------------
# summarize reads retained at each step per sample; saved as TSV
cat("\n--- Step 8b: Generating read tracking table ---\n")
get_n <- function(x) sum(getUniques(x))
track <- data.frame(
  sample    = sample_names,
  reads_in  = out[kept, "reads.in"],
  filtered  = out[kept, "reads.out"],
  denoised  = sapply(dadaFs, get_n),
  merged    = sapply(mergers, get_n),
  nonchim   = rowSums(seqtab_nochim)
)
track$pct_retained <- round(track$nonchim / track$reads_in * 100, 1)
track_file <- paste0(output_prefix, "_read_tracking.tsv")
write.table(track, track_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Read tracking table:\n")
print(track)
cat("Read tracking written to:", track_file, "\n")

# [9] ASV tables -------------------------------------------------------
# write counts tables with sequence-based and numbered ASV row names
cat("\n--- Step 9: Writing ASV tables ---\n")
asv_seqs <- colnames(seqtab_nochim)
rownames(seqtab_nochim) <- sample_names

# Sequence-based ASVs
asv_tab_seq <- t(seqtab_nochim)
rownames(asv_tab_seq) <- asv_seqs
write.table(asv_tab_seq,
            paste0(output_prefix, "_Counts_seqASV_b.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

# Numbered ASVs
asv_ids     <- paste0("ASV_", seq_along(asv_seqs))
asv_tab_num <- asv_tab_seq
rownames(asv_tab_num) <- asv_ids
write.table(asv_tab_num,
            paste0(output_prefix, "_Counts_numASV.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

# [10] FASTA output ----------------------------------------------------
# write numbered ASV sequences to FASTA for use in classification step
cat("\n--- Step 10: Writing FASTA ---\n")
asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
asv_fasta      <- c(rbind(paste0(">", asv_ids), asv_seqs))
write(asv_fasta, asv_fasta_path)
cat("FASTA written to:", asv_fasta_path, "\n")

# [11] phyloseq object -------------------------------------------------
# create phyloseq object with ASV counts and sample metadata
# taxonomy is added in the classification step (vc_hsp60_classify.R)
cat("\n--- Step 11: Creating phyloseq object ---\n")
sample_metadata <- data.frame(
  SampleID   = colnames(asv_tab_num),
  TotalReads = colSums(asv_tab_num),
  row.names  = colnames(asv_tab_num),
  stringsAsFactors = FALSE
)

ps <- phyloseq(
  otu_table(as.matrix(asv_tab_num), taxa_are_rows=TRUE),
  sample_data(sample_metadata)
)

ps_file <- paste0(output_prefix, "_phyloseq.rds")
saveRDS(ps, ps_file)
cat("Phyloseq object written to:", ps_file, "\n")

# [12] summary ---------------------------------------------------------
# print output file paths and next-step instructions
cat("\n=== Pipeline finished ===\n")
cat("\nOutputs:\n")
cat("1. Sequence-based ASV table:", paste0(output_prefix, "_Counts_seqASV_b.tsv"), "\n")
cat("2. Numbered ASV table:      ", paste0(output_prefix, "_Counts_numASV.tsv"), "\n")
cat("3. FASTA for classifier:    ", asv_fasta_path, "\n")
cat("4. Phyloseq object:         ", ps_file, "\n")
cat("5. Quality profiles PDF:    ", paste0(output_prefix, "_quality_profiles.pdf"), "\n")
cat("6. Read tracking table:     ", track_file, "\n")
if (length(dropped) > 0) {
  cat("7. Dropped samples log:     ", paste0(output_prefix, "_dropped_samples.txt"), "\n")
}
cat("\nNext step: run taxonomy classification using run_vchsp60_classify.slurm\n")
cat("or: Rscript vc_hsp60_classify.R --output_prefix", output_prefix, "--classifier /path/to/cpn60_classifier_v11.qza\n")
