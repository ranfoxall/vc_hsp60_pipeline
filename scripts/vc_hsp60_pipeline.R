# ==========================================
# 16s_v4v5_pipeline: End-to-end 16S V4-V5 microbiome pipeline
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
option_list <- list(
  make_option(c("--reads_path"),      type="character",
              help="Path to trimmed paired FASTQ reads"),
  make_option(c("--output_prefix"),   type="character",
              help="Prefix for output files"),
  make_option(c("--error_model"),     type="character", default="loess",
              help="Error model: 'loess' (default), 'default', or 'compare'"),
  make_option(c("--error_model_rds"), type="character", default=NULL,
              help="(Optional) Path to a previously saved error model .rds — skips error learning"),
  make_option(c("--pool"),            type="character", default="pseudo",
              help="DADA2 pooling method: 'pseudo' (default), 'FALSE' (fastest), 'TRUE' (most sensitive)"),
  make_option(c("--trim_left_f"),     type="integer",   default=0,
              help="Bases to trim from 5' end of forward reads (default: 0)"),
  make_option(c("--trim_left_r"),     type="integer",   default=0,
              help="Bases to trim from 5' end of reverse reads (default: 0)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$reads_path) | is.null(opt$output_prefix)) {
  stop("Please provide --reads_path and --output_prefix")
}

if (!opt$error_model %in% c("loess", "default", "compare")) {
  stop("--error_model must be one of: loess, default, compare")
}

reads_path      <- opt$reads_path
output_prefix   <- opt$output_prefix
error_model     <- opt$error_model
error_model_rds <- opt$error_model_rds
pool_method     <- opt$pool
trim_left_f     <- opt$trim_left_f
trim_left_r     <- opt$trim_left_r
if (pool_method == "TRUE")  pool_method <- TRUE
if (pool_method == "FALSE") pool_method <- FALSE

cat("Reading FASTQ files from:", reads_path, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Error model:", error_model, "\n")
if (!is.null(error_model_rds)) cat("Error model RDS:", error_model_rds, "\n")
cat("Pooling method:", pool_method, "\n")
cat("Trim left F:", trim_left_f, "\n")
cat("Trim left R:", trim_left_r, "\n")

# [2] read files -------------------------------------------------------
fnFs <- sort(list.files(reads_path, pattern="_R1_001\\.fastq\\.gz$", full.names=TRUE))
fnRs <- sub("_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz", fnFs)

paired <- file.exists(fnRs)
fnFs <- fnFs[paired]
fnRs <- fnRs[paired]

stopifnot(length(fnFs) == length(fnRs))

sample_names <- sub("_S[0-9]+_L[0-9]+_R1_001\\.fastq\\.gz$", "", basename(fnFs))
names(fnFs) <- sample_names
names(fnRs) <- sample_names

cat("Found", length(fnFs), "samples:\n")
print(sample_names)

# exclude very small files before any processing
min_size <- 10000
big_enough <- file.size(fnFs) >= min_size
if (any(!big_enough)) {
  small <- sample_names[!big_enough]
  cat("WARNING: The following samples have very small input files and will be excluded before processing:\n")
  cat(paste("  ", small, collapse="\n"), "\n")
  fnFs         <- fnFs[big_enough]
  fnRs         <- fnRs[big_enough]
  sample_names <- sample_names[big_enough]
  cat(length(fnFs), "samples remain after size filter.\n")
}

# [3] quality profiles and truncation lengths --------------------------
# truncation lengths are calculated empirically at median Q25 position
# across representative samples, then clamped to floors that guarantee
# enough overlap for the V4-V5 amplicon (~411bp after primer removal).
#
# Floors: forward >= 220bp, reverse >= 180bp  (total >= 400bp)
# Hard cap: 251bp (standard Illumina 2x250 run)
#
# If median quality never drops below Q25 the hard cap is used.
# If the amplicon is shorter on your run, check the quality profiles PDF
# and adjust via --trim_left_f / --trim_left_r if primers remain on reads.

repFs <- fnFs[1:min(3, length(fnFs))]
repRs <- fnRs[1:min(3, length(fnRs))]

pdf(file = paste0(output_prefix, "_quality_profiles.pdf"))
plotQualityProfile(repFs)
plotQualityProfile(repRs)
dev.off()
cat("Quality profiles written to:", paste0(output_prefix, "_quality_profiles.pdf"), "\n")

calc_trunc_len <- function(fastq_files, min_len=200, max_len=251, q_threshold=25) {
  qual_list <- lapply(fastq_files, function(f) {
    fq <- ShortRead::readFastq(f)
    as(quality(fq), "matrix")
  })
  qual_mat <- do.call(rbind, qual_list)
  med_qual <- apply(qual_mat, 2, median)
  passing  <- which(med_qual >= q_threshold)
  if (length(passing) == 0) {
    cat("WARNING: No positions met Q >=", q_threshold, "-- using min_len floor.\n")
    return(min_len)
  }
  trunc <- max(passing)
  trunc <- max(min_len, min(trunc, max_len))
  return(trunc)
}

truncLenF <- calc_trunc_len(repFs, min_len=220, max_len=251)
truncLenR <- calc_trunc_len(repRs, min_len=180, max_len=251)

cat("Auto-detected truncation lengths:\n")
cat("  Forward:", truncLenF, "\n")
cat("  Reverse:", truncLenR, "\n")
cat("  Combined:", truncLenF + truncLenR, "bp  (minimum overlap target: 430bp)\n")

if (truncLenF + truncLenR < 430) {
  cat("WARNING: Combined truncation length", truncLenF + truncLenR,
      "bp may be insufficient for V4-V5 overlap (~20bp needed).\n")
  cat("  Check quality profiles and consider relaxing --trim_left values or re-sequencing.\n")
}

# [4] filter and trim --------------------------------------------------
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
  trimLeft    = c(trim_left_f, trim_left_r),
  truncLen    = c(truncLenF, truncLenR),
  maxN        = 0,
  maxEE       = c(2, 2),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = TRUE
)
cat("Read counts after filtering:\n")
print(out)

kept    <- out[, "reads.out"] > 0 & file.exists(filtFs) & file.size(filtFs) > 0
dropped <- sample_names[!kept]
if (length(dropped) > 0) {
  cat("WARNING: The following samples had zero reads after filtering and will be excluded:\n")
  cat(paste("  ", dropped, collapse="\n"), "\n")
  write(dropped, paste0(output_prefix, "_dropped_samples.txt"))
  cat("Dropped sample names written to:", paste0(output_prefix, "_dropped_samples.txt"), "\n")
}
filtFs       <- filtFs[kept]
filtRs       <- filtRs[kept]
sample_names <- sample_names[kept]

if (length(filtFs) == 0) stop("No reads survived filtering. Check truncation lengths and input quality.")
cat(length(filtFs), "samples remain after filtering.\n")

# [5] error model ------------------------------------------------------
# cached to .rds after first run -- delete to force relearning
# loess model (default) fits a smoothed curve with span=0.95 and generally
# outperforms the DADA2 default on 16S data

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

error_fit_metric <- function(errObj) {
  trans      <- errObj$trans
  trans_norm <- sweep(trans, 2, colSums(trans), "/")
  err_mat    <- errObj$err_out
  err_mat    <- err_mat[rownames(trans_norm), colnames(trans_norm)]
  return(sum((err_mat - trans_norm)^2))
}

err_cache_file <- paste0(output_prefix, "_error_model.rds")

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
    cat("Learning error rates using LOESS model...\n")
    errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e9,
                        errorEstimationFunction=loessErrfun_mod)
    errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e9,
                        errorEstimationFunction=loessErrfun_mod)
    cat("LOESS error model complete.\n\n")

  } else if (error_model == "default") {
    cat("Learning error rates using DADA2 default model...\n")
    errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e9)
    errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e9)
    cat("Default error model complete.\n\n")

  } else if (error_model == "compare") {
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

  cat("Saving error model to:", err_cache_file, "\n")
  saveRDS(list(errF=errF, errR=errR), err_cache_file)
}

# [6] dereplicate and denoise ------------------------------------------
cat("Dereplicating reads...\n")
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample_names
names(derepRs) <- sample_names
cat("Dereplication complete.\n")

# denoising cached to .rds -- reused automatically on reruns
dada_cache_file <- paste0(output_prefix, "_dada_objects.rds")

if (file.exists(dada_cache_file)) {
  cat("Loading saved DADA2 objects from:", dada_cache_file, "\n")
  dada_cache <- readRDS(dada_cache_file)
  dadaFs <- dada_cache$dadaFs
  dadaRs <- dada_cache$dadaRs
  cat("DADA2 objects loaded.\n")
} else {
  cat("DADA2 denoising (this is the slow step)...\n")
  cat("Denoising forward reads...\n")
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=pool_method)
  cat("Denoising reverse reads...\n")
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=pool_method)
  cat("Saving DADA2 objects to:", dada_cache_file, "\n")
  saveRDS(list(dadaFs=dadaFs, dadaRs=dadaRs), dada_cache_file)
  cat("Denoising complete.\n")
}

# [7] merge paired reads -----------------------------------------------
cat("Merging paired reads...\n")
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names
merge_failed <- c()

for (s in sample_names) {
  mergers[[s]] <- tryCatch({
    mergePairs(dadaFs[[s]], derepFs[[s]], dadaRs[[s]], derepRs[[s]],
               minOverlap=20, maxMismatch=0, verbose=FALSE)
  }, error = function(e) {
    cat("WARNING: Merging failed for sample", s, "--", conditionMessage(e), "\n")
    merge_failed <<- c(merge_failed, s)
    NULL
  })
}

if (length(merge_failed) > 0) {
  cat("WARNING:", length(merge_failed), "sample(s) failed merging and will be excluded:\n")
  cat(paste("  ", merge_failed, collapse="\n"), "\n")
  write(merge_failed, paste0(output_prefix, "_merge_failed_samples.txt"))
  mergers      <- mergers[!names(mergers) %in% merge_failed]
  dadaFs       <- dadaFs[!names(dadaFs) %in% merge_failed]
  sample_names <- sample_names[!sample_names %in% merge_failed]
}

if (length(mergers) == 0) stop("No samples survived merging.")
cat("Merging complete.", length(mergers), "samples retained.\n")

cat("Making sequence table...\n")
seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")

# [8] length filter ----------------------------------------------------
# expected V4-V5 amplicon range after primer removal: 370-450 bp
# (515F/926R target is ~411bp; range accommodates natural length variation)
cat("Filtering by amplicon length (370-450 bp)...\n")
target_range <- 370:450
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% target_range]
cat("Sequences in target length range (370-450 bp):", ncol(seqtab), "\n")

# [9] chimera removal --------------------------------------------------
# chimera removal is enabled by default for 16S V4-V5 data.
# uses consensus method -- calls chimeras only when the majority of samples
# agree, reducing false positive removal.
# if >50% of reads are removed, primers may still be present on reads --
# check that cutadapt trimming ran successfully.

cat("Removing chimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
n_before <- sum(seqtab)
n_after  <- sum(seqtab_nochim)
pct_chim <- round(100 * (1 - n_after / n_before), 1)
cat("Reads before chimera removal:", n_before, "\n")
cat("Reads after  chimera removal:", n_after, "(", pct_chim, "% removed)\n")
cat("ASVs before:", ncol(seqtab), "  after:", ncol(seqtab_nochim), "\n")

if (pct_chim > 50) {
  cat("WARNING: >50% of reads removed as chimeras.\n")
  cat("  This may indicate primer sequences remain on reads.\n")
  cat("  Check that cutadapt trimming used the correct primer sequences:\n")
  cat("    515F: GTGYCAGCMGCCGCGGTAA\n")
  cat("    926R: CCGYCAATTYMTTTRAGTTT\n")
}

# [9b] read tracking ---------------------------------------------------
cat("Generating read tracking table...\n")
get_n <- function(x) sum(getUniques(x))
track <- data.frame(
  sample    = sample_names,
  reads_in  = out[kept, "reads.in"],
  filtered  = out[kept, "reads.out"],
  denoised  = sapply(dadaFs, get_n)[sample_names],
  merged    = sapply(mergers, get_n)[sample_names],
  nonchim   = rowSums(seqtab_nochim)[sample_names]
)
track[is.na(track)] <- 0
track$pct_retained <- round(track$nonchim / track$reads_in * 100, 1)
track_file <- paste0(output_prefix, "_read_tracking.tsv")
write.table(track, track_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Read tracking table:\n")
print(track)
cat("Read tracking written to:", track_file, "\n")

# [10] ASV tables ------------------------------------------------------
cat("Writing ASV tables...\n")
asv_seqs <- colnames(seqtab_nochim)
rownames(seqtab_nochim) <- sample_names

asv_tab_seq <- t(seqtab_nochim)
rownames(asv_tab_seq) <- asv_seqs
write.table(asv_tab_seq,
            paste0(output_prefix, "_Counts_seqASV_b.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

asv_ids     <- paste0("ASV_", seq_along(asv_seqs))
asv_tab_num <- asv_tab_seq
rownames(asv_tab_num) <- asv_ids
write.table(asv_tab_num,
            paste0(output_prefix, "_Counts_numASV.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

# [11] FASTA output ----------------------------------------------------
cat("Writing FASTA...\n")
asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
asv_fasta      <- c(rbind(paste0(">", asv_ids), asv_seqs))
write(asv_fasta, asv_fasta_path)
cat("FASTA written to:", asv_fasta_path, "\n")

# [12] phyloseq object -------------------------------------------------
cat("Creating phyloseq object...\n")
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

# [13] summary ---------------------------------------------------------
cat("\nPipeline finished.\n")
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
cat("\nNext step: run taxonomy classification using run_16s_classify.slurm\n")
cat("or: Rscript scripts/16s_classify.R --output_prefix", output_prefix,
    "--db both --silva_classifier /path/to/silva-classifier.qza",
    "--gg2_classifier /path/to/gg2-classifier.qza\n")
