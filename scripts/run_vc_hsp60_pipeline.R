# ==========================================
# vc_hsp60_pipeline: End-to-end Hsp60 microbiome pipeline
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

# ----------------------
# 1. Parse command line arguments
# ----------------------
option_list <- list(
  make_option(c("--reads_path"), type="character", help="Path to raw FASTQ reads"),
  make_option(c("--output_prefix"), type="character", help="Prefix for output files")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$reads_path) | is.null(opt$output_prefix)) {
  stop("Please provide --reads_path and --output_prefix")
}

reads_path <- opt$reads_path
output_prefix <- opt$output_prefix

cat("Reading FASTQ files from:", reads_path, "\n")
cat("Output prefix:", output_prefix, "\n")

# ----------------------
# 2. List forward and reverse reads
# ----------------------
fnFs <- sort(list.files(reads_path, pattern="_R1_001.fastq.gz", full.names=TRUE))
fnRs <- sub("_R1_001.fastq.gz", "_R2_001.fastq.gz", fnFs)

paired <- file.exists(fnRs)
fnFs <- fnFs[paired]
fnRs <- fnRs[paired]

stopifnot(length(fnFs) == length(fnRs))

sample_names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
names(fnFs) <- sample_names
names(fnRs) <- sample_names

# ----------------------
# 3. Adaptive truncation with minimum length safeguard
# ----------------------

get_trunc <- function(fn_list, min_len = 200, qual_threshold = 30) {

  cat("Checking first FASTQ file:\n")
  cat("  ", basename(fn_list[1]), "\n")

  fq <- readFastq(fn_list[1])
  n_reads <- length(fq)
  read_width <- width(fq)[1]
  cat("  Number of reads:", n_reads, "; read width:", read_width, "\n")

  qs <- as(quality(fq), "matrix")
  med_q <- apply(qs, 2, median)

  pos <- which(med_q < qual_threshold)
  trunc_pos <- if(length(pos) == 0) ncol(qs) else pos[1] - 1

  # SAFEGUARD: do not truncate below minimum length
  if(trunc_pos < min_len) {
    cat("  Quality drops early, but enforcing minimum length:", min_len, "\n")
    trunc_pos <- min_len
  }

  cat("  Final truncation position:", trunc_pos, "\n")
  return(trunc_pos)
}

# Your enforced minimums
truncLenF <- get_trunc(fnFs, min_len = 240)
truncLenR <- get_trunc(fnRs, min_len = 220)

cat("Truncation lengths: forward =", truncLenF,
    ", reverse =", truncLenR, "\n")

# ----------------------
# 4. Filter reads
# ----------------------
filt_path <- file.path(reads_path, paste0(output_prefix, "_filtered"))
if(!dir.exists(filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(truncLenF, truncLenR),
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# ----------------------
# 5. Adaptive error model selection
# ----------------------
cat("Learning error rates using multiple models...\n")

# Default
errF_default <- learnErrors(filtFs, multithread=TRUE, nbases=1e8)
errR_default <- learnErrors(filtRs, multithread=TRUE, nbases=1e8)

# Loess-modified function
loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=16, ncol=length(qq))
  rownames(est) <- paste0(rep(c("A","C","G","T"), each=4), "2", rep(c("A","C","G","T"), 4))
  colnames(est) <- colnames(trans)
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        df <- data.frame(q=qq, rlogp=log10((errs+1)/(tot+1)))
        mod.lo <- loess(rlogp ~ q, df, span=0.95, degree=1)
        pred <- predict(mod.lo, qq)
        pred[is.na(pred)] <- 0
        est[paste0(nti,"2",ntj), ] <- 10^pred
      }
    }
  }
  est[est > 0.25] <- 0.25
  est[est < 1e-7] <- 1e-7
  return(est)
}

errF_loess <- learnErrors(filtFs, multithread=TRUE, nbases=1e8, errorEstimationFunction=loessErrfun_mod)
errR_loess <- learnErrors(filtRs, multithread=TRUE, nbases=1e8, errorEstimationFunction=loessErrfun_mod)

# Fit metric
error_fit_metric <- function(errObj) {
  obs <- errObj$err_out$obs
  exp <- errObj$err_out$est
  sum((obs - exp)^2, na.rm=TRUE)
}

fits <- c(
  default = error_fit_metric(errF_default) + error_fit_metric(errR_default),
  loess = error_fit_metric(errF_loess) + error_fit_metric(errR_loess)
)
best_model <- names(which.min(fits))
cat("Best error model chosen:", best_model, "\n")

if(best_model == "default") {
  errF <- errF_default
  errR <- errR_default
} else {
  errF <- errF_loess
  errR <- errR_loess
}

# ----------------------
# 6. Dereplicate and DADA2 denoise
# ----------------------
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

# ----------------------
# 7. Filter by expected Hsp60 length
# ----------------------
target_range <- 390:475
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% target_range]

# ----------------------
# 8. Remove chimeras
# ----------------------
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)

# ----------------------
# 9. Create ASV tables
# ----------------------
asv_seqs <- colnames(seqtab_nochim)
rownames(seqtab_nochim) <- sample_names

# Sequence-based ASVs
asv_tab_seq <- t(seqtab_nochim)
rownames(asv_tab_seq) <- asv_seqs
write.table(asv_tab_seq,
            paste0(output_prefix, "_Counts_seqASV_b.tsv"),
            sep="\t", quote=F, col.names=NA)

# Numbered ASVs
asv_ids <- paste0("ASV_", seq_along(asv_seqs))
asv_tab_num <- asv_tab_seq
rownames(asv_tab_num) <- asv_ids
write.table(asv_tab_num,
            paste0(output_prefix, "_Counts_numASV.tsv"),
            sep="\t", quote=F, col.names=NA)

# ----------------------
# 10. Write FASTA of numbered ASVs
# ----------------------
asv_fasta <- c(rbind(paste0(">", asv_ids), asv_seqs))
write(asv_fasta, paste0(output_prefix, "_ASVs.fa"))

# ----------------------
# 11. Run QIIME2 classifier from R
# ----------------------
cat("Running QIIME2 classification...\n")

# Define paths
asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
queries_qza <- paste0(output_prefix, "_ASVs.qza")
taxonomy_qza <- paste0(output_prefix, "_taxonomy.qza")
taxonomy_tsv <- paste0(output_prefix, "_taxonomy.tsv")
cpn60_classifier <- "/path/to/cpn60_classifier_v11.qza"  # <-- adjust this

# Import ASVs to QIIME2 artifact
system2("qiime", args=c("tools", "import",
                        "--type", "'FeatureData[Sequence]'",
                        "--input-path", asv_fasta_path,
                        "--output-path", queries_qza))

# Run classification
system2("qiime", args=c("feature-classifier", "classify-sklearn",
                        "--i-classifier", cpn60_classifier,
                        "--i-reads", queries_qza,
                        "--o-classification", taxonomy_qza))

# Export classification to TSV
taxonomy_export_dir <- paste0(output_prefix, "_taxonomy_export")
if(!dir.exists(taxonomy_export_dir)) dir.create(taxonomy_export_dir)

system2("qiime", args=c("tools", "export",
                        "--input-path", taxonomy_qza,
                        "--output-path", taxonomy_export_dir))

file.rename(file.path(taxonomy_export_dir, "taxonomy.tsv"), taxonomy_tsv)
cat("Taxonomy written to:", taxonomy_tsv, "\n")

# ----------------------
# 12. Merge counts + taxonomy
# ----------------------
counts_num <- asv_tab_num
taxonomy <- read.table(taxonomy_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)

# Ensure ASV IDs match
taxonomy <- taxonomy[rownames(counts_num), , drop=FALSE]

counts_tax <- cbind(counts_num, taxonomy)
write.table(counts_tax, paste0(output_prefix, "_ASVs_counts_taxonomy.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

cat("Pipeline finished. Outputs:\n")
cat("1. Sequence-based ASV table:", paste0(output_prefix, "_Counts_seqASV_b.tsv"), "\n")
cat("2. Numbered ASV table:", paste0(output_prefix, "_Counts_numASV.tsv"), "\n")
cat("3. FASTA for classifier:", asv_fasta_path, "\n")
cat("4. Taxonomy TSV:", taxonomy_tsv, "\n")
cat("5. Counts + taxonomy table:", paste0(output_prefix, "_ASVs_counts_taxonomy.tsv"), "\n")

# ----------------------
# 13. Create phyloseq object with numbered ASVs
# ----------------------
sample_metadata <- data.frame(
  SampleID = colnames(asv_tab_num),
  TotalReads = colSums(asv_tab_num),
  row.names = colnames(asv_tab_num),
  stringsAsFactors = FALSE
)

ps <- phyloseq(
  otu_table(as.matrix(asv_tab_num), taxa_are_rows=TRUE),
  sample_data(sample_metadata)
)

# Save phyloseq object
ps_file <- paste0(output_prefix, "_phyloseq.rds")
saveRDS(ps, ps_file)
cat("Phyloseq object written to:", ps_file, "\n")
