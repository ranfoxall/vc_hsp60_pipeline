# vc_hsp60_classify.R
# Taxonomy classification for Vibrio-centric hsp60 ASVs
#
# Two classification methods:
#   vibrio  (default) — two-step approach matching King et al. 2019:
#                       1. BLAST screen against Vibrio hsp60 reference set
#                       2. Vibrio ASVs classified with Vibrio-only sklearn classifier
#                       3. Non-Vibrio ASVs classified with universal cpn60 classifier
#   sklearn — universal cpn60 sklearn classifier only
#
# Run AFTER vc_hsp60_pipeline.R.
# On UNH Premise, use the SLURM script (run_vchsp60_classify.slurm).
# Manually:
#   conda activate vc_hsp60_pipeline
#   Rscript vc_hsp60_classify.R \
#     --output_prefix Hsp60_MyRun \
#     --method vibrio \
#     --blast_db /path/to/repset_final_130219.fas \
#     --vibrio_classifier /path/to/vc_hsp60_classifier_sklearn022.qza \
#     --classifier /path/to/cpn60_classifier_v11.qza \
#     --qiime2_env qiime2-2020.2 \
#     --confidence 0.7
#
# UNH Premise reference files:
#   /mnt/home/whistler/foxall/hsp60_ref/vc_hsp60_ref/repset_final_130219.fas
#   /mnt/home/whistler/foxall/hsp60_ref/vc_hsp60_ref/vc_hsp60_classifier_sklearn022.qza
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
  result <- list(
    output_prefix     = NULL,
    method            = "vibrio",
    blast_db          = NULL,
    vibrio_classifier = NULL,
    classifier        = NULL,
    confidence        = "0.7",
    qiime2_env        = "qiime2-2020.2",
    phyloseq_rds      = NULL
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--output_prefix" && i < length(args)) {
      result$output_prefix <- args[i+1]; i <- i + 2
    } else if (args[i] == "--method" && i < length(args)) {
      result$method <- args[i+1]; i <- i + 2
    } else if (args[i] == "--blast_db" && i < length(args)) {
      result$blast_db <- args[i+1]; i <- i + 2
    } else if (args[i] == "--vibrio_classifier" && i < length(args)) {
      result$vibrio_classifier <- args[i+1]; i <- i + 2
    } else if (args[i] == "--classifier" && i < length(args)) {
      result$classifier <- args[i+1]; i <- i + 2
    } else if (args[i] == "--confidence" && i < length(args)) {
      result$confidence <- args[i+1]; i <- i + 2
    } else if (args[i] == "--qiime2_env" && i < length(args)) {
      result$qiime2_env <- args[i+1]; i <- i + 2
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
if (!opt$method %in% c("vibrio", "sklearn")) stop("--method must be 'vibrio' or 'sklearn'")

output_prefix     <- opt$output_prefix
method            <- opt$method
phyloseq_rds      <- if (!is.null(opt$phyloseq_rds)) opt$phyloseq_rds else paste0(output_prefix, "_phyloseq.rds")
confidence        <- opt$confidence
qiime2_env        <- opt$qiime2_env

cat("qiime2_env:    ", qiime2_env, "\n")

asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
counts_table   <- paste0(output_prefix, "_Counts_numASV.tsv")

for (f in c(asv_fasta_path, counts_table)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

cat("output_prefix: ", output_prefix, "\n")
cat("method:        ", method, "\n")
cat("confidence:    ", confidence, "\n")

# taxonomy rank prefixes — standard phyloseq convention
tax_ranks   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_prefix <- c("k__",     "p__",    "c__",   "o__",   "f__",    "g__",   "s__")

# helper: build tax_mat from semicolon-delimited taxonomy strings
build_tax_mat <- function(taxon_strings, asv_ids) {
  tax_split <- strsplit(taxon_strings, ";")
  tax_mat <- t(sapply(tax_split, function(x) {
    x <- x[1:min(7, length(x))]
    x <- mapply(function(val, pfx) {
      val <- trimws(val)
      if (is.na(val) || val == "" || val == "-") return(pfx)
      if (!grepl("^[kpcofgs]__", val)) val <- paste0(pfx, gsub(" ", "_", val))
      return(val)
    }, x, rank_prefix[1:length(x)], SIMPLIFY=TRUE)
    length(x) <- 7
    x
  }))
  rownames(tax_mat) <- asv_ids
  colnames(tax_mat) <- tax_ranks
  for (i in seq_along(tax_ranks)) {
    empty <- is.na(tax_mat[, i]) | tax_mat[, i] == "" | tax_mat[, i] == rank_prefix[i]
    tax_mat[empty, i] <- rank_prefix[i]
  }
  tax_mat
}

# helper: run QIIME2 sklearn classification on a FASTA file
# qiime must be available in PATH — activate your QIIME2 conda environment before running,
# or use the SLURM scripts which handle environment activation automatically.
run_sklearn <- function(fasta_path, classifier, confidence, prefix_tag) {
  queries_qza         <- paste0(prefix_tag, "_ASVs.qza")
  taxonomy_qza        <- paste0(prefix_tag, "_taxonomy.qza")
  taxonomy_tsv        <- paste0(prefix_tag, "_taxonomy.tsv")
  taxonomy_export_dir <- paste0(prefix_tag, "_taxonomy_export")

  ret <- system2("qiime", args=c(
    "tools", "import",
    "--type", "FeatureData[Sequence]",
    "--input-path",  fasta_path,
    "--output-path", queries_qza
  ))
  if (ret != 0) stop(
    "QIIME2 import failed for: ", fasta_path, "\n",
    "Make sure your QIIME2 environment is active and qiime is in PATH."
  )

  ret <- system2("qiime", args=c(
    "feature-classifier", "classify-sklearn",
    "--i-classifier",     classifier,
    "--i-reads",          queries_qza,
    "--p-confidence",     confidence,
    "--o-classification", taxonomy_qza
  ))
  if (ret != 0) stop(
    "Classification failed for: ", fasta_path, "\n",
    "Common cause: scikit-learn version mismatch between classifier and QIIME2 environment.\n",
    "See README for retraining instructions."
  )

  if (!dir.exists(taxonomy_export_dir)) dir.create(taxonomy_export_dir)
  system2("qiime", args=c(
    "tools", "export",
    "--input-path",  taxonomy_qza,
    "--output-path", taxonomy_export_dir
  ))
  file.rename(file.path(taxonomy_export_dir, "taxonomy.tsv"), taxonomy_tsv)

  read.table(taxonomy_tsv, sep="\t", header=TRUE,
             stringsAsFactors=FALSE, row.names=1)
}

# =============================================================================
# VIBRIO METHOD (default) — BLAST filter + Vibrio sklearn + universal sklearn
# =============================================================================
if (method == "vibrio") {

  blast_db             <- opt$blast_db
  vibrio_classifier    <- opt$vibrio_classifier
  universal_classifier <- opt$classifier

  if (is.null(blast_db))          stop("Please provide --blast_db for vibrio method")
  if (is.null(vibrio_classifier)) stop("Please provide --vibrio_classifier for vibrio method")
  if (!file.exists(blast_db))          stop("BLAST database FASTA not found: ", blast_db)
  if (!file.exists(vibrio_classifier)) stop("Vibrio classifier not found: ", vibrio_classifier)

  cat("blast_db:           ", blast_db, "\n")
  cat("vibrio_classifier:  ", vibrio_classifier, "\n")
  if (!is.null(universal_classifier))
    cat("universal_classifier:", universal_classifier, "\n")

  # check BLAST
  blast_check <- system2("which", args="blastn", stdout=TRUE, stderr=TRUE)
  if (length(blast_check) == 0 || !grepl("blastn", blast_check[1])) {
    stop("blastn not found in PATH.\nActivate vc_hsp60_pipeline conda environment.")
  }
  cat("blastn found at:", blast_check[1], "\n")

  # build BLAST database if needed
  db_prefix <- sub("\\.fas$", "", blast_db)
  if (!file.exists(paste0(db_prefix, ".nhr"))) {
    cat("Building BLAST database...\n")
    ret <- system2("makeblastdb", args=c(
      "-in", blast_db, "-dbtype", "nucl", "-out", db_prefix
    ))
    if (ret != 0) stop("makeblastdb failed.")
    cat("BLAST database built.\n")
  } else {
    cat("BLAST database already exists — skipping makeblastdb.\n")
  }

  # run BLAST — screen for Vibrio sequences (>=90% identity, >=90% coverage)
  blast_out <- paste0(output_prefix, "_blast_screen.tsv")
  cat("Running BLAST screen...\n")
  blast_cmd <- paste0(
    "blastn",
    " -query ",   asv_fasta_path,
    " -db ",      db_prefix,
    " -out ",     blast_out,
    " -outfmt '6 qseqid sseqid pident qcovs evalue bitscore'",
    " -perc_identity 90",
    " -max_target_seqs 1",
    " -num_threads 4"
  )
  ret <- system(blast_cmd)
  if (ret != 0) stop("BLAST failed.")
  cat("BLAST screen complete.\n")

  # get all ASV IDs from FASTA
  fasta_lines <- readLines(asv_fasta_path)
  all_asvs    <- sub("^>", "", fasta_lines[grepl("^>", fasta_lines)])

  # determine which ASVs are Vibrio (have a BLAST hit)
  blast_cols <- c("qseqid", "sseqid", "pident", "qcovs", "evalue", "bitscore")
  if (file.exists(blast_out) && file.size(blast_out) > 0) {
    blast_raw   <- read.table(blast_out, sep="\t", header=FALSE,
                               col.names=blast_cols, stringsAsFactors=FALSE)
    # keep only hits with qcovs >= 90
    blast_raw   <- blast_raw[blast_raw$qcovs >= 90, ]
    vibrio_asvs <- unique(blast_raw$qseqid)
  } else {
    blast_raw   <- data.frame()
    vibrio_asvs <- character(0)
  }
  non_vibrio_asvs <- all_asvs[!all_asvs %in% vibrio_asvs]

  cat("Vibrio ASVs:     ", length(vibrio_asvs), "\n")
  cat("Non-Vibrio ASVs: ", length(non_vibrio_asvs), "\n")

  # write BLAST screen summary
  blast_screen_file <- paste0(output_prefix, "_blast_screen.csv")
  screen_summary <- data.frame(
    ASV_ID   = all_asvs,
    is_vibrio = all_asvs %in% vibrio_asvs,
    stringsAsFactors = FALSE
  )
  if (nrow(blast_raw) > 0) {
    screen_summary <- merge(screen_summary,
                             blast_raw[, c("qseqid","sseqid","pident","qcovs")],
                             by.x="ASV_ID", by.y="qseqid", all.x=TRUE)
  }
  write.csv(screen_summary, blast_screen_file, row.names=FALSE, quote=FALSE)
  cat("BLAST screen results written to:", blast_screen_file, "\n")

  # split FASTA into Vibrio and non-Vibrio
  fasta_seqs <- split(fasta_lines, cumsum(grepl("^>", fasta_lines)))

  all_tax_raw <- NULL

  # classify Vibrio ASVs with Vibrio-only classifier
  if (length(vibrio_asvs) > 0) {
    vibrio_fasta_path <- paste0(output_prefix, "_ASVs_vibrio.fa")
    vibrio_fasta <- unlist(lapply(fasta_seqs, function(x) {
      id <- sub("^>", "", x[1])
      if (id %in% vibrio_asvs) x else NULL
    }))
    writeLines(vibrio_fasta, vibrio_fasta_path)
    cat("Classifying", length(vibrio_asvs), "Vibrio ASVs...\n")
    vibrio_tax  <- run_sklearn(vibrio_fasta_path, vibrio_classifier,
                                confidence, paste0(output_prefix, "_vibrio"))
    all_tax_raw <- vibrio_tax
    cat("Vibrio classification complete.\n")
  }

  # classify non-Vibrio ASVs with universal classifier
  if (length(non_vibrio_asvs) > 0) {
    if (!is.null(universal_classifier) && file.exists(universal_classifier)) {
      non_vibrio_fasta_path <- paste0(output_prefix, "_ASVs_non_vibrio.fa")
      non_vibrio_fasta <- unlist(lapply(fasta_seqs, function(x) {
        id <- sub("^>", "", x[1])
        if (id %in% non_vibrio_asvs) x else NULL
      }))
      writeLines(non_vibrio_fasta, non_vibrio_fasta_path)
      cat("Classifying", length(non_vibrio_asvs), "non-Vibrio ASVs with universal classifier...\n")
      non_vibrio_tax <- run_sklearn(non_vibrio_fasta_path, universal_classifier,
                                     confidence, paste0(output_prefix, "_non_vibrio"))
      all_tax_raw    <- rbind(all_tax_raw, non_vibrio_tax)
      cat("Universal classification complete.\n")
    } else {
      cat("No universal classifier provided — non-Vibrio ASVs will be unassigned.\n")
      unassigned <- data.frame(
        Taxon      = rep("k__;p__;c__;o__;f__;g__;s__", length(non_vibrio_asvs)),
        Confidence = rep(0, length(non_vibrio_asvs)),
        row.names  = non_vibrio_asvs,
        stringsAsFactors = FALSE
      )
      all_tax_raw <- rbind(all_tax_raw, unassigned)
    }
  }

  # reorder to match original ASV order
  present_asvs <- all_asvs[all_asvs %in% rownames(all_tax_raw)]
  all_tax_raw  <- all_tax_raw[present_asvs, , drop=FALSE]

  taxon_strings   <- all_tax_raw$Taxon
  tax_mat         <- build_tax_mat(taxon_strings, rownames(all_tax_raw))
  confidence_file <- paste0(output_prefix, "_taxonomy_confidence.csv")
  write.csv(data.frame(Confidence=all_tax_raw$Confidence,
                        row.names=rownames(all_tax_raw)),
            confidence_file, quote=FALSE)
  cat("Confidence scores written to:", confidence_file, "\n")

# =============================================================================
# SKLEARN METHOD — universal classifier only
# =============================================================================
} else {

  cpn60_classifier <- opt$classifier
  if (is.null(cpn60_classifier)) stop("Please provide --classifier for sklearn method")
  if (!file.exists(cpn60_classifier)) stop("Classifier not found: ", cpn60_classifier)
  cat("classifier:    ", cpn60_classifier, "\n")

  tax_raw <- run_sklearn(asv_fasta_path, cpn60_classifier, confidence, output_prefix)
  tax_mat <- build_tax_mat(tax_raw$Taxon, rownames(tax_raw))

  confidence_file <- paste0(output_prefix, "_taxonomy_confidence.csv")
  write.csv(data.frame(Confidence=tax_raw$Confidence, row.names=rownames(tax_raw)),
            confidence_file, quote=FALSE)
  cat("Confidence scores written to:", confidence_file, "\n")
}

# =============================================================================
# SHARED — taxonomy table, counts merge, summary
# =============================================================================

taxonomy_table      <- as.data.frame(tax_mat, stringsAsFactors=FALSE)
taxonomy_table_file <- paste0(output_prefix, "_taxonomy_table.csv")
write.csv(taxonomy_table, taxonomy_table_file, quote=FALSE)
cat("Taxonomy table written to:", taxonomy_table_file, "\n")

counts_num    <- read.table(counts_table, sep="\t", header=TRUE,
                            stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
tax_for_merge <- taxonomy_table[rownames(counts_num), , drop=FALSE]
counts_tax      <- cbind(counts_num, tax_for_merge)
counts_tax_file <- paste0(output_prefix, "_ASVs_counts_taxonomy.tsv")
write.table(counts_tax, counts_tax_file, sep="\t", quote=FALSE, col.names=NA)
cat("Counts + taxonomy written to:", counts_tax_file, "\n")

cat("\nClassification complete. Outputs:\n")
cat("  ", taxonomy_table_file, "\n")
cat("  ", confidence_file, "\n")
cat("  ", counts_tax_file, "\n")
if (method == "vibrio") cat("  ", blast_screen_file, "\n")
cat("\nNext step: vc_hsp60_add_taxonomy.R will merge taxonomy and refseq into phyloseq RDS.\n")
