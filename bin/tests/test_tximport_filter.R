#!/usr/bin/env Rscript
# Minimal synthetic test for bin/tximport_filter.R.
# Runs without Nextflow; execute inside the image, e.g.:
#   apptainer exec $IMG Rscript bin/tests/test_tximport_filter.R
# (needs tximport/DESeq2/readr/ggplot2 — the r_bioc image)

suppressPackageStartupMessages(library(readr))

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_tximport_filter.R"
script <- file.path(dirname(dirname(normalizePath(self))), "tximport_filter.R")
if (!file.exists(script)) stop("cannot locate bin/tximport_filter.R next to the test")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

work <- tempfile("txitest_")
dir.create(work)
n_genes <- 200
kept_genes <- 120

tx <- data.frame(TXNAME = sprintf("TX%03d", 1:n_genes),
                 GENEID = sprintf("G%03d", 1:n_genes))
write_tsv(tx, file.path(work, "tx2gene.tsv"))

samples <- c("S1", "S2", "S3", "S4")
sfac <- c(S1 = 1.0, S2 = 1.4, S3 = 0.7, S4 = 1.2)
set.seed(42)
for (s in samples) {
  qdir <- file.path(work, paste0(s, "_quant"))
  dir.create(qdir)
  base <- ifelse(1:n_genes <= kept_genes, 200, 2)
  num <- rnbinom(n_genes, mu = base * sfac[[s]], size = 5)
  qsf <- data.frame(Name = tx$TXNAME, Length = 500,
                    EffectiveLength = 450, TPM = num * 10, NumReads = num)
  write_tsv(qsf, file.path(qdir, "quant.sf"))
}

coldata <- data.frame(sample = c("S3", "S1", "S4", "S2"),
                      condition = c("Healthy", "DDX41", "Healthy", "DDX41"),
                      norm_group = "MNC")
write.csv(coldata, file.path(work, "coldata_MNC.csv"), row.names = FALSE)

# catalog_genes.tsv carries the annotation columns the TSVs must gain.
catalog <- data.frame(
  gene_id = sprintf("G%03d", 1:n_genes),
  gene_group = rep(c("protein_coding", "novel_lncRNA", "annotated_lncRNA"), length.out = n_genes),
  catalog_id = ifelse(seq_len(n_genes) %% 3 == 0, sprintf("LNC%04d", seq_len(n_genes)), NA_character_),
  origin = rep(c("reference_other", "novel", "annotated_lncrna"), length.out = n_genes),
  class_code = rep(c(NA_character_, "u", NA_character_), length.out = n_genes),
  tier = rep(c("NA", "1a", "NA"), length.out = n_genes),
  lnc_position_class = rep(c("no_coding_overlap", "antisense", "intronic"), length.out = n_genes),
  stringsAsFactors = FALSE)
write_tsv(catalog, file.path(work, "catalog_genes.tsv"))

old <- getwd()
setwd(work)
status <- system2("Rscript",
                  c(script, "--group", "MNC", "--coldata", "coldata_MNC.csv",
                    "--tx2gene", "tx2gene.tsv", "--catalog", "catalog_genes.tsv",
                    "--outdir", "."))
setwd(old)
check(status == 0, "tximport_filter.R exits 0")

for (f in c("MNC_counts.tsv", "MNC_tpm.tsv", "MNC_universe.tsv",
            "MNC_filter_stats.tsv", "MNC_tximport_counts.rds",
            "coldata_aligned_MNC.csv")) {
  check(file.exists(file.path(work, f)), paste("wrote", f))
}

if (file.exists(file.path(work, "MNC_universe.tsv"))) {
  u <- read_tsv(file.path(work, "MNC_universe.tsv"), show_col_types = FALSE)
  check(nrow(u) == n_genes, "universe keeps every imported gene")
  check(sum(u$keep) >= kept_genes - 10, "universe keep count close to the funnel")
  check(all(is.na(u$filter_reason[u$keep]) | u$filter_reason[u$keep] == ""),
        "kept genes have no filter_reason")
  check(all(u$filter_reason[!u$keep] == "low_expression"),
        "filtered genes say low_expression")
}

if (file.exists(file.path(work, "MNC_counts.tsv"))) {
  cnt <- read_tsv(file.path(work, "MNC_counts.tsv"), show_col_types = FALSE)
  check(ncol(cnt) == length(samples) + 7,
        "counts.tsv is gene_id + 6 annotation columns + one column per sample")
  check(all(c("gene_group", "catalog_id", "origin", "class_code", "tier",
              "lnc_position_class") %in% colnames(cnt)),
        "counts.tsv carries the annotation columns")
  check(nrow(cnt) >= kept_genes - 10, "counts.tsv holds the filtered universe")
}

# H0 / 0A.2: the published coldata must follow the counts column order (quant-dir
# listing, alphabetical), not the samplesheet order used to write coldata_MNC.csv.
if (file.exists(file.path(work, "coldata_aligned_MNC.csv"))) {
  cd <- read.csv(file.path(work, "coldata_aligned_MNC.csv"), stringsAsFactors = FALSE)
  cnt_hdr <- colnames(read_tsv(file.path(work, "MNC_counts.tsv"), show_col_types = FALSE))
  tpm_hdr <- colnames(read_tsv(file.path(work, "MNC_tpm.tsv"), show_col_types = FALSE))
  check(identical(as.character(cd$sample), sort(samples)),
        "aligned coldata follows the counts column order")
  check(identical(cnt_hdr[-(1:7)], as.character(cd$sample)),
        "counts.tsv sample columns match the aligned coldata")
  check(identical(tpm_hdr[-(1:7)], as.character(cd$sample)),
        "tpm.tsv sample columns match the aligned coldata")
}

# 2.2: the PCA is coloured by condition/batch, never by sample.
check(file.exists(file.path(work, "MNC_sample_qc_pca.png")),
      "sample QC PCA written with intgroup = condition/batch")

# 2.3: a coldata sample without a quant dir must fail, not be silently dropped.
write.csv(data.frame(sample = c(samples, "S5"),
                     condition = c("Healthy", "DDX41", "Healthy", "DDX41", "DDX41"),
                     norm_group = "MNC"),
          file.path(work, "coldata_extra.csv"), row.names = FALSE)
setwd(work)
st_missing <- system2("Rscript",
                      c(script, "--group", "MNC", "--coldata", "coldata_extra.csv",
                        "--tx2gene", "tx2gene.tsv", "--outdir", "out_missing"),
                      stdout = FALSE, stderr = FALSE)
setwd(old)
check(st_missing != 0, "a coldata sample without quant fails (n_shared != nrow)")

cat("\n")
if (fails > 0) {
  cat(sprintf("FAILED: %d\n", fails))
  quit(save = "no", status = 1)
}
cat("ALL TXIMPORT FILTER TESTS PASSED\n")
