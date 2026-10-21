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
for (s in samples) {
  qdir <- file.path(work, paste0(s, "_quant"))
  dir.create(qdir)
  num <- ifelse(1:n_genes <= kept_genes, 100, 1)
  qsf <- data.frame(Name = tx$TXNAME, Length = 500,
                    EffectiveLength = 450, TPM = num * 10, NumReads = num)
  write_tsv(qsf, file.path(qdir, "quant.sf"))
}

coldata <- data.frame(sample = samples,
                      condition = c("DDX41", "DDX41", "Healthy", "Healthy"),
                      norm_group = "MNC")
write.csv(coldata, file.path(work, "coldata_MNC.csv"), row.names = FALSE)

old <- getwd()
setwd(work)
status <- system2("Rscript",
                  c(script, "--group", "MNC", "--coldata", "coldata_MNC.csv",
                    "--tx2gene", "tx2gene.tsv", "--outdir", "."))
setwd(old)
check(status == 0, "tximport_filter.R exits 0")

for (f in c("MNC_counts.tsv", "MNC_tpm.tsv", "MNC_universe.tsv",
            "MNC_filter_stats.tsv", "MNC_tximport_counts.rds")) {
  check(file.exists(file.path(work, f)), paste("wrote", f))
}

if (file.exists(file.path(work, "MNC_universe.tsv"))) {
  u <- read_tsv(file.path(work, "MNC_universe.tsv"), show_col_types = FALSE)
  check(nrow(u) == n_genes, "universe keeps every imported gene")
  check(sum(u$keep) == kept_genes, "universe keep count matches the funnel")
  check(all(u$filter_reason[u$keep] == ""), "kept genes have no filter_reason")
  check(all(u$filter_reason[!u$keep] == "low_expression"),
        "filtered genes say low_expression")
}

if (file.exists(file.path(work, "MNC_counts.tsv"))) {
  cnt <- read_tsv(file.path(work, "MNC_counts.tsv"), show_col_types = FALSE)
  check(ncol(cnt) == length(samples) + 1, "counts.tsv is gene_id + one column per sample")
  check(nrow(cnt) == kept_genes, "counts.tsv holds the filtered universe")
}

cat("\n")
if (fails > 0) {
  cat(sprintf("FAILED: %d\n", fails))
  quit(save = "no", status = 1)
}
cat("ALL TXIMPORT FILTER TESTS PASSED\n")
