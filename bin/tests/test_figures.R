#!/usr/bin/env Rscript
# F5.1 test: figures.R writes PNGs + figures_manifest.tsv without failing.
#   apptainer exec $REPORTING_IMG Rscript bin/tests/test_figures.R

suppressPackageStartupMessages(library(readr))

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_figures.R"
script <- file.path(dirname(dirname(normalizePath(self))), "figures.R")
if (!file.exists(script)) stop("cannot locate bin/figures.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

d <- tempfile("fig_"); dir.create(d)
write_tsv(data.frame(
  gene_id = sprintf("G%02d", 1:10),
  origin = rep(c("novel", "annotated_lncrna"), 5),
  class_code = rep(c("u", "j", "x", "o", "i"), 2),
  chromosome = rep(c("chr1", "chr2"), 5),
  lnc_position_class = rep(c("antisense", "intronic", "overlapping", "intergenic"), length.out = 10),
  n_tier1_transcripts = c(1, 2, 1, 3, 1, 2, 2, 1, 1, 4),
  stringsAsFactors = FALSE), file.path(d, "catalog_genes.tsv"))
write_tsv(data.frame(
  tier1_allow_single_exon = c(TRUE, FALSE),
  unit = c("gene", "gene"),
  class_set = c("u", "u"),
  n = c(10, 20), stringsAsFactors = FALSE), file.path(d, "tier_final_by_class.tsv"))
write_tsv(data.frame(
  gene_id = sprintf("G%02d", 1:10),
  S1 = rep(c(5, 0), 5), S2 = rep(c(5, 0), 5),
  stringsAsFactors = FALSE), file.path(d, "grpA_tpm.tsv"))
write_tsv(data.frame(
  sample = c("S1", "S2"), norm_group = c("grpA", "grpA"),
  condition = c("case", "control"), stringsAsFactors = FALSE),
  file.path(d, "sample_metadata.tsv"))

o <- file.path(d, "out"); dir.create(o)
st <- system2("Rscript", c(script, "--catalog", file.path(d, "catalog_genes.tsv"),
                           "--by-class", file.path(d, "tier_final_by_class.tsv"),
                           "--tpm", file.path(d, "grpA_tpm.tsv"),
                           "--sample-metadata", file.path(d, "sample_metadata.tsv"),
                           "--outdir", o))
check(st == 0, "figures.R exits 0 (non-blocking)")
manifest <- file.path(o, "figures_manifest.tsv")
check(file.exists(manifest), "figures_manifest.tsv written")
if (file.exists(manifest)) {
  m <- read_tsv(manifest, show_col_types = FALSE)
  check(all(c("figure", "path", "description", "status") %in% colnames(m)),
        "manifest has figure/path/description/status")
  check(any(m$status == "ok"), "at least one figure produced")
  for (fig in c("fig_catalog_chromosome", "fig_lnc_position_class", "fig_lnc_position_expr")) {
    row <- m[m$figure == fig, , drop = FALSE]
    check(nrow(row) == 1 && row$status[1] == "ok", sprintf("%s produced", fig))
  }
  pngs <- list.files(o, pattern = "\\.png$")
  check(length(pngs) >= 1, "at least one PNG written")
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL FIGURE TESTS PASSED\n")
