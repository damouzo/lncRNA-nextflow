#!/usr/bin/env Rscript
# F3.1 test: cis_associations.R, both normalization paths.
#   apptainer exec $R_BIOC_IMG Rscript bin/tests/test_cis_associations.R

suppressPackageStartupMessages(library(readr))

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_cis_associations.R"
script <- file.path(dirname(dirname(normalizePath(self))), "cis_associations.R")
if (!file.exists(script)) stop("cannot locate bin/cis_associations.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

d <- tempfile("cis_"); dir.create(d)
write_tsv(data.frame(
  gene_id = c("L1", "L2", "P1", "P2"),
  chromosome = "1", start = c(1000, 50000, 2500, 52000), end = c(2000, 51000, 3500, 53000),
  strand = c("+", "-", "+", "-"),
  gene_biotype = c(NA, NA, "protein_coding", "protein_coding"),
  origin = c("novel", "annotated_lncrna", "reference_other", "reference_other"),
  catalog_id = c("LNC.1", "LNC.2", NA, NA),
  class_code = c("u", "NA", NA, NA),
  gene_name = c(NA, NA, "P1", "P2"),
  stringsAsFactors = FALSE), file.path(d, "gene_catalog.tsv"))

set.seed(2)
counts <- matrix(rnbinom(4 * 6, mu = 300, size = 5), nrow = 4,
                 dimnames = list(c("L1", "L2", "P1", "P2"), paste0("S", 1:6)))
counts["P1", ] <- counts["L1", ] + 5   # correlated pair
saveRDS(counts, file.path(d, "counts.rds"))
write.csv(data.frame(sample = paste0("S", 1:6),
                     condition = c("A", "A", "A", "B", "B", "B"),
                     batch = c("b1", "b1", "b2", "b1", "b2", "b2")),
          file.path(d, "coldata.csv"), row.names = FALSE)

for (norm in c("false", "true")) {
  o <- file.path(d, paste0("out_", norm)); dir.create(o)
  st <- system2("Rscript", c(script, "--group", "G", "--gene-catalog", file.path(d, "gene_catalog.tsv"),
                             "--counts", file.path(d, "counts.rds"),
                             "--coldata", file.path(d, "coldata.csv"),
                             "--cis-window", "100000", "--fdr", "0.05",
                             "--normalize", norm, "--outdir", o))
  check(st == 0, sprintf("cis_associations.R exits 0 (normalize=%s)", norm))
  pairs <- file.path(o, "G_cis_pairs.tsv")
  summ <- file.path(o, "G_cis_associations_summary.tsv")
  check(file.exists(pairs) && file.exists(summ),
        sprintf("outputs written (normalize=%s)", norm))
  if (file.exists(pairs)) {
    p <- read_tsv(pairs, show_col_types = FALSE)
    check(all(c("lncrna_id", "target_gene_id", "spearman_rho", "padj") %in% colnames(p)),
          sprintf("cis table keeps its columns (normalize=%s)", norm))
  }
  if (file.exists(summ)) {
    hdr <- strsplit(readLines(summ, n = 1), "\t")[[1]]
    check(identical(hdr, c("contrast", "n_pairs", "n_significant", "median_cor")),
          sprintf("summary format unchanged (normalize=%s)", norm))
  }
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL CIS ASSOCIATION TESTS PASSED\n")
