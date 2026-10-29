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
  lnc_position_class = c("antisense", "intronic", NA, NA),
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
    check("lncrna_position_class" %in% colnames(p),
          sprintf("lncrna_position_class column present (normalize=%s)", norm))
    # 4a.1: the p-value is a real t-approximation with recorded method, never the
    # masked .Machine$double.xmin floor.
    check(all(c("spearman_pvalue", "df_effective", "p_method") %in% colnames(p)),
          sprintf("4a.1 p_method/df_effective columns present (normalize=%s)", norm))
    check(all(c("loo_min_abs_rho", "loo_sign_stable", "loo_max_delta") %in% colnames(p)),
          sprintf("4a.1b leave-one-out stability columns present (normalize=%s)", norm))
    if (nrow(p) > 0) {
      check(!any(is.finite(p$spearman_pvalue) & p$spearman_pvalue <= .Machine$double.xmin),
            sprintf("4a.1 no p-value clamped to double.xmin (normalize=%s)", norm))
      check(all(is.na(p$spearman_pvalue) | p$spearman_pvalue > 1e-12),
            sprintf("4a.1 p-values are realistic for n=6, not 2e-308 (normalize=%s)", norm))
    }
  }
  if (file.exists(summ)) {
    hdr <- strsplit(readLines(summ, n = 1), "\t")[[1]]
    check(identical(hdr, c("contrast", "n_pairs", "n_significant", "median_cor",
                           "correlation_input", "min_samples", "min_abs_rho")),
          sprintf("summary carries the effective thresholds (normalize=%s)", norm))
  }
}

# 2.4: two loci that track the condition only must not look associated once the
# within-condition residuals are used.
d2 <- tempfile("cis24_"); dir.create(d2)
write_tsv(data.frame(
  gene_id = c("L1", "P1"), chromosome = "1", start = c(1000, 1500), end = c(2000, 2500),
  strand = c("+", "-"), gene_biotype = c(NA, "protein_coding"),
  origin = c("novel", "reference_other"), catalog_id = c("LNC.1", NA),
  class_code = c("u", NA), lnc_position_class = c("antisense", NA),
  gene_name = c(NA, "P1"), stringsAsFactors = FALSE),
  file.path(d2, "gene_catalog.tsv"))
cnt <- matrix(c(c(10, 12, 11, 200, 210, 190),
                c(20, 18, 22, 300, 290, 310)), nrow = 2, byrow = TRUE,
              dimnames = list(c("L1", "P1"), paste0("S", 1:6)))
saveRDS(cnt, file.path(d2, "counts.rds"))
write.csv(data.frame(sample = paste0("S", 1:6),
                     condition = c("A", "A", "A", "B", "B", "B"), batch = ""),
          file.path(d2, "coldata.csv"), row.names = FALSE)

rho_of <- function(use_coldata) {
  o <- file.path(d2, paste0("o_", use_coldata)); dir.create(o, showWarnings = FALSE)
  args <- c(script, "--group", "G", "--gene-catalog", file.path(d2, "gene_catalog.tsv"),
            "--counts", file.path(d2, "counts.rds"), "--cis-window", "100000",
            "--fdr", "0.05", "--normalize", "false", "--outdir", o)
  if (use_coldata) args <- c(args, "--coldata", file.path(d2, "coldata.csv"))
  st <- system2("Rscript", args)
  check(st == 0, sprintf("2.4 run exits 0 (coldata=%s)", use_coldata))
  p <- file.path(o, "G_cis_pairs.tsv")
  if (!file.exists(p)) return(NA_real_)
  r <- read_tsv(p, show_col_types = FALSE)
  row <- r[which(r$target_gene_id == "P1"), , drop = FALSE]
  if (nrow(row) == 0) return(0)
  abs(as.numeric(row$spearman_rho[1]))
}
rho_expr <- rho_of(FALSE)
rho_resid <- rho_of(TRUE)
cat(sprintf("2.4 |rho|: expression=%.3f residuals=%.3f\n", rho_expr, rho_resid))
check(rho_expr > 0.5, "expression correlation is high (shared condition effect)")
check(is.na(rho_resid) || rho_resid < 0.5 || rho_resid < rho_expr,
      "residual correlation is not the condition artefact")

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL CIS ASSOCIATION TESTS PASSED\n")
