#!/usr/bin/env Rscript
# Test: contrast_summary_plots.R writes per-subset plots + tables without failing.
#   apptainer exec $REPORTING_IMG Rscript bin/tests/test_contrast_summary_plots.R

suppressPackageStartupMessages(library(readr))

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_contrast_summary_plots.R"
script <- file.path(dirname(dirname(normalizePath(self))), "contrast_summary_plots.R")
if (!file.exists(script)) stop("cannot locate bin/contrast_summary_plots.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

d <- tempfile("ccp_"); dir.create(d)

# 6 contrasts (> VENN_MAX_SETS) so the Venn path takes the text-placeholder branch
# without needing ggvenn; upsets are guarded so their absence is not an error.
cts <- sprintf("grp__c%d", 1:6)
ac <- data.frame(
  gene_id = sprintf("G%02d", 1:8),
  gene_id_stable = sprintf("G%02d", 1:8),
  gene_name = sprintf("n%02d", 1:8),
  gene_group = c("novel_lncRNA", "novel_lncRNA", "annotated_lncRNA", "annotated_lncRNA",
                 "novel_lncRNA", "annotated_lncRNA", "novel_lncRNA", "annotated_lncRNA"),
  origin = c("novel", "novel", "annotated_lncrna", "annotated_lncrna",
             "novel", "annotated_lncrna", "novel", "annotated_lncrna"),
  tier = "1a", lnc_position_class = "antisense", n_contrasts_significant = 0,
  stringsAsFactors = FALSE)
for (i in seq_along(cts)) {
  up   <- c(0.001, 0.01, 0.5, NA, 0.2, 0.03, NA, 0.4)[(seq_len(8) + i) %% 8 + 1]
  lfc  <- c(-2.5, 1.8, -0.4, NA, 2.2, -1.5, NA, 0.3)[(seq_len(8) + i) %% 8 + 1]
  ac[[paste0(cts[i], ".log2FoldChange")]] <- lfc
  ac[[paste0(cts[i], ".padj")]] <- up
}
write_tsv(ac, file.path(d, "lncrna_across_contrasts.tsv"))

o <- file.path(d, "out"); dir.create(o)
st <- system2("Rscript", c(script, "--across", file.path(d, "lncrna_across_contrasts.tsv"),
                           "--outdir", o))
check(st == 0, "contrast_summary_plots.R exits 0 (non-blocking)")

for (sub in c("combined", "annotated", "novel")) {
  check(file.exists(file.path(o, "plots", sub, "sig_counts_by_contrast.png")),
        sprintf("plots/%s/sig_counts_by_contrast.png written", sub))
  check(file.exists(file.path(o, "plots", sub, "sig_counts_by_contrast.pdf")),
        sprintf("plots/%s/sig_counts_by_contrast.pdf written", sub))
  check(file.exists(file.path(o, "tables", sprintf("membership_%s.tsv", sub))),
        sprintf("tables/membership_%s.tsv written", sub))
  # 6 contrasts > VENN_MAX_SETS: the Venn is a text placeholder (no ggvenn needed).
  check(file.exists(file.path(o, "plots", sub, "venn_lncrna.png")),
        sprintf("plots/%s/venn_lncrna.png written (placeholder)", sub))
}
check(file.exists(file.path(o, "tables", "lncrna_across_contrasts.tsv")),
      "tables/lncrna_across_contrasts.tsv written")

if (file.exists(file.path(o, "tables", "membership_novel.tsv"))) {
  m <- read_tsv(file.path(o, "tables", "membership_novel.tsv"), show_col_types = FALSE)
  check(all(c("gene_id", "origin", "n_contrasts_significant") %in% colnames(m)),
        "membership table carries gene_id/origin/n_contrasts_significant")
  check(all(m$origin == "novel"), "novel membership only contains novel lncRNAs")
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL CONTRAST SUMMARY PLOT TESTS PASSED\n")
