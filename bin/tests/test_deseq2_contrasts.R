#!/usr/bin/env Rscript
# F3.3 / F3.5 test: deseq2_contrasts.R on synthetic counts.
#   apptainer exec $R_BIOC_IMG Rscript bin/tests/test_deseq2_contrasts.R

suppressPackageStartupMessages({ library(DESeq2); library(readr) })

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_deseq2_contrasts.R"
script <- file.path(dirname(dirname(normalizePath(self))), "deseq2_contrasts.R")
if (!file.exists(script)) stop("cannot locate bin/deseq2_contrasts.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

set.seed(1)
ng <- 60
cond <- c("A", "A", "A", "B", "B", "B")
counts <- matrix(rnbinom(ng * 6, mu = 200, size = 5), nrow = ng,
                 dimnames = list(sprintf("G%03d", 1:ng), paste0("S", 1:6)))
# make the first 10 genes clearly DE
counts[1:10, cond == "B"] <- counts[1:10, cond == "B"] * 8

d <- tempfile("dx_"); dir.create(d)
saveRDS(counts, file.path(d, "counts.rds"))
coldata <- data.frame(sample = colnames(counts), condition = cond,
                      batch = "", norm_group = "G", stringsAsFactors = FALSE)
write.csv(coldata, file.path(d, "coldata.csv"), row.names = FALSE)
writeLines("A_vs_B\tA\tB\t", file.path(d, "contrasts.tsv"))

st <- system2("Rscript", c(script, "--group", "G", "--counts", file.path(d, "counts.rds"),
                           "--coldata", file.path(d, "coldata.csv"),
                           "--contrasts", file.path(d, "contrasts.tsv"),
                           "--design-formula", "~ condition", "--fdr", "0.05",
                           "--outdir", d))
check(st == 0, "deseq2_contrasts.R exits 0")

de_path <- file.path(d, "G_deseq2_results", "A_vs_B_DE.tsv")
check(file.exists(de_path), "DE table written")
if (file.exists(de_path)) {
  de <- read_tsv(de_path, show_col_types = FALSE)
  check(all(c("gene_id", "log2FoldChange", "padj") %in% colnames(de)),
        "DE table has gene_id/log2FoldChange/padj")
}

# F3.5: padj must match results() WITHOUT lfcThreshold (and differ from with).
cd <- coldata; rownames(cd) <- colnames(counts); cd$condition <- factor(cd$condition)
dds <- DESeq(DESeqDataSetFromMatrix(counts, cd, ~ condition))
ref_no <- as.data.frame(results(dds, contrast = c("condition", "A", "B"), alpha = 0.05))
ref_thr <- as.data.frame(results(dds, contrast = c("condition", "A", "B"),
                                 lfcThreshold = 1, alpha = 0.05))
if (file.exists(de_path)) {
  de <- as.data.frame(read_tsv(de_path, show_col_types = FALSE))
  rownames(de) <- de$gene_id
  common <- intersect(rownames(de), rownames(ref_no))
  same_no <- isTRUE(all.equal(de[common, "padj"], ref_no[common, "padj"], tolerance = 1e-8))
  diff_thr <- !isTRUE(all.equal(de[common, "padj"], ref_thr[common, "padj"], tolerance = 1e-8))
  check(same_no, "padj equals results() without lfcThreshold")
  check(diff_thr, "padj differs from results(lfcThreshold=1) -> cutoff not applied in DESeq2")
}

na_path <- file.path(d, "G_deseq2_results", "padj_na_summary.tsv")
check(file.exists(na_path), "padj_na_summary.tsv written")
if (file.exists(na_path)) {
  na <- read_tsv(na_path, show_col_types = FALSE)
  check(all(c("contrast", "n_tested", "n_padj_na", "n_padj_na_baseMean_zero",
              "n_padj_na_independent_filtering") %in% colnames(na)),
        "padj_na_summary has the F3.3 columns")
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL DESEQ2 CONTRASTS TESTS PASSED\n")
