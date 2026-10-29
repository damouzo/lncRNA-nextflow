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
                           "--design-formula", shQuote("~ condition"), "--fdr", "0.05", "--allow-no-offset", "true",
                           "--outdir", d))
check(st == 0, "deseq2_contrasts.R exits 0")

de_path <- file.path(d, "G_deseq2_results", "A_vs_B_DE.tsv")
check(file.exists(de_path), "DE table written")
if (file.exists(de_path)) {
  de <- read_tsv(de_path, show_col_types = FALSE)
  check(all(c("gene_id", "log2FoldChange", "padj") %in% colnames(de)),
        "DE table has gene_id/log2FoldChange/padj")
  # H4 / 1.2: unshrunken statistics + shrink type are published next to the LFC.
  check(all(c("stat", "lfc_unshrunken", "shrink_type") %in% colnames(de)),
        "DE table has stat/lfc_unshrunken/shrink_type")
  tested <- !is.na(de$pvalue)
  check(any(tested) && all(is.finite(de$stat[tested])),
        "stat is finite for the tested genes")
  check(all(de$shrink_type %in% c("apeglm", "ashr")),
        "shrink_type records which shrinkage ran")
  if (any(tested)) {
    dt <- de[tested, , drop = FALSE]
    check(!isTRUE(all.equal(dt$log2FoldChange, dt$lfc_unshrunken)),
          "shrunken LFC differs from the unshrunken LFC")
  }
}

# F3.5: padj must match results() WITHOUT lfcThreshold (and differ from with).
# Same relevel as the script (ref = denominator "B"): dispersion estimation is
# not bit-identical across reference levels, so the comparison would otherwise
# differ at ~1e-6 for reasons unrelated to lfcThreshold.
cd <- coldata; rownames(cd) <- colnames(counts)
cd$condition <- relevel(factor(cd$condition), ref = "B")
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

# H0 / 0A.1: coldata permuted (and inverted) relative to the counts columns must
# give byte-identical results, because the script aligns by sample name. With the
# old positional assignment this scrambling changes the DE table.
if (file.exists(de_path)) {
  de1 <- as.data.frame(read_tsv(de_path, show_col_types = FALSE))
  cd_perm <- coldata[rev(seq_len(nrow(coldata))), , drop = FALSE]
  d2 <- tempfile("dxp_"); dir.create(d2)
  write.csv(cd_perm, file.path(d2, "coldata.csv"), row.names = FALSE)
  st2 <- system2("Rscript", c(script, "--group", "G", "--counts", file.path(d, "counts.rds"),
                              "--coldata", file.path(d2, "coldata.csv"),
                              "--contrasts", file.path(d, "contrasts.tsv"),
                              "--design-formula", shQuote("~ condition"), "--fdr", "0.05", "--allow-no-offset", "true",
                              "--outdir", d2))
  p2 <- file.path(d2, "G_deseq2_results", "A_vs_B_DE.tsv")
  check(st2 == 0 && file.exists(p2), "permuted coldata still runs")
  if (file.exists(p2)) {
    de2 <- as.data.frame(read_tsv(p2, show_col_types = FALSE))
    m <- match(de1$gene_id, de2$gene_id)
    check(isTRUE(all.equal(de1$log2FoldChange, de2$log2FoldChange[m], tolerance = 1e-10)),
          "permuted coldata gives identical log2FoldChange")
    check(isTRUE(all.equal(de1$padj, de2$padj[m], tolerance = 1e-10)),
          "permuted coldata gives identical padj")
  }

  # A coldata missing a sample name must abort instead of aligning silently.
  cd_miss <- coldata[-1, , drop = FALSE]
  d3 <- tempfile("dxm_"); dir.create(d3)
  write.csv(cd_miss, file.path(d3, "coldata.csv"), row.names = FALSE)
  st3 <- system2("Rscript", c(script, "--group", "G", "--counts", file.path(d, "counts.rds"),
                              "--coldata", file.path(d3, "coldata.csv"),
                              "--contrasts", file.path(d, "contrasts.tsv"),
                              "--design-formula", shQuote("~ condition"), "--fdr", "0.05", "--allow-no-offset", "true",
                              "--outdir", d3), stdout = FALSE, stderr = FALSE)
  check(st3 != 0, "coldata missing a sample name fails the run")
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL DESEQ2 CONTRASTS TESTS PASSED\n")
