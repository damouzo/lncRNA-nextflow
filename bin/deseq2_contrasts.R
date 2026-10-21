#!/usr/bin/env Rscript
# F3.3 / F3.5: DESeq2 per-contrast differential expression.
#
# Extracted from the DESEQ2_DE module so the logic is testable without Nextflow
# (bin/tests/test_deseq2_contrasts.R).
#
# F3.5: the LFC cutoff is NOT passed to results() as lfcThreshold. lfcThreshold
# changes the null hypothesis (|LFC| <= thr) and therefore pvalue/padj/stat, and
# the cutoff would be applied twice (here and in build_contrast_tables.py). The
# cutoff lives only in the `significant` flag, computed upstream on the SHUNKEN
# log2FoldChange (lfcShrink below). results() keeps alpha = fdr so independent
# filtering targets the configured FDR.
#
# F3.3: writes padj_na_summary.tsv — one row per contrast with the padj = NA
# count split by reason, plus the independent-filtering threshold.
#
# usage:
#   deseq2_contrasts.R --group MNC --counts counts.rds --coldata coldata.csv \
#       --contrasts contrasts.tsv --design-formula '~ condition' \
#       --fdr 0.05 --lfc 0.0

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

group <- get_arg("--group")
counts_rds <- get_arg("--counts")
coldata_csv <- get_arg("--coldata")
contrasts_tsv <- get_arg("--contrasts")
formula <- get_arg("--design-formula", "~ condition")
fdr <- as.numeric(get_arg("--fdr", "0.05"))
outdir <- get_arg("--outdir", ".")
if (is.null(group) || is.null(counts_rds) || is.null(coldata_csv) || is.null(contrasts_tsv)) {
  stop("usage: deseq2_contrasts.R --group G --counts counts.rds --coldata coldata.csv --contrasts contrasts.tsv")
}
res_dir <- file.path(outdir, paste0(group, "_deseq2_results"))
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

counts <- readRDS(counts_rds)
coldata <- read.csv(coldata_csv, stringsAsFactors = FALSE)
rownames(coldata) <- colnames(counts)

# Collapse ~ batch + condition to ~ condition when batch has a single level.
formula_used <- formula
if (grepl("batch", formula_used) && "batch" %in% colnames(coldata)) {
  batch_levels <- unique(coldata[["batch"]])
  batch_levels <- batch_levels[nchar(batch_levels) > 0]
  if (length(batch_levels) <= 1) {
    formula_used <- "~ condition"
    cat(sprintf("Batch has %d unique non-empty level(s) - using ~ condition\n",
                length(batch_levels)))
  }
}

# Global design diagnostic (rank of the whole-group design matrix).
diag_file <- file.path(res_dir, "design_diagnostics.txt")
sink(diag_file)
n_coef <- ncol(model.matrix(as.formula(formula_used), data = coldata))
rank_mm <- qr(model.matrix(as.formula(formula_used), data = coldata))$rank
cat(sprintf("formula: %s\n", formula_used))
cat(sprintf("rank: %d / %d coefficients\n", rank_mm, n_coef))
if (rank_mm < n_coef) {
  cat("WARNING: design matrix is rank-deficient.\n")
  cat("This usually means a covariate is confounded with condition or batch.\n")
  if ("batch" %in% colnames(coldata)) {
    cat("\nCondition x batch contingency table:\n")
    print(table(coldata$condition, coldata$batch))
  }
} else {
  cat("OK\n")
}
sink()
writeLines(readLines(diag_file))

contrasts <- read_tsv(contrasts_tsv, show_col_types = FALSE,
                      col_names = c("contrast_name", "numerator", "denominator", "batch"))

all_results <- list()
na_summary <- list()

for (i in seq_len(nrow(contrasts))) {
  ct_name <- contrasts$contrast_name[i]
  num <- contrasts$numerator[i]
  den <- contrasts$denominator[i]
  ct_batch <- contrasts$batch[i]
  ct_batch <- if (is.na(ct_batch) || nchar(ct_batch) == 0) NULL else ct_batch

  cat(sprintf("\n=== Contrast: %s (%s vs %s)", ct_name, num, den))
  if (!is.null(ct_batch)) cat(sprintf(" [batch: %s]", ct_batch))
  cat("\n")

  if (!is.null(ct_batch)) {
    if (!("batch" %in% colnames(coldata))) {
      cat("  SKIP: 'batch' column not found in coldata\n"); next
    }
    idx <- coldata[["batch"]] == ct_batch
    if (sum(idx) < 4) {
      cat(sprintf("  SKIP: only %d samples in batch '%s'\n", sum(idx), ct_batch)); next
    }
    sub_counts <- counts[, idx, drop = FALSE]
    sub_coldata <- coldata[idx, , drop = FALSE]
    design <- as.formula("~ condition")
  } else {
    sub_counts <- counts
    sub_coldata <- coldata
    design <- as.formula(formula_used)
  }

  sub_coldata$condition <- factor(sub_coldata$condition, levels = unique(sub_coldata$condition))
  if (!(num %in% levels(sub_coldata$condition))) {
    cat(sprintf("  SKIP: numerator '%s' not in condition levels\n", num)); next
  }
  if (!(den %in% levels(sub_coldata$condition))) {
    cat(sprintf("  SKIP: denominator '%s' not in condition levels\n", den)); next
  }
  sub_coldata$condition <- relevel(sub_coldata$condition, ref = den)

  # F3.5: no lfcThreshold here — alpha only.
  res <- tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_coldata, design = design)
    dds$condition <- droplevels(dds$condition)
    dds <- DESeq(dds)
    results(dds, contrast = c("condition", num, den), alpha = fdr)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message)); NULL
  })
  if (is.null(res)) next

  # lfcShrink(res = res) keeps res' pvalue/padj and gives the shrunken LFC that
  # build_contrast_tables.py reads via the log2FoldChange column.
  res_shrunk <- tryCatch({
    lfcShrink(dds, res = res, coef = resultsNames(dds)[grep(
      paste0("condition_", make.names(num)), resultsNames(dds))],
      type = "apeglm", quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  WARNING: apeglm failed (%s), using ashr\n", e$message))
    lfcShrink(dds, res = res, coef = resultsNames(dds)[grep(
      paste0("condition_", make.names(num)), resultsNames(dds))],
      type = "ashr", quiet = TRUE)
  })

  res_df <- as.data.frame(res_shrunk)
  res_df$contrast <- ct_name
  res_df$gene_id <- rownames(res_df)
  res_df$batch_subset <- if (is.null(ct_batch)) "all" else ct_batch

  fname <- file.path(res_dir, paste0(ct_name, "_DE.tsv"))
  write_tsv(res_df, fname)
  cat(sprintf("  Wrote %s (%d genes)\n", fname, nrow(res_df)))
  all_results[[ct_name]] <- res_df

  # F3.3: padj = NA diagnostics for this contrast.
  tested <- !is.na(res_df$pvalue)
  na_idx <- tested & is.na(res_df$padj)
  bm <- res_df$baseMean
  n_zero <- sum(na_idx & !is.na(bm) & bm == 0)
  ft <- metadata(res)$filterThreshold
  na_summary[[ct_name]] <- data.frame(
    contrast = ct_name,
    numerator = num, denominator = den,
    n_tested = sum(tested),
    n_padj_na = sum(na_idx),
    n_padj_na_baseMean_zero = n_zero,
    n_padj_na_independent_filtering = sum(na_idx) - n_zero,
    independent_filtering = isTRUE(metadata(res)$independentFiltering),
    filter_threshold = if (is.null(ft)) NA_real_ else as.numeric(ft),
    alpha = fdr,
    stringsAsFactors = FALSE
  )

  if (i == 1 || !exists("saved_dds")) saved_dds <- dds
}

if (length(all_results) > 0) {
  write_tsv(as.data.frame(do.call(rbind, all_results)),
            file.path(res_dir, "contrasts_summary.tsv"))
} else {
  write_tsv(data.frame(), file.path(res_dir, "contrasts_summary.tsv"))
}

if (length(na_summary) > 0) {
  write_tsv(as.data.frame(do.call(rbind, na_summary)),
            file.path(res_dir, "padj_na_summary.tsv"))
}

if (exists("saved_dds")) saveRDS(saved_dds, file.path(res_dir, "dds.rds"))

cat("DONE\n")
