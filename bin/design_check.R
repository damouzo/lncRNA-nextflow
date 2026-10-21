#!/usr/bin/env Rscript
# F3.2: per-contrast design rank check.
#
# For each contrast, on the samples that contrast actually uses (batch subset when
# the contrast pins a batch, otherwise the whole group), the design matrix must be
# full rank: qr(mm)$rank == ncol(mm). A rank-deficient design means a covariate is
# confounded with condition (e.g. every DDX41 sample in batch 1 and every Healthy
# sample in batch 2) and DESeq2 would silently produce unstable coefficients.
#
# Controlled by params.design_check: 'fail' (default) aborts, 'warn' records the
# status and continues. Always writes design_check.tsv for the record.
#
# usage:
#   design_check.R --group MSC --coldata coldata_MSC.csv --contrasts contrasts.tsv \
#       --design-formula '~ batch + condition' --mode fail --outdir .

suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

group <- get_arg("--group")
coldata_csv <- get_arg("--coldata")
contrasts_tsv <- get_arg("--contrasts")
formula <- get_arg("--design-formula", "~ condition")
mode <- get_arg("--mode", "fail")
outdir <- get_arg("--outdir", ".")
if (is.null(group) || is.null(coldata_csv) || is.null(contrasts_tsv)) {
  stop("usage: design_check.R --group G --coldata coldata.csv --contrasts contrasts.tsv")
}
if (!mode %in% c("fail", "warn")) stop("--mode must be fail|warn")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

coldata <- read.csv(coldata_csv, stringsAsFactors = FALSE)

# Same collapse rule as deseq2_contrasts.R: drop batch from the formula when it
# has a single level.
formula_used <- formula
if (grepl("batch", formula_used) && "batch" %in% colnames(coldata)) {
  batch_levels <- unique(coldata[["batch"]])
  if (length(batch_levels[nchar(batch_levels) > 0]) <= 1) formula_used <- "~ condition"
}

contrasts <- read_tsv(contrasts_tsv, show_col_types = FALSE,
                      col_names = c("contrast_name", "numerator", "denominator", "batch"))

rows <- list()
for (i in seq_len(nrow(contrasts))) {
  ct_name <- contrasts$contrast_name[i]
  num <- contrasts$numerator[i]
  den <- contrasts$denominator[i]
  ct_batch <- contrasts$batch[i]
  ct_batch <- if (is.na(ct_batch) || nchar(ct_batch) == 0) NULL else ct_batch

  if (!is.null(ct_batch) && "batch" %in% colnames(coldata)) {
    idx <- coldata[["batch"]] == ct_batch
    sub <- coldata[idx, , drop = FALSE]
    design <- "~ condition"
  } else {
    sub <- coldata
    design <- formula_used
  }

  conds <- unique(sub[["condition"]])
  applies <- num %in% conds && den %in% conds
  if (!applies) {
    rows[[length(rows) + 1L]] <- data.frame(
      group = group, contrast = ct_name, n_samples = nrow(sub),
      design = design, n_coef = NA_integer_, rank = NA_integer_,
      status = "not_applicable", stringsAsFactors = FALSE)
    next
  }

  if ("condition" %in% colnames(sub)) sub[["condition"]] <- factor(sub[["condition"]])
  mm <- model.matrix(as.formula(design), data = sub)
  rk <- qr(mm)$rank
  ncoef <- ncol(mm)
  ok <- rk == ncoef
  rows[[length(rows) + 1L]] <- data.frame(
    group = group, contrast = ct_name, n_samples = nrow(sub),
    design = design, n_coef = ncoef, rank = rk,
    status = if (ok) "ok" else "rank_deficient", stringsAsFactors = FALSE)

  if (!ok) {
    msg <- sprintf("design rank-deficient for %s %s: rank %d < %d coefficients (covariate confounded with condition)",
                   group, ct_name, rk, ncoef)
    if (mode == "fail") stop(msg) else cat("WARNING:", msg, "\n")
  }
}

out <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
write_tsv(out, file.path(outdir, "design_check.tsv"))
cat(sprintf("design_check [%s]: %d contrast(s), %d rank-deficient\n",
            group, nrow(out), sum(out$status == "rank_deficient")))
