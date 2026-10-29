#!/usr/bin/env Rscript
# F3.1: cis-regulatory candidate associations (gene-level, per norm_group).
#
# Extracted from the CIS_ASSOCIATIONS module so the logic is testable without
# Nextflow (bin/tests/test_cis_associations.R).
#
# Correlation is Spearman over a normalized expression matrix:
#   cis_normalize = true (default): VST (blind = FALSE), plus limma::removeBatchEffect
#     keeping condition when the group has >1 batch (MSCline/NB4). Falls back to
#     log2(counts + 1) if VST or limma are unavailable.
#   cis_normalize = false: log2(counts + 1) (previous behaviour).
# Outputs keep their format (contrast/n_pairs/n_significant/median_cor summary).
#
# usage:
#   cis_associations.R --group MNC --gene-catalog gene_catalog.tsv --counts counts.rds \
#       --coldata coldata_MNC.csv --cis-window 100000 --fdr 0.05 \
#       --normalize true --outdir .

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

group <- get_arg("--group")
gene_catalog_path <- get_arg("--gene-catalog")
counts_rds <- get_arg("--counts")
coldata_csv <- get_arg("--coldata")
window <- as.numeric(get_arg("--cis-window", "100000"))
fdr_cutoff <- as.numeric(get_arg("--fdr", "0.05"))
normalize <- tolower(get_arg("--normalize", "true")) %in% c("true", "t", "1", "yes")
min_samples <- as.integer(get_arg("--min-samples", "4"))
min_abs_rho <- as.numeric(get_arg("--min-abs-rho", "0"))
outdir <- get_arg("--outdir", ".")
if (is.null(group) || is.null(gene_catalog_path) || is.null(counts_rds)) {
  stop("usage: cis_associations.R --group G --gene-catalog catalog.tsv --counts counts.rds --coldata coldata.csv")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(outdir, f)

empty_out <- function(msg) {
  cat("WARNING:", msg, "\n")
  write_tsv(data.frame(), out(paste0(group, "_cis_pairs.tsv")))
  write_tsv(data.frame(), out(paste0(group, "_cis_pairs_significant.tsv")))
  writeLines(c("contrast\tn_pairs\tn_significant\tmedian_cor\tcorrelation_input\tmin_samples\tmin_abs_rho",
               "all\t0\t0\tNA\tNA\tNA\tNA"),
             out(paste0(group, "_cis_associations_summary.tsv")))
}

# ---- 1. Catalog split ----
catalog <- read_tsv(gene_catalog_path, show_col_types = FALSE)
lnc_df <- subset(catalog, origin %in% c("novel", "annotated_lncrna"))
pc_df <- subset(catalog, gene_biotype == "protein_coding")
if (nrow(pc_df) == 0) {
  pc_df <- subset(catalog, origin == "reference_other")
  cat("WARNING: no protein_coding biotype found, using all non-lncRNA genes\n")
}
cat("Catalog:", nrow(lnc_df), "lncRNAs (novel+annotated),", nrow(pc_df), "protein-coding targets\n")
lnc_df$strand <- ifelse(lnc_df$strand %in% c("+", "-", "*"), lnc_df$strand, "*")
pc_df$strand <- ifelse(pc_df$strand %in% c("+", "-", "*"), pc_df$strand, "*")
# lnc_position_class (antisense/intronic/overlapping/no_coding_overlap) describes
# every lncRNA, unlike class_code which gffcompare only assigns to novel ones.
# catalog_genes.tsv carries it; degrade to NA for a bare gene_catalog.tsv.
lnc_df$position_class <- if ("lnc_position_class" %in% colnames(lnc_df)) {
  as.character(lnc_df$lnc_position_class)
} else {
  NA_character_
}
if (nrow(lnc_df) == 0 || nrow(pc_df) == 0) {
  empty_out("empty lncRNA or target set")
  quit(save = "no", status = 0)
}

# ---- 2. Genomic window search ----
lnc_gr <- GRanges(seqnames = lnc_df$chromosome,
                  ranges = IRanges(start = lnc_df$start, end = lnc_df$end),
                  strand = lnc_df$strand)
mcols(lnc_gr) <- DataFrame(gene_id = lnc_df$gene_id, catalog_id = lnc_df$catalog_id,
                           origin = lnc_df$origin, class_code = lnc_df$class_code,
                           position_class = lnc_df$position_class)
pc_gr <- GRanges(seqnames = pc_df$chromosome,
                 ranges = IRanges(start = pc_df$start, end = pc_df$end),
                 strand = pc_df$strand)
mcols(pc_gr) <- DataFrame(gene_id = pc_df$gene_id, gene_name = pc_df$gene_name)

ov <- findOverlaps(lnc_gr + window, pc_gr, ignore.strand = TRUE)
qh <- queryHits(ov)
sh <- subjectHits(ov)
# 2.4: host/overlap marker — 0 bp gap means the spans overlap (putative host or
# overlapping locus), otherwise the pair is flanking within the window.
dist_bp <- as.integer(distance(lnc_gr[qh], pc_gr[sh], ignore.strand = TRUE))
pairs <- unique(data.frame(
  lncrna_id = lnc_df$gene_id[qh],
  lncrna_catalog_id = lnc_df$catalog_id[qh],
  lncrna_origin = lnc_df$origin[qh],
  lncrna_class_code = lnc_df$class_code[qh],
  lncrna_position_class = lnc_df$position_class[qh],
  target_gene_id = pc_df$gene_id[sh],
  target_gene_name = pc_df$gene_name[sh],
  overlaps_target = dist_bp == 0,
  distance_bp = dist_bp,
  stringsAsFactors = FALSE))
cat("Found", nrow(pairs), "lncRNA-gene pairs within +/-", window / 1000, "kb\n")
if (nrow(pairs) == 0) {
  empty_out("no cis pairs found")
  quit(save = "no", status = 0)
}

# ---- 3. Expression matrix (normalized) ----
counts <- readRDS(counts_rds)
coldata <- NULL
if (!is.null(coldata_csv) && file.exists(coldata_csv)) {
  coldata <- tryCatch(read.csv(coldata_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(coldata)) {
    # H0: align by sample name, never positionally (colnames(counts) is the
    # quant-dir listing order, not necessarily the coldata order).
    if (!("sample" %in% colnames(coldata))) {
      stop("coldata must have a 'sample' column to align it against the counts matrix")
    }
    coldata <- coldata[match(colnames(counts), coldata$sample), , drop = FALSE]
    stopifnot(!anyNA(coldata$sample),
              identical(as.character(coldata$sample), colnames(counts)))
    rownames(coldata) <- coldata$sample
  }
}

mat <- log2(counts + 1)
norm_used <- "log2_counts_plus1"
if (normalize) {
  vsd <- tryCatch(DESeq2::vst(as.matrix(counts), blind = FALSE), error = function(e) NULL)
  if (!is.null(vsd)) {
    mat <- vsd
    norm_used <- "vst"
    if (!is.null(coldata) && "batch" %in% colnames(coldata)) {
      batch <- coldata[["batch"]]
      if (length(unique(batch[nchar(as.character(batch)) > 0])) > 1 &&
          requireNamespace("limma", quietly = TRUE)) {
        design <- model.matrix(~ condition, data = coldata)
        mat <- limma::removeBatchEffect(mat, batch = batch, design = design)
        norm_used <- "vst_removeBatchEffect"
      }
    }
  } else {
    cat("WARNING: VST failed, falling back to log2(counts + 1)\n")
  }
}
cat("Normalization:", norm_used, "\n")

if (ncol(mat) < min_samples) {
  empty_out(sprintf("only %d samples (< min_samples %d)", ncol(mat), min_samples))
  quit(save = "no", status = 0)
}

# 2.4: correlate within-condition residuals so a shared condition effect cannot
# create a spurious association. Needs the per-sample condition; without it the
# raw normalized expression is used (previous behaviour).
mat_res <- mat
cor_input <- "expression"
n_conditions <- 0L
if (!is.null(coldata) && "condition" %in% colnames(coldata)) {
  condvec <- as.character(coldata$condition)
  mat_res <- mat
  for (lvl in unique(condvec)) {
    idx <- which(condvec == lvl)
    if (length(idx) >= 2) {
      mat_res[, idx] <- mat[, idx, drop = FALSE] - rowMeans(mat[, idx, drop = FALSE])
    } else {
      mat_res[, idx] <- 0
    }
  }
  cor_input <- "residuals_intra_condition"
  n_conditions <- length(unique(condvec))
}
cat("Correlation input:", cor_input, "\n")

# 4a.1: p-value for Spearman on the (possibly residualised) matrix.
# A t-approximation is used on the residual degrees of freedom instead of
# cor.test's exact AS 89 path: AS 89 assumes independent samples and, with the
# extreme rho values produced here, underflowed to .Machine$double.xmin, which a
# `max(p, .Machine$double.xmin)` floor then masked as a fake 2.2e-308 value.
# Removing one degree of freedom per condition mean removes the between-condition
# signal the residuals already discarded. The p-value is never floored: a
# non-finite or degenerate value becomes NA and the method is recorded, so
# "no p-value" is visible instead of a fabricated number.
spearman_pvalue <- function(rho, df) {
  if (is.na(rho) || !is.finite(rho)) return(list(p = NA_real_, method = "not_computed"))
  if (df < 1) return(list(p = NA_real_, method = "df_lt_1"))
  if (abs(rho) >= 1) return(list(p = NA_real_, method = "degenerate_rho"))
  tstat <- rho * sqrt(df / (1 - rho^2))
  p <- 2 * stats::pt(-abs(tstat), df = df)
  if (!is.finite(p) || p <= 0) return(list(p = NA_real_, method = "p_not_finite"))
  list(p = p, method = if (n_conditions > 0) "t_approx_residual_df" else "t_approx_df")
}

# 4a.1b: leave-one-sample-out stability. With n=10 residual samples a rank
# correlation can be driven by a single influential sample (or by rank ties in a
# low-expression vector), which produces the same rho across unrelated pairs.
# Nothing is dropped: the pair keeps its rho/p and gains these columns so a
# robust candidate can be recognised as one that survives removing any sample.
loo_stability <- function(x, y, full_rho) {
  n <- length(x)
  if (n < 4) return(list(min_abs = NA_real_, sign_stable = NA, max_delta = NA_real_))
  rhos <- numeric(0)
  for (i in seq_len(n)) {
    r <- suppressWarnings(stats::cor(x[-i], y[-i], method = "spearman"))
    if (is.finite(r)) rhos <- c(rhos, r)
  }
  if (length(rhos) == 0) return(list(min_abs = NA_real_, sign_stable = NA, max_delta = NA_real_))
  list(
    min_abs = min(abs(rhos)),
    sign_stable = all(rhos != 0) && all(sign(rhos) == sign(full_rho)),
    max_delta = max(abs(full_rho - rhos))
  )
}

expressed_ids <- rownames(mat)
shared_lnc <- intersect(unique(pairs$lncrna_id), expressed_ids)
shared_gene <- intersect(unique(pairs$target_gene_id), expressed_ids)
cat("Expression data for", length(shared_lnc), "lncRNAs and", length(shared_gene), "target genes\n")

pairs <- pairs[pairs$lncrna_id %in% expressed_ids & pairs$target_gene_id %in% expressed_ids, , drop = FALSE]

# 4a.1: residual degrees of freedom — n samples minus one per condition mean
# (residual design) minus the two correlation parameters.
df_effective <- ncol(mat_res) - n_conditions - 2L

cis_list <- list()
for (i in seq_len(nrow(pairs))) {
  lnc <- pairs$lncrna_id[i]
  gene <- pairs$target_gene_id[i]
  lnc_expr <- as.numeric(mat_res[lnc, ])
  gene_expr <- as.numeric(mat_res[gene, ])
  if (sd(lnc_expr) < 0.1 || sd(gene_expr) < 0.1) next
  rho <- tryCatch(suppressWarnings(stats::cor(lnc_expr, gene_expr, method = "spearman")),
                  error = function(e) NA_real_)
  if (is.na(rho) || !is.finite(rho)) next
  if (abs(rho) < min_abs_rho) next
  sp <- spearman_pvalue(rho, df_effective)
  # LOO is O(n) extra correlations per pair; only pairs that could still be
  # significant after BH (raw p < 0.05) need it, so the rest stay NA.
  loo <- if (is.finite(sp$p) && sp$p < 0.05) {
    loo_stability(lnc_expr, gene_expr, rho)
  } else {
    list(min_abs = NA_real_, sign_stable = NA, max_delta = NA_real_)
  }
  cis_list[[length(cis_list) + 1L]] <- data.frame(
    lncrna_id = lnc,
    lncrna_catalog_id = pairs$lncrna_catalog_id[i],
    lncrna_origin = pairs$lncrna_origin[i],
    lncrna_class_code = pairs$lncrna_class_code[i],
    lncrna_position_class = pairs$lncrna_position_class[i],
    target_gene_id = gene,
    target_gene_name = pairs$target_gene_name[i],
    overlaps_target = pairs$overlaps_target[i],
    distance_bp = pairs$distance_bp[i],
    spearman_rho = rho,
    spearman_pvalue = sp$p,
    n_samples = length(lnc_expr),
    df_effective = df_effective,
    p_method = sp$method,
    loo_min_abs_rho = loo$min_abs,
    loo_sign_stable = loo$sign_stable,
    loo_max_delta = loo$max_delta,
    stringsAsFactors = FALSE)
}
cis_results <- if (length(cis_list) > 0) do.call(rbind, cis_list) else data.frame()
if (nrow(cis_results) == 0) {
  empty_out("no valid correlation tests")
  quit(save = "no", status = 0)
}

# ---- 4. BH-FDR across the full test family ----
cis_results$padj <- p.adjust(cis_results$spearman_pvalue, method = "BH")
cis_results <- cis_results[order(cis_results$padj), ]
cis_sig <- subset(cis_results, padj < fdr_cutoff)

cat(sprintf("Cis associations: %d total, %d significant (FDR < %.2f)\n",
            nrow(cis_results), nrow(cis_sig), fdr_cutoff))
write_tsv(cis_results, out(paste0(group, "_cis_pairs.tsv")))
write_tsv(cis_sig, out(paste0(group, "_cis_pairs_significant.tsv")))
writeLines(c("contrast\tn_pairs\tn_significant\tmedian_cor\tcorrelation_input\tmin_samples\tmin_abs_rho",
             sprintf("all\t%d\t%d\t%.4f\t%s\t%d\t%g", nrow(cis_results), nrow(cis_sig),
                     median(cis_results$spearman_rho, na.rm = TRUE), cor_input,
                     min_samples, min_abs_rho)),
           out(paste0(group, "_cis_associations_summary.tsv")))
cat("DONE\n")
