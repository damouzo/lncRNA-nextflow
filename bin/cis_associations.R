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
  writeLines(c("contrast\tn_pairs\tn_significant\tmedian_cor", "all\t0\t0\tNA"),
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
if (nrow(lnc_df) == 0 || nrow(pc_df) == 0) {
  empty_out("empty lncRNA or target set")
  quit(save = "no", status = 0)
}

# ---- 2. Genomic window search ----
lnc_gr <- GRanges(seqnames = lnc_df$chromosome,
                  ranges = IRanges(start = lnc_df$start, end = lnc_df$end),
                  strand = lnc_df$strand)
mcols(lnc_gr) <- DataFrame(gene_id = lnc_df$gene_id, catalog_id = lnc_df$catalog_id,
                           origin = lnc_df$origin, class_code = lnc_df$class_code)
pc_gr <- GRanges(seqnames = pc_df$chromosome,
                 ranges = IRanges(start = pc_df$start, end = pc_df$end),
                 strand = pc_df$strand)
mcols(pc_gr) <- DataFrame(gene_id = pc_df$gene_id, gene_name = pc_df$gene_name)

ov <- findOverlaps(lnc_gr + window, pc_gr, ignore.strand = TRUE)
pairs <- unique(data.frame(
  lncrna_id = lnc_df$gene_id[queryHits(ov)],
  lncrna_catalog_id = lnc_df$catalog_id[queryHits(ov)],
  lncrna_origin = lnc_df$origin[queryHits(ov)],
  lncrna_class_code = lnc_df$class_code[queryHits(ov)],
  target_gene_id = pc_df$gene_id[subjectHits(ov)],
  target_gene_name = pc_df$gene_name[subjectHits(ov)],
  stringsAsFactors = FALSE))
cat("Found", nrow(pairs), "lncRNA-gene pairs within +/-", window / 1000, "kb\n")
if (nrow(pairs) == 0) {
  empty_out("no cis pairs found")
  quit(save = "no", status = 0)
}

# ---- 3. Expression matrix (normalized) ----
counts <- readRDS(counts_rds)
mat <- log2(counts + 1)
norm_used <- "log2_counts_plus1"
if (normalize) {
  coldata <- tryCatch(read.csv(coldata_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  vsd <- tryCatch(DESeq2::vst(as.matrix(counts), blind = FALSE), error = function(e) NULL)
  if (!is.null(vsd)) {
    mat <- vsd
    norm_used <- "vst"
    if (!is.null(coldata) && "batch" %in% colnames(coldata)) {
      rownames(coldata) <- colnames(counts)
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

expressed_ids <- rownames(mat)
shared_lnc <- intersect(unique(pairs$lncrna_id), expressed_ids)
shared_gene <- intersect(unique(pairs$target_gene_id), expressed_ids)
cat("Expression data for", length(shared_lnc), "lncRNAs and", length(shared_gene), "target genes\n")

pairs <- pairs[pairs$lncrna_id %in% expressed_ids & pairs$target_gene_id %in% expressed_ids, , drop = FALSE]

cis_list <- list()
for (i in seq_len(nrow(pairs))) {
  lnc <- pairs$lncrna_id[i]
  gene <- pairs$target_gene_id[i]
  lnc_expr <- as.numeric(mat[lnc, ])
  gene_expr <- as.numeric(mat[gene, ])
  if (sd(lnc_expr) < 0.1 || sd(gene_expr) < 0.1) next
  cor_test <- tryCatch(cor.test(lnc_expr, gene_expr, method = "spearman"), error = function(e) NULL)
  if (is.null(cor_test)) next
  cis_list[[length(cis_list) + 1L]] <- data.frame(
    lncrna_id = lnc,
    lncrna_catalog_id = pairs$lncrna_catalog_id[i],
    lncrna_origin = pairs$lncrna_origin[i],
    lncrna_class_code = pairs$lncrna_class_code[i],
    target_gene_id = gene,
    target_gene_name = pairs$target_gene_name[i],
    spearman_rho = cor_test$estimate,
    spearman_pvalue = max(cor_test$p.value, .Machine$double.xmin),
    n_samples = length(lnc_expr),
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
writeLines(c("contrast\tn_pairs\tn_significant\tmedian_cor",
             sprintf("all\t%d\t%d\t%.4f", nrow(cis_results), nrow(cis_sig),
                     median(cis_results$spearman_rho, na.rm = TRUE))),
           out(paste0(group, "_cis_associations_summary.tsv")))
cat("DONE\n")
