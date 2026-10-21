#!/usr/bin/env Rscript
# F4.4 / 5f: tximport + condition-blind filter, now also emitting the per-gene
# expression tables and the pre-filter universe with a reason per gene.
#
# Extracted from the TXIMPORT_AND_FILTER module so the logic can be exercised
# without Nextflow (see bin/tests/test_tximport_filter.R). All previous outputs
# are kept byte-compatible: the filtered counts/TPM RDS, the PCA plot and the
# filter stats. New:
#   <group>_counts.tsv    filtered counts (gene_id + one column per sample)
#   <group>_tpm.tsv       filtered TPM
#   <group>_universe.tsv  every imported gene with keep/filter_status/filter_reason
#
# usage:
#   tximport_filter.R --group MNC --coldata coldata_MNC.csv --tx2gene tx2gene.tsv
#                     [--quant-pattern "_quant$"] [--min-count 10]

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

group <- get_arg("--group")
coldata_path <- get_arg("--coldata")
tx2gene_path <- get_arg("--tx2gene")
quant_pattern <- get_arg("--quant-pattern", "_quant$")
min_count <- as.integer(get_arg("--min-count", "10"))
outdir <- get_arg("--outdir", ".")
if (is.null(group) || is.null(coldata_path) || is.null(tx2gene_path)) {
  stop("usage: tximport_filter.R --group G --coldata coldata.csv --tx2gene tx2gene.tsv")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(outdir, f)

coldata <- read.csv(coldata_path, stringsAsFactors = FALSE)
rownames(coldata) <- coldata[["sample"]]

tx2gene <- read_tsv(tx2gene_path, show_col_types = FALSE)[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")

quant_dirs <- list.files(pattern = quant_pattern, full.names = TRUE)
if (length(quant_dirs) == 0) stop("No _quant directories found in working directory")
quant_files <- file.path(quant_dirs, "quant.sf")
names(quant_files) <- gsub(quant_pattern, "", basename(quant_dirs))
cat("Found", length(quant_files), "quantification files\n")

shared <- intersect(names(quant_files), rownames(coldata))
if (length(shared) == 0) {
  cat("Quant sample names:", names(quant_files), "\n")
  cat("Coldata sample names:", rownames(coldata), "\n")
  stop("No matching samples between quant files and coldata")
}
quant_files <- quant_files[shared]
coldata <- coldata[shared, , drop = FALSE]
cat("Matched", length(shared), "samples\n")

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE)

# Condition-blind filter: >= min_count counts in >= max(3, floor(n/4)) samples
min_n <- max(3, floor(ncol(txi$counts) / 4))
all_counts <- txi$counts
n_detected <- rowSums(all_counts >= min_count)
keep <- n_detected >= min_n

dds <- DESeqDataSetFromTximport(txi, coldata, ~1)
dds <- dds[keep, ]

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
n_intgroup <- min(3, ncol(coldata))
pca <- plotPCA(vsd, intgroup = colnames(coldata)[1:n_intgroup])
ggsave(out(paste0(group, "_sample_qc_pca.png")), pca, width = 10, height = 8)

keep_counts <- counts(dds)
saveRDS(keep_counts, out(paste0(group, "_tximport_counts.rds")))
saveRDS(txi$abundance, out(paste0(group, "_tximport_tpm.rds")))

stats <- data.frame(
  total_genes = nrow(txi$counts),
  kept_genes  = sum(keep),
  removed     = nrow(txi$counts) - sum(keep),
  n_samples   = ncol(dds)
)
write_tsv(stats, out(paste0(group, "_filter_stats.tsv")))

# TSV expression tables (filtered), gene_id + one column per sample
counts_df <- data.frame(gene_id = rownames(keep_counts), keep_counts,
                        check.names = FALSE)
write_tsv(counts_df, out(paste0(group, "_counts.tsv")))

tpm <- txi$abundance[rownames(keep_counts), , drop = FALSE]
tpm_df <- data.frame(gene_id = rownames(tpm), tpm, check.names = FALSE)
write_tsv(tpm_df, out(paste0(group, "_tpm.tsv")))

# Pre-filter universe with a reason per gene (nothing is lost)
universe <- data.frame(
  gene_id        = rownames(all_counts),
  total_count    = rowSums(all_counts),
  n_samples_ge_min = n_detected,
  keep           = keep,
  filter_status  = ifelse(keep, "kept", "filtered"),
  filter_reason  = ifelse(keep, "", "low_expression"),
  stringsAsFactors = FALSE
)
write_tsv(universe, out(paste0(group, "_universe.tsv")))

cat("Filter stats:", stats$total_genes, "genes,", stats$kept_genes, "kept\n")
