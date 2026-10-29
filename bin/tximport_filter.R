#!/usr/bin/env Rscript
# F4.4 / 5f: tximport + condition-blind filter, now also emitting the per-gene
# expression tables and the pre-filter universe with a reason per gene.
#
# Extracted from the TXIMPORT_AND_FILTER module so the logic can be exercised
# without Nextflow (see bin/tests/test_tximport_filter.R). All previous outputs
# are kept byte-compatible: the filtered counts/TPM RDS, the PCA plot and the
# filter stats. New:
#   <group>_counts.tsv    filtered counts (gene_id + gene_group/catalog_id + samples)
#   <group>_tpm.tsv       filtered TPM (same layout)
#   <group>_universe.tsv  every imported gene with keep/filter_status/filter_reason
#
# With --catalog (catalog_genes.tsv) the TSVs gain two annotation columns right
# after gene_id: `gene_group` (protein_coding / annotated_lncRNA / novel_lncRNA /
# other) and `catalog_id` (stable id; NA for reference genes). The numeric RDS
# matrices stay untouched so DESeq2/cis keep consuming a clean matrix.
#
# usage:
#   tximport_filter.R --group MNC --coldata coldata_MNC.csv --tx2gene tx2gene.tsv
#                     [--catalog catalog_genes.tsv] [--quant-pattern "_quant$"] [--min-count 10]

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
catalog_path <- get_arg("--catalog")
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
# H0/2.3: a coldata sample without a quantification directory is a hard error.
# Silently intersecting would drop it and shift the sample set without warning.
if (length(shared) != nrow(coldata)) {
  missing <- setdiff(rownames(coldata), names(quant_files))
  cat("Quant sample names:", names(quant_files), "\n")
  cat("Coldata sample names:", rownames(coldata), "\n")
  stop(sprintf("n_shared (%d) != nrow(coldata) (%d); missing quant for: %s",
               length(shared), nrow(coldata), paste(missing, collapse = ", ")))
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
# H/2.2: colour the PCA by the biological metadata, never by `sample`. `batch`
# joins `condition` only when the samplesheet actually has it. plotPCA pastes
# every intgroup value into the legend label, so a batch column that is empty/NA
# for the whole group (single-batch groups) would render as "<condition>:NA" —
# drop it unless it actually separates samples.
intgroup <- intersect(c("condition", "batch"), colnames(coldata))
has_levels <- function(x) {
  x <- as.character(x)
  length(unique(x[!is.na(x) & nzchar(x)])) > 1
}
if ("batch" %in% intgroup && !has_levels(coldata[["batch"]])) {
  intgroup <- setdiff(intgroup, "batch")
}
if (length(intgroup) == 0) intgroup <- colnames(coldata)[1]
pca <- plotPCA(vsd, intgroup = intgroup) +
  ggplot2::labs(
    title    = sprintf("%s — variance-stabilized PCA (VST)", group),
    subtitle = sprintf("%d samples; colour = %s", ncol(vsd), paste(intgroup, collapse = " / ")),
    colour   = paste(tools::toTitleCase(intgroup), collapse = " / ")
  )
ggsave(out(paste0(group, "_sample_qc_pca.png")), pca, width = 10, height = 8)

keep_counts <- counts(dds)
saveRDS(keep_counts, out(paste0(group, "_tximport_counts.rds")))
saveRDS(txi$abundance, out(paste0(group, "_tximport_tpm.rds")))
# 2.1: persist the tximport-derived DESeqDataSet (filtered to the kept genes) so
# DESEQ2_DE can run on it and keep the transcript-length offset (avgTxLength).
saveRDS(dds, out(paste0(group, "_tximport_dds.rds")))

# H0: publish the coldata aligned to the counts columns. Quant dirs are listed
# alphabetically, so colnames(keep_counts) is NOT the samplesheet order; consumers
# (DESeq2/cis/heterogeneity) must read this file instead of the samplesheet-order
# coldata so a label can never be paired positionally with the wrong sample.
# The name must NOT equal the input coldata_<group>.csv: Nextflow stages that input
# as a symlink in this task dir, so writing the same basename would follow it back
# into the upstream work dir and corrupt the cached SPLIT_BY_NORMGROUP output.
aligned_coldata <- coldata[colnames(keep_counts), , drop = FALSE]
stopifnot(!anyNA(aligned_coldata[["sample"]]),
          identical(as.character(aligned_coldata[["sample"]]), colnames(keep_counts)))
write.csv(aligned_coldata, out(paste0("coldata_aligned_", group, ".csv")),
          row.names = FALSE, quote = FALSE)

stats <- data.frame(
  total_genes = nrow(txi$counts),
  kept_genes  = sum(keep),
  removed     = nrow(txi$counts) - sum(keep),
  n_samples   = ncol(dds)
)
write_tsv(stats, out(paste0(group, "_filter_stats.tsv")))

# Annotation for the expression tables: gene_id -> gene_group / catalog_id.
# catalog_genes.tsv already carries gene_group; if it is missing, derive it from
# origin + biotype with the same 4-value rule used across the deliverables
# (protein_coding / annotated_lncRNA / novel_lncRNA / other).
gene_group_map <- NULL
catalog_id_map <- NULL
if (!is.null(catalog_path) && file.exists(catalog_path)) {
  cat_df <- tryCatch(read_tsv(catalog_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(cat_df) || !("gene_id" %in% colnames(cat_df))) {
    cat("WARNING: catalog unreadable or missing gene_id; expression tables stay unannotated\n")
  } else {
    if ("gene_group" %in% colnames(cat_df)) {
      gg <- as.character(cat_df$gene_group)
    } else {
      origin  <- cat_df$origin
      biotype <- if ("gene_biotype" %in% colnames(cat_df)) cat_df$gene_biotype else NA_character_
      gg <- ifelse(origin == "novel", "novel_lncRNA",
             ifelse(origin == "annotated_lncrna", "annotated_lncRNA",
              ifelse(!is.na(biotype) & biotype == "protein_coding", "protein_coding", "other")))
    }
    cid <- if ("catalog_id" %in% colnames(cat_df)) as.character(cat_df$catalog_id) else NA_character_
    gene_group_map <- setNames(gg, cat_df$gene_id)
    catalog_id_map <- setNames(cid, cat_df$gene_id)
  }
}

# Insert the annotation columns right after gene_id; samples keep their order.
with_annotation <- function(df) {
  if (is.null(gene_group_map)) return(df)
  gid <- df[["gene_id"]]
  ann <- data.frame(gene_group = unname(gene_group_map[gid]),
                    catalog_id = unname(catalog_id_map[gid]),
                    stringsAsFactors = FALSE, check.names = FALSE)
  cbind(df["gene_id"], ann, df[setdiff(colnames(df), "gene_id")])
}

# TSV expression tables (filtered), gene_id (+ annotation) + one column per sample
counts_df <- with_annotation(data.frame(gene_id = rownames(keep_counts), keep_counts,
                                        check.names = FALSE))
write_tsv(counts_df, out(paste0(group, "_counts.tsv")))

tpm <- txi$abundance[rownames(keep_counts), , drop = FALSE]
tpm_df <- with_annotation(data.frame(gene_id = rownames(tpm), tpm, check.names = FALSE))
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
