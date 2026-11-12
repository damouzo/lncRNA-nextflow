#!/usr/bin/env Rscript
# F4.4 / 5f: tximport + condition-blind filter, now also emitting the per-gene
# expression tables and the pre-filter universe with a reason per gene.
#
# Extracted from the TXIMPORT_AND_FILTER module so the logic can be exercised
# without Nextflow (see bin/tests/test_tximport_filter.R). All previous outputs
# are kept byte-compatible: the filtered counts/TPM RDS, the PCA plot and the
# filter stats. New:
#   <group>_counts.tsv    filtered counts (gene_id + annotation + samples)
#   <group>_tpm.tsv       filtered TPM (same layout)
#   <group>_universe.tsv  every imported gene with keep/filter_status/filter_reason
#
# With --catalog (catalog_genes.tsv) the TSVs gain annotation columns right after
# gene_id, using the catalog's own names so every gen-centric table reads the same
# way: `gene_group` (protein_coding / annotated_lncRNA / novel_lncRNA / other),
# `catalog_id` (stable id; NA for reference genes), `origin`, `class_code`,
# `tier` and `lnc_position_class`. The numeric RDS matrices stay untouched so
# DESeq2/cis keep consuming a clean matrix.
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

# ---- D4 invariants: the quantified transcripts must match tx2gene exactly ----
# tximport only `stop()`s when NOTHING matches; a partial mismatch (e.g. versioned
# quant names vs stripped tx2gene ids, which the `--gencode` index flag does not
# strip) is a message and a silent drop — the reference genes then vanish and the
# DE universe is quietly wrong. Assert both directions here, before any filtering.
if (!identical(txi$countsFromAbundance, "no")) {
  stop(sprintf("expected countsFromAbundance='no', got '%s'", txi$countsFromAbundance))
}
quant_tx <- character(0)
num_reads_total <- 0
for (f in quant_files) {
  q <- read_tsv(f, show_col_types = FALSE,
                col_types = cols(Name = col_character(), .default = col_double()))
  quant_tx <- c(quant_tx, q$Name)
  num_reads_total <- num_reads_total + sum(q$NumReads, na.rm = TRUE)
}
quant_tx <- unique(quant_tx)
unmatched <- setdiff(quant_tx, tx2gene$TXNAME)
if (length(unmatched) > 0) {
  stop(sprintf("%d of %d quantified transcripts are absent from tx2gene (e.g. %s). The index FASTA ids and tx2gene ids diverged (version suffixes?); check that BUILD_TX2GENE built tx2gene from the same universe as the index.",
               length(unmatched), length(quant_tx), paste(head(unmatched, 5), collapse = ", ")))
}
if (length(tx2gene$TXNAME) != length(unique(tx2gene$TXNAME))) {
  stop("tx2gene has duplicate transcript ids; tximport keeps only the first match, silently misassigning reads")
}
# Every gene count must equal the sum of its transcripts' NumReads (~sum over the
# matched universe). Float counts need a tolerance, not an equality test.
read_sum <- sum(txi$counts)
if (abs(read_sum - num_reads_total) > 1e-3 * max(1, num_reads_total)) {
  stop(sprintf("sum(gene counts)=%.3f != sum(NumReads)=%.3f: tx2gene aggregation lost or duplicated reads",
               read_sum, num_reads_total))
}

# Condition-blind filter: >= min_count counts in >= max(3, floor(n/4)) samples
min_n <- max(3, floor(ncol(txi$counts) / 4))
all_counts <- txi$counts
n_detected <- rowSums(all_counts >= min_count)
keep <- n_detected >= min_n

dds <- DESeqDataSetFromTximport(txi, coldata, ~1)
dds <- dds[keep, ]

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
# H/2.2: colour the PCA by the biological metadata, never by `sample`. `seq_batch`
# (sequencing configuration: flowcell/lane) is drawn as shape so a technical
# structure is visible without hiding the condition split; `batch` is the fallback
# shape for cell-line experiments. A column that is empty/NA for the whole group
# (single level) is dropped, and with no shape column left the plot is colour-only.
has_levels <- function(x) {
  x <- as.character(x)
  length(unique(x[!is.na(x) & nzchar(x)])) > 1
}
color_var <- if ("condition" %in% colnames(coldata)) "condition" else colnames(coldata)[1]
shape_var <- NULL
for (cand in c("seq_batch", "batch")) {
  if (cand %in% colnames(coldata) && has_levels(coldata[[cand]])) {
    shape_var <- cand
    break
  }
}

# plotPCA maps a single aesthetic (several intgroup collapse into one colour), so
# take the projected coordinates and build the plot here to keep colour and shape
# separate.
pca_data <- plotPCA(vsd, intgroup = unique(c(color_var, shape_var)), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$.color_group <- factor(pca_data[[color_var]])
mapping <- ggplot2::aes(x = PC1, y = PC2, color = .color_group)
if (!is.null(shape_var)) {
  pca_data$.shape_group <- factor(pca_data[[shape_var]])
  mapping <- ggplot2::aes(x = PC1, y = PC2, color = .color_group, shape = .shape_group)
}
pca <- ggplot2::ggplot(pca_data, mapping) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = name), size = 2.5, vjust = -0.8,
                     check_overlap = TRUE) +
  ggplot2::labs(
    title    = sprintf("%s — variance-stabilized PCA (VST)", group),
    subtitle = sprintf("%d samples; colour = %s%s", ncol(vsd), color_var,
                       if (is.null(shape_var)) "" else paste0("; shape = ", shape_var)),
    color    = tools::toTitleCase(color_var),
    x        = sprintf("PC1: %d%% variance", percentVar[1]),
    y        = sprintf("PC2: %d%% variance", percentVar[2])
  )
if (!is.null(shape_var)) {
  pca <- pca + ggplot2::labs(shape = tools::toTitleCase(shape_var))
}
ggsave(out(paste0(group, "_sample_qc_pca.png")), pca, width = 10, height = 8)
ggsave(out(paste0(group, "_sample_qc_pca.pdf")), pca, width = 10, height = 8)

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

# Annotation for the expression tables: gene_id -> gene_group / catalog_id / origin /
# class_code / tier / lnc_position_class, using the catalog's own column names.
# catalog_genes.tsv already carries them; if gene_group is missing, derive it from
# origin + biotype with the same 4-value rule used across the deliverables
# (protein_coding / annotated_lncRNA / novel_lncRNA / other).
ANNOTATION_COLS <- c("gene_group", "catalog_id", "origin", "class_code", "tier",
                     "lnc_position_class")
ann_maps <- NULL
if (!is.null(catalog_path) && file.exists(catalog_path)) {
  cat_df <- tryCatch(read_tsv(catalog_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(cat_df) || !("gene_id" %in% colnames(cat_df))) {
    cat("WARNING: catalog unreadable or missing gene_id; expression tables stay unannotated\n")
  } else {
    gid <- cat_df$gene_id
    get_col <- function(name) {
      if (name %in% colnames(cat_df)) as.character(cat_df[[name]]) else rep(NA_character_, length(gid))
    }
    gg <- if ("gene_group" %in% colnames(cat_df)) {
      get_col("gene_group")
    } else {
      origin  <- get_col("origin")
      biotype <- get_col("gene_biotype")
      ifelse(origin == "novel", "novel_lncRNA",
             ifelse(origin == "annotated_lncrna", "annotated_lncRNA",
                    ifelse(!is.na(biotype) & biotype == "protein_coding", "protein_coding", "other")))
    }
    ann_vals <- list(gene_group = gg, catalog_id = get_col("catalog_id"),
                     origin = get_col("origin"), class_code = get_col("class_code"),
                     tier = get_col("tier"), lnc_position_class = get_col("lnc_position_class"))
    ann_maps <- setNames(lapply(ann_vals[ANNOTATION_COLS], function(v) setNames(v, gid)),
                         ANNOTATION_COLS)
  }
}

# Insert the annotation columns right after gene_id; samples keep their order.
with_annotation <- function(df) {
  if (is.null(ann_maps)) return(df)
  gid <- df[["gene_id"]]
  ann <- as.data.frame(lapply(ann_maps, function(m) unname(m[gid])),
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

# ---- 4b.6: normalized + VST matrices (same annotation layout as counts/TPM) ----
# counts(normalized = TRUE) are DESeq2's median-of-ratios values; the VST matrix
# is what the PCA/correlation/heatmaps use. Both mirror counts.tsv. Size factors
# must exist before the normalized counts are read (vsd did not store them in dds).
if (is.null(DESeq2::sizeFactors(dds))) dds <- DESeq2::estimateSizeFactors(dds)
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
norm_df <- with_annotation(data.frame(gene_id = rownames(norm_counts), norm_counts,
                                      check.names = FALSE))
write_tsv(norm_df, out(paste0(group, "_normalized_counts.tsv")))

vst_mat <- SummarizedExperiment::assay(vsd)
vst_df <- with_annotation(data.frame(gene_id = rownames(vst_mat), vst_mat,
                                     check.names = FALSE))
write_tsv(vst_df, out(paste0(group, "_vst.tsv")))

# ---- 4b.7: enriched sample metadata (size factors + per-sample detection) ----
if (is.null(DESeq2::sizeFactors(dds))) dds <- DESeq2::estimateSizeFactors(dds)
size_factors <- DESeq2::sizeFactors(dds)
samples <- colnames(keep_counts)
# lncRNA gene ids known to the annotation (novel + annotated), for the detection split
lncrna_ids <- if (!is.null(ann_maps$gene_group)) {
  names(ann_maps$gene_group)[ann_maps$gene_group %in% c("novel_lncRNA", "annotated_lncRNA")]
} else character(0)
lncrna_kept <- intersect(rownames(keep_counts), lncrna_ids)
meta_col <- function(name) {
  if (name %in% colnames(aligned_coldata)) as.character(aligned_coldata[[name]]) else NA_character_
}
sample_meta <- data.frame(
  sample            = samples,
  norm_group        = meta_col("norm_group"),
  condition         = meta_col("condition"),
  batch             = meta_col("batch"),
  seq_batch         = meta_col("seq_batch"),
  size_factor       = as.numeric(size_factors[samples]),
  total_reads       = as.numeric(colSums(txi$counts)[samples]),
  n_genes_detected  = as.integer(colSums(keep_counts > 0)),
  n_lncrna_detected = as.integer(colSums(keep_counts[lncrna_kept, , drop = FALSE] > 0)),
  stringsAsFactors  = FALSE, check.names = FALSE
)
write_tsv(sample_meta, out(paste0(group, "_sample_metadata.tsv")))

# ---- 4b.8: lncRNA type composition per sample (annotated + novel pooled) ----
# Type from lnc_position_class mapped to friendly names; detection at TPM >= 1.
TYPE_MAP <- c(no_coding_overlap = "lincRNA", antisense = "antisense",
              intronic = "intronic", overlapping = "overlapping")
comp_cols <- c("sample", "condition", "lncRNA_type", "n_detected",
               "fraction_detected", "tpm_sum", "fraction_tpm")
comp <- if (length(lncrna_kept) > 0) {
  tpm_lnc <- tpm[lncrna_kept, , drop = FALSE]
  pos <- ann_maps$lnc_position_class[lncrna_kept]
  ltype <- ifelse(is.na(pos) | pos %in% c("", "NA"), "unclassified",
                  ifelse(pos %in% names(TYPE_MAP), unname(TYPE_MAP[pos]), "other"))
  comp_rows <- list()
  for (s in samples) {
    v <- as.numeric(tpm_lnc[, s])
    det <- v >= 1
    total_tpm <- sum(v, na.rm = TRUE)
    for (ty in sort(unique(ltype))) {
      idx <- ltype == ty
      comp_rows[[length(comp_rows) + 1L]] <- data.frame(
        sample = s,
        condition = sample_meta$condition[match(s, sample_meta$sample)],
        lncRNA_type = ty,
        n_detected = sum(det[idx], na.rm = TRUE),
        fraction_detected = sum(det[idx], na.rm = TRUE) / length(v),
        tpm_sum = sum(v[idx], na.rm = TRUE),
        fraction_tpm = if (total_tpm > 0) sum(v[idx], na.rm = TRUE) / total_tpm else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, comp_rows)
} else {
  as.data.frame(setNames(
    lapply(comp_cols, function(x) if (x == "n_detected") integer(0) else
      if (x %in% c("fraction_detected", "tpm_sum", "fraction_tpm")) double(0) else character(0)),
    comp_cols))
}
write_tsv(comp, out(paste0(group, "_lncrna_type_composition.tsv")))

# ---- QC figures (dual format: PDF for print, PNG for slides) ----
# Non-blocking: a plotting hiccup logs a warning and never fails the task.
tryCatch({
  save_panel <- function(name, plot, width, height) {
    ggplot2::ggsave(out(paste0(group, "_", name, ".png")), plot,
                    width = width, height = height, dpi = 150)
    ggplot2::ggsave(out(paste0(group, "_", name, ".pdf")), plot,
                    width = width, height = height)
  }
  lev <- sample_meta$sample[order(sample_meta$condition, sample_meta$sample)]

  # Sample-to-sample correlation (VST), ordered by condition then sample.
  cor_mat <- stats::cor(vst_mat)
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("sample_x", "sample_y", "cor")
  cor_df$sample_x <- factor(cor_df$sample_x, levels = lev)
  cor_df$sample_y <- factor(cor_df$sample_y, levels = rev(lev))
  p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(sample_x, sample_y, fill = cor)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#1565c0", mid = "white", high = "#b71c1c",
                                  midpoint = stats::median(cor_mat)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(title = sprintf("%s - sample correlation (VST)", group),
                  x = NULL, y = NULL, fill = "rho")
  save_panel("sample_correlation", p_cor, 8, 7)

  # Features detected per sample (genes vs lncRNAs).
  det_df <- data.frame(
    sample = rep(sample_meta$sample, 2),
    feature = rep(c("genes", "lncRNAs"), each = nrow(sample_meta)),
    n = c(sample_meta$n_genes_detected, sample_meta$n_lncrna_detected),
    stringsAsFactors = FALSE)
  det_df$sample <- factor(det_df$sample, levels = lev)
  p_det <- ggplot2::ggplot(det_df, ggplot2::aes(sample, n, fill = feature)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(title = sprintf("%s - features detected per sample", group),
                  x = NULL, y = "n detected (count > 0)", fill = NULL)
  save_panel("detected_features", p_det, 9, 5)

  # TPM distribution per sample.
  tpm_long <- do.call(rbind, lapply(samples, function(s) data.frame(
    sample = s, log2tpm = log2(as.numeric(tpm[, s]) + 1), stringsAsFactors = FALSE)))
  tpm_long$sample <- factor(tpm_long$sample, levels = lev)
  p_tpm <- ggplot2::ggplot(tpm_long, ggplot2::aes(sample, log2tpm)) +
    ggplot2::geom_violin(fill = "#4c72b0", alpha = 0.6, scale = "width") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(title = sprintf("%s - TPM distribution per sample", group),
                  x = NULL, y = "log2(TPM + 1)")
  save_panel("tpm_distribution", p_tpm, 9, 5)

  # lncRNA type composition (fraction of lncRNA TPM per type), samples by condition.
  if (nrow(comp) > 0 && "fraction_tpm" %in% colnames(comp)) {
    p_type <- ggplot2::ggplot(comp, ggplot2::aes(lncRNA_type, fraction_tpm, fill = condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.15),
                          size = 1.2, alpha = 0.7) +
      ggplot2::theme_minimal() +
      labs(title = sprintf("%s - lncRNA type composition", group),
           x = "lncRNA type", y = "fraction of lncRNA TPM", fill = "condition")
    save_panel("lncrna_type_composition", p_type, 8, 5)
  }
}, error = function(e) cat("WARNING: QC figures skipped:", conditionMessage(e), "\n"))

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
