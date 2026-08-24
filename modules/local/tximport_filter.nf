// tximport + independent filtering (condition-blind)
// Imports Salmon quantifications, builds colData from samplesheet CSV,
// runs PCA sample QC, applies condition-blind low-count filter

process TXIMPORT_AND_FILTER {

    input:
    tuple val(group), path(quant_dirs), path(coldata_csv)  // keyed: no positional pairing downstream
    path tx2gene_detailed

    output:
    tuple val(group), path("${group}_tximport_counts.rds"), emit: counts_rds
    // Keyed by group (unlike counts_rds's sibling outputs below): needed so downstream
    // consumers can .join() it per norm_group instead of risking the F1-style hazard of
    // pairing single-execution/queue channels positionally across >1 norm_group.
    tuple val(group), path("${group}_tximport_tpm.rds"),    emit: tpm_rds
    path "${group}_sample_qc_pca.png",                      emit: pca_plot
    path "${group}_filter_stats.tsv",                       emit: filter_stats

    script:
    """
    #!/usr/bin/env Rscript
    library(tximport)
    library(DESeq2)
    library(readr)
    library(ggplot2)

    # Read coldata from the samplesheet CSV
    coldata <- read.csv("${coldata_csv}", stringsAsFactors = FALSE)
    rownames(coldata) <- coldata\$sample

    # Read tx2gene mapping (first two columns, drop category)
    tx2gene <- read_tsv("${tx2gene_detailed}", show_col_types = FALSE)
    tx2gene <- tx2gene[, 1:2]
    colnames(tx2gene) <- c("TXNAME", "GENEID")

    # Find all quant.sf files inside staged {sample}_quant/ directories
    quant_dirs <- list.files(pattern = "_quant\$", full.names = TRUE)
    if (length(quant_dirs) == 0) {
        stop("No _quant directories found in working directory")
    }
    quant_files <- file.path(quant_dirs, "quant.sf")
    names(quant_files) <- gsub("_quant\$", "", basename(quant_dirs))

    cat("Found", length(quant_files), "quantification files\\n")

    # Match samples between quant files and coldata
    shared <- intersect(names(quant_files), rownames(coldata))
    if (length(shared) == 0) {
        cat("Quant sample names:", names(quant_files), "\\n")
        cat("Coldata sample names:", rownames(coldata), "\\n")
        stop("No matching samples between quant files and coldata")
    }
    quant_files <- quant_files[shared]
    coldata <- coldata[shared, , drop = FALSE]

    cat("Matched", length(shared), "samples\\n")

    txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                    ignoreTxVersion = TRUE)

    # Condition-blind filtering: keep genes with >= 10 counts in >= 3 samples or 1/4 of samples
    min_n <- max(3, floor(ncol(txi\$counts) / 4))
    dds <- DESeqDataSetFromTximport(txi, coldata, ~1)
    keep <- rowSums(counts(dds) >= 10) >= min_n
    dds <- dds[keep, ]

    # PCA for sample-level QC
    # varianceStabilizingTransformation fits the trend on all genes instead of
    # sampling nsub=1000 like vst(), so it works on small lncRNA catalogs
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    n_intgroup <- min(3, ncol(coldata))
    pca <- plotPCA(vsd, intgroup = colnames(coldata)[1:n_intgroup])
    ggsave("${group}_sample_qc_pca.png", pca, width = 10, height = 8)

    saveRDS(counts(dds), "${group}_tximport_counts.rds")
    saveRDS(txi\$abundance, "${group}_tximport_tpm.rds")

    stats <- data.frame(
        total_genes  = nrow(txi\$counts),
        kept_genes   = sum(keep),
        removed      = nrow(txi\$counts) - sum(keep),
        n_samples    = ncol(dds)
    )
    write_tsv(stats, "${group}_filter_stats.tsv")
    cat("Filter stats:", nrow(stats), "rows written\\n")
    """
}
