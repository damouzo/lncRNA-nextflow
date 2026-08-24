// Cis-regulatory candidate associations
// GenomicRanges-based window search for lncRNA–protein_coding gene pairs,
// Spearman correlation (global + within-condition), BH-FDR correction
// across the full test family per norm_group.
// lncRNA side: novel + annotated lncRNAs from the gene catalog.
// Targets: protein-coding genes from the same catalog.

process CIS_ASSOCIATIONS {

    input:
    tuple val(group), path(counts_rds), path(deseq2_summary)  // keyed: no positional pairing downstream
    path gene_catalog       // gene_catalog.tsv: origin + coordinates for every gene

    output:
    tuple val(group), path("${group}_cis_pairs.tsv"),              emit: cis_pairs
    tuple val(group), path("${group}_cis_pairs_significant.tsv"),  emit: cis_pairs_sig
    path "${group}_cis_associations_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env Rscript
    library(GenomicRanges)
    library(readr)
    library(dplyr)

    WINDOW <- ${params.cis_window}
    FDR_CUTOFF <- ${params.fdr_threshold}

    # ---- 1. Load gene catalog and split by role ----
    catalog <- read_tsv("${gene_catalog}", show_col_types = FALSE)
    lnc_df <- subset(catalog, origin %in% c("novel", "annotated_lncrna"))
    pc_df  <- subset(catalog, gene_biotype == "protein_coding")
    if (nrow(pc_df) == 0) {
        # Fallback: any non-lncRNA reference gene
        pc_df <- subset(catalog, origin == "reference_other")
        cat("WARNING: no protein_coding biotype found, using all non-lncRNA genes\\n")
    }
    cat("Catalog:", nrow(lnc_df), "lncRNAs (novel+annotated),", nrow(pc_df), "protein-coding targets\\n")

    if (nrow(lnc_df) == 0 || nrow(pc_df) == 0) {
        cat("WARNING: empty lncRNA or target set; writing empty output\\n")
        write_tsv(data.frame(), "${group}_cis_pairs.tsv")
        write_tsv(data.frame(), "${group}_cis_pairs_significant.tsv")
        sink("${group}_cis_associations_summary.tsv")
        cat("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n")
        cat("all\\t0\\t0\\tNA\\n")
        sink()
        quit(save = "no", status = 0)
    }

    # ---- 2. Genomic window search (gene-level) ----
    lnc_gr <- GRanges(
        seqnames = lnc_df\$chromosome,
        ranges   = IRanges(start = lnc_df\$start, end = lnc_df\$end),
        strand   = lnc_df\$strand
    )
    mcols(lnc_gr) <- DataFrame(
        gene_id  = lnc_df\$gene_id,
        origin   = lnc_df\$origin,
        class_code = lnc_df\$class_code
    )
    lnc_window <- lnc_gr + WINDOW

    pc_gr <- GRanges(
        seqnames = pc_df\$chromosome,
        ranges = IRanges(start = pc_df\$start, end = pc_df\$end),
        strand = pc_df\$strand
    )
    mcols(pc_gr) <- DataFrame(
        gene_id   = pc_df\$gene_id,
        gene_name = pc_df\$gene_name
    )

    ov <- findOverlaps(lnc_window, pc_gr)
    pairs <- data.frame(
        lncrna_id        = lnc_df\$gene_id[queryHits(ov)],
        lncrna_origin    = lnc_df\$origin[queryHits(ov)],
        lncrna_class_code = lnc_df\$class_code[queryHits(ov)],
        target_gene_id   = pc_df\$gene_id[subjectHits(ov)],
        target_gene_name = pc_df\$gene_name[subjectHits(ov)],
        stringsAsFactors = FALSE
    )
    pairs <- unique(pairs)

    n_pairs <- nrow(pairs)
    cat("Found", n_pairs, "lncRNA-gene pairs within +/-", WINDOW/1000, "kb\\n")

    if (n_pairs == 0) {
        cat("WARNING: no cis pairs found; writing empty output\\n")
        write_tsv(data.frame(), "${group}_cis_pairs.tsv")
        write_tsv(data.frame(), "${group}_cis_pairs_significant.tsv")
        sink("${group}_cis_associations_summary.tsv")
        cat("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n")
        cat("all\\t0\\t0\\tNA\\n")
        sink()
        quit(save = "no", status = 0)
    }

    # ---- 3. Load expression data (gene-level, full transcriptome) ----
    counts <- readRDS("${counts_rds}")
    lnc_ids <- unique(pairs\$lncrna_id)
    gene_ids <- unique(pairs\$target_gene_id)
    shared_lnc <- intersect(lnc_ids, rownames(counts))
    shared_gene <- intersect(gene_ids, rownames(counts))
    cat("Expression data available for", length(shared_lnc), "lncRNAs and",
        length(shared_gene), "target genes\\n")

    # ---- 5. Compute Spearman correlations ----
    log_cpm <- log2(counts + 1)
    expressed_ids <- rownames(log_cpm)

    # Pre-filter pairs to those with expression data on both sides
    pairs <- pairs[
        pairs\$lncrna_id %in% expressed_ids &
        pairs\$target_gene_id %in% expressed_ids, , drop = FALSE
    ]

    cis_list <- list()
    for (i in seq_len(nrow(pairs))) {
        lnc  <- pairs\$lncrna_id[i]
        gene <- pairs\$target_gene_id[i]

        lnc_expr  <- as.numeric(log_cpm[lnc, ])
        gene_expr <- as.numeric(log_cpm[gene, ])

        if (sd(lnc_expr) < 0.1 || sd(gene_expr) < 0.1) next

        cor_test <- tryCatch({
            cor.test(lnc_expr, gene_expr, method = "spearman")
        }, error = function(e) NULL)

        if (is.null(cor_test)) next

        cis_list[[length(cis_list) + 1L]] <- data.frame(
            lncrna_id      = lnc,
            lncrna_origin  = pairs\$lncrna_origin[i],
            lncrna_class_code = pairs\$lncrna_class_code[i],
            target_gene_id   = gene,
            target_gene_name = pairs\$target_gene_name[i],
            spearman_rho     = cor_test\$estimate,
            spearman_pvalue  = cor_test\$p.value,
            n_samples        = length(lnc_expr),
            stringsAsFactors = FALSE
        )
    }
    cis_results <- if (length(cis_list) > 0) {
        do.call(rbind, cis_list)
    } else {
        data.frame()
    }

    if (nrow(cis_results) == 0) {
        cat("WARNING: no valid correlation tests; writing empty output\\n")
        write_tsv(data.frame(), "${group}_cis_pairs.tsv")
        write_tsv(data.frame(), "${group}_cis_pairs_significant.tsv")
        sink("${group}_cis_associations_summary.tsv")
        cat("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n")
        cat("all\\t0\\t0\\tNA\\n")
        sink()
        quit(save = "no", status = 0)
    }

    # ---- 6. BH-FDR correction across full test family ----
    cis_results\$padj <- p.adjust(cis_results\$spearman_pvalue, method = "BH")
    cis_results <- cis_results[order(cis_results\$padj), ]
    cis_sig <- subset(cis_results, padj < FDR_CUTOFF)

    n_sig <- nrow(cis_sig)
    cat(sprintf("Cis associations: %d total, %d significant (FDR < %.2f)\\n",
                nrow(cis_results), n_sig, FDR_CUTOFF))

    # Annotate with DESeq2 results if available
    de <- tryCatch({
        read_tsv("${deseq2_summary}", show_col_types = FALSE)
    }, error = function(e) NULL)

    if (!is.null(de) && nrow(de) > 0 && "gene_id" %in% colnames(de)) {
        # DE is per-contrast; annotate the lncRNA side on the first match
        de_idx <- match(cis_sig\$lncrna_id, de\$gene_id)
        cis_sig\$lncrna_DE_log2FC <- de\$log2FoldChange[de_idx]
        cis_sig\$lncrna_DE_padj   <- de\$padj[de_idx]
    }

    write_tsv(cis_results, "${group}_cis_pairs.tsv")
    write_tsv(cis_sig, "${group}_cis_pairs_significant.tsv")

    sink("${group}_cis_associations_summary.tsv")
    cat(sprintf("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n"))
    cat(sprintf("all\\t%d\\t%d\\t%.4f\\n",
                nrow(cis_results), n_sig, median(cis_results\$spearman_rho, na.rm = TRUE)))
    sink()

    cat("DONE\\n")
    """
}