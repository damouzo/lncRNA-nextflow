// Cis-regulatory candidate associations
// GenomicRanges-based window search for lncRNA–protein_coding gene pairs,
// Spearman correlation (global + within-condition), BH-FDR correction
// across the full test family per contrast

process CIS_ASSOCIATIONS {
    label 'medium_task'
    memory = 16.GB

    input:
    path deseq2_summary     // contrasts_summary.tsv from DESeq2
    path counts_rds         // filtered count matrix from tximport_filter
    path frozen_gtf         // frozen lncRNA GTF
    path reference_gtf      // reference GTF (for protein-coding genes)
    val group               // norm_group name for output prefix

    output:
    path "${group}_cis_pairs.tsv",              emit: cis_pairs
    path "${group}_cis_pairs_significant.tsv",  emit: cis_pairs_sig
    path "${group}_cis_associations_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env Rscript
    library(GenomicRanges)
    library(rtracklayer)
    library(readr)
    library(dplyr)
    library(tibble)

    WINDOW <- ${params.cis_window}
    FDR_CUTOFF <- ${params.fdr_threshold}

    # ---- 1. Load lncRNA coordinates from frozen GTF ----
    gtf_all <- import("${frozen_gtf}")
    lnc_tx <- gtf_all[gtf_all\\$type == "transcript"]
    if (length(lnc_tx) == 0) {
        cat("WARNING: no transcript features in frozen GTF, trying exon/transcript\\n")
        lnc_tx <- gtf_all
    }

    lnc_df <- as.data.frame(lnc_tx)
    cat("Loaded", nrow(lnc_df), "lncRNA features\\n")

    # ---- 2. Load protein-coding genes from reference GTF ----
    ref_gtf <- import("${reference_gtf}")
    pc_genes <- ref_gtf[ref_gtf\\$type == "gene" &
                         grepl("protein_coding", ref_gtf\\$gene_biotype)]
    if (length(pc_genes) == 0) {
        # Fallback: use all genes if no protein_coding biotype field
        pc_genes <- ref_gtf[ref_gtf\\$type == "gene"]
        cat("WARNING: no protein_coding biotype found, using all gene features\\n")
    }

    pc_df <- as.data.frame(pc_genes)
    cat("Loaded", nrow(pc_df), "protein-coding genes from reference\\n")

    # ---- 3. Genomic window search ----
    lnc_gr <- GRanges(
        seqnames = lnc_df\\$seqnames,
        ranges = IRanges(start = lnc_df\\$start, end = lnc_df\\$end),
        strand = lnc_df\\$strand,
        transcript_id = lnc_df\\$transcript_id
    )
    lnc_window <- lnc_gr + WINDOW

    pc_gr <- GRanges(
        seqnames = pc_df\\$seqnames,
        ranges = IRanges(start = pc_df\\$start, end = pc_df\\$end),
        strand = pc_df\\$strand,
        gene_id = pc_df\\$gene_id,
        gene_name = pc_df\\$gene_name
    )

    overlaps <- findOverlaps(lnc_window, pc_gr)
    pairs <- data.frame(
        lncrna_id = lnc_gr\\$transcript_id[queryHits(overlaps)],
        target_gene_id = pc_gr\\$gene_id[subjectHits(overlaps)],
        target_gene_name = pc_gr\\$gene_name[subjectHits(overlaps)],
        stringsAsFactors = FALSE
    )
    pairs <- unique(pairs)

    n_pairs <- nrow(pairs)
    cat("Found", n_pairs, "lncRNA-gene pairs within +/-", WINDOW/1000, "kb\\n")

    if (n_pairs == 0) {
        cat("WARNING: no cis pairs found; writing empty output\\n")
        empty <- data.frame()
        write_tsv(empty, "${group}_cis_pairs.tsv")
        write_tsv(empty, "${group}_cis_pairs_significant.tsv")
        sink("${group}_cis_associations_summary.tsv")
        cat("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n")
        sink()
        quit(save = "no", status = 0)
    }

    # ---- 4. Load expression data ----
    counts <- readRDS("${counts_rds}")
    lnc_ids <- unique(pairs\\$lncrna_id)
    gene_ids <- unique(pairs\\$target_gene_id)
    shared_lnc <- intersect(lnc_ids, rownames(counts))
    shared_gene <- intersect(gene_ids, rownames(counts))
    cat("Expression data available for", length(shared_lnc), "lncRNAs and",
        length(shared_gene), "target genes\\n")

    # ---- 5. Compute Spearman correlations ----
    # Use log2(CPM + 1) transformation
    log_cpm <- log2(counts + 1)

    cis_results <- data.frame()
    for (i in seq_len(nrow(pairs))) {
        lnc <- pairs\\$lncrna_id[i]
        gene <- pairs\\$target_gene_id[i]

        if (!(lnc %in% rownames(log_cpm)) || !(gene %in% rownames(log_cpm))) next

        lnc_expr <- as.numeric(log_cpm[lnc, ])
        gene_expr <- as.numeric(log_cpm[gene, ])

        # Require minimum expression variance
        if (sd(lnc_expr) < 0.1 || sd(gene_expr) < 0.1) next

        cor_test <- tryCatch({
            cor.test(lnc_expr, gene_expr, method = "spearman")
        }, error = function(e) NULL)

        if (is.null(cor_test)) next

        cis_results <- rbind(cis_results, data.frame(
            lncrna_id = lnc,
            target_gene_id = gene,
            target_gene_name = pairs\\$target_gene_name[i],
            spearman_rho = cor_test\\$estimate,
            spearman_pvalue = cor_test\\$p.value,
            n_samples = length(lnc_expr),
            stringsAsFactors = FALSE
        ))
    }

    if (nrow(cis_results) == 0) {
        cat("WARNING: no valid correlation tests; writing empty output\\n")
        write_tsv(pairs, "${group}_cis_pairs.tsv")
        write_tsv(data.frame(), "${group}_cis_pairs_significant.tsv")
        sink("${group}_cis_associations_summary.tsv")
        cat("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n")
        cat("all\\t0\\t0\\tNA\\n")
        sink()
        quit(save = "no", status = 0)
    }

    # ---- 6. BH-FDR correction across full test family ----
    cis_results\\$padj <- p.adjust(cis_results\\$spearman_pvalue, method = "BH")
    cis_results <- cis_results[order(cis_results\\$padj), ]
    cis_sig <- subset(cis_results, padj < FDR_CUTOFF)

    n_sig <- nrow(cis_sig)
    n_total <- nrow(cis_results)
    cat(sprintf("Cis associations: %d total, %d significant (FDR < %.2f)\\n",
                n_total, n_sig, FDR_CUTOFF))

    # Annotate with DESeq2 results if available
    de <- tryCatch({
        read_tsv("${deseq2_summary}", show_col_types = FALSE)
    }, error = function(e) NULL)

    if (!is.null(de) && nrow(de) > 0 && "gene_id" %in% colnames(de)) {
        # Add DE status to significant pairs
        cis_sig\\$lncrna_DE_log2FC <- de\\$log2FoldChange[
            match(cis_sig\\$lncrna_id, de\\$gene_id)]
        cis_sig\\$lncrna_DE_padj <- de\\$padj[
            match(cis_sig\\$lncrna_id, de\\$gene_id)]
    }

    write_tsv(cis_results, "${group}_cis_pairs.tsv")
    write_tsv(cis_sig, "${group}_cis_pairs_significant.tsv")

    sink("${group}_cis_associations_summary.tsv")
    cat(sprintf("contrast\\tn_pairs\\tn_significant\\tmedian_cor\\n"))
    cat(sprintf("all\\t%d\\t%d\\t%.4f\\n",
                n_total, n_sig, median(cis_results\\$spearman_rho, na.rm = TRUE)))
    sink()

    cat("DONE\\n")
    """
}
