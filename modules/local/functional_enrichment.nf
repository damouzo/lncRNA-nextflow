// B6: Functional enrichment (ORA + GSEA via clusterProfiler)
// Uses the experiment's own expressed-gene universe as background

process FUNCTIONAL_ENRICHMENT {
    label 'medium_task'
    memory = 16.GB

    input:
    path cis_pairs_sig
    path deseq2_summary

    output:
    path "enrichment_ora.tsv", emit: ora_results
    path "enrichment_gsea.tsv", emit: gsea_results
    path "enrichment_plots/*", emit: plots

    script:
    """
    #!/usr/bin/env Rscript
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(readr)
    library(dplyr)
    library(ggplot2)

    dir.create("enrichment_plots", showWarnings = FALSE)

    de_results <- read_tsv("${deseq2_summary}", show_col_types = FALSE)

    # Separate lncRNA targets by contrast for enrichment
    if ("contrast" %in% colnames(de_results)) {
        contrasts <- unique(de_results\\$contrast)
    } else {
        contrasts <- "all"
    }

    ora_all <- list()
    gsea_all <- list()

    for (ct in contrasts) {
        ct_data <- de_results
        if (ct != "all") {
            ct_data <- subset(de_results, contrast == ct)
        }

        # Use ranked list for GSEA
        if ("stat" %in% colnames(ct_data)) {
            gene_list <- ct_data\\$stat
            names(gene_list) <- ct_data\\$gene_id
            gene_list <- sort(gene_list, decreasing = TRUE)
        } else {
            gene_list <- ct_data\\$log2FoldChange
            names(gene_list) <- ct_data\\$gene_id
            gene_list <- sort(gene_list, decreasing = TRUE)
        }

        # ORA: enriched GO terms
        sig_genes <- subset(ct_data, padj < 0.05)\\$gene_id
        if (length(sig_genes) > 0) {
            ego <- tryCatch({
                enrichGO(gene = sig_genes,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         universe = unique(ct_data\\$gene_id))
            }, error = function(e) NULL)

            if (!is.null(ego)) {
                ora_all[[ct]] <- as.data.frame(ego)
            }
        }

        # GSEA
        gsea <- tryCatch({
            gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL",
                  ont = "BP",
                  pAdjustMethod = "BH")
        }, error = function(e) NULL)

        if (!is.null(gsea)) {
            gsea_df <- as.data.frame(gsea)
            gsea_df\\$contrast <- ct
            gsea_all[[ct]] <- gsea_df

            # Save dotplot
            png(file.path("enrichment_plots", paste0(ct, "_gsea_dotplot.png")),
                width = 800, height = 600)
            print(dotplot(gsea, showCategory = 20))
            dev.off()
        }
    }

    ora_df <- if (length(ora_all) > 0) do.call(rbind, ora_all) else data.frame()
    gsea_df <- if (length(gsea_all) > 0) do.call(rbind, gsea_all) else data.frame()

    write_tsv(as.data.frame(ora_df), "enrichment_ora.tsv")
    write_tsv(as.data.frame(gsea_df), "enrichment_gsea.tsv")
    """
}
