// B6: Functional enrichment (ORA + GSEA via clusterProfiler)
// Runs over protein-coding targets of significant cis pairs (guilt-by-association).
// MSTRG genes have no functional annotation, so the lncRNA signal is mapped
// through the coding neighbors they are correlated with.
// Background = experiment's expressed protein-coding universe (per contrast).

process FUNCTIONAL_ENRICHMENT {

    input:
    tuple val(group), path(cis_pairs_sig), path(deseq2_summary)  // keyed: no positional pairing downstream
    path gene_catalog

    output:
    tuple val(group), path("${group}_enrichment_ora.tsv"),  emit: ora_results
    tuple val(group), path("${group}_enrichment_gsea.tsv"), emit: gsea_results
    path "${group}_enrichment_plots/*", optional: true, emit: plots

    script:
    """
    #!/usr/bin/env Rscript
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(readr)
    library(dplyr)
    library(ggplot2)

    dir.create("${group}_enrichment_plots", showWarnings = FALSE)

    catalog <- read_tsv("${gene_catalog}", show_col_types = FALSE)
    pc_ids <- unique(catalog\$gene_id[catalog\$gene_biotype == "protein_coding"])

    cis_sig <- tryCatch({
        read_tsv("${cis_pairs_sig}", show_col_types = FALSE)
    }, error = function(e) NULL)
    if (is.null(cis_sig) || nrow(cis_sig) == 0 ||
        !("target_gene_id" %in% colnames(cis_sig))) {
        cis_sig <- NULL
    } else {
        cis_sig <- cis_sig[!is.na(cis_sig\$target_gene_id), ]
    }
    targets <- if (!is.null(cis_sig)) unique(cis_sig\$target_gene_id) else character()

    de_results <- read_tsv("${deseq2_summary}", show_col_types = FALSE)
    contrasts <- if ("contrast" %in% colnames(de_results)) {
        unique(de_results\$contrast)
    } else {
        "all"
    }

    ora_all <- list()
    gsea_all <- list()

    for (ct in contrasts) {
        ct_data <- de_results
        if (ct != "all") {
            ct_data <- subset(de_results, contrast == ct)
        }
        if (nrow(ct_data) == 0) next

        # Restrict to protein-coding genes tested in this contrast (universe)
        pc_data <- ct_data[ct_data\$gene_id %in% pc_ids, ]
        universe <- unique(pc_data\$gene_id)
        if (length(universe) < 10) {
            cat("Contrast", ct, ": too few protein-coding genes tested, skipping\\n")
            next
        }

        # Ranked list for GSEA (protein-coding only)
        stat_col <- if ("stat" %in% colnames(pc_data)) "stat" else "log2FoldChange"
        gene_list <- as.numeric(pc_data[[stat_col]])
        names(gene_list) <- pc_data\$gene_id
        gene_list <- sort(gene_list[is.finite(gene_list)], decreasing = TRUE)

        # ---- ORA: DE-significant cis targets (guilt by association) ----
        de_sig <- ct_data\$gene_id[!is.na(ct_data\$padj) & ct_data\$padj < 0.05]
        foreground <- unique(intersect(targets, de_sig))
        if (length(foreground) < 5) {
            # fallback: all cis targets expressed in this contrast
            foreground <- unique(intersect(targets, universe))
        }
        cat("Contrast", ct, ": ORA foreground =", length(foreground), "cis targets\\n")

        if (length(foreground) >= 5) {
            ego <- tryCatch({
                enrichGO(gene = foreground,
                         universe = universe,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,
                         minGSSize = 5)
            }, error = function(e) NULL)
            if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                ego_df <- as.data.frame(ego)
                ego_df\$contrast <- ct
                ora_all[[ct]] <- ego_df
                cat(sprintf("  ORA: %d GO-BP terms (p<0.1)\n", nrow(ego_df)))
            } else {
                cat(sprintf("  ORA: enrichGO returned 0 terms (foreground=%d, universe=%d)\n",
                            length(foreground), length(universe)))
            }
        } else {
            cat(sprintf("  ORA: insufficient foreground (%d genes, need >= 5)\n",
                        length(foreground)))
        }

        # ---- GSEA: all protein-coding genes ranked by DE stat ----
        gsea <- tryCatch({
            gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL",
                  ont = "BP",
                  pAdjustMethod = "BH")
        }, error = function(e) NULL)

        if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
            gsea_df <- as.data.frame(gsea)
            gsea_df\$contrast <- ct
            gsea_all[[ct]] <- gsea_df

            png(file.path("${group}_enrichment_plots", paste0(ct, "_gsea_dotplot.png")),
                width = 800, height = 600)
            print(dotplot(gsea, showCategory = 20))
            dev.off()
        }
    }

    ora_df <- if (length(ora_all) > 0) do.call(rbind, ora_all) else data.frame()
    gsea_df <- if (length(gsea_all) > 0) do.call(rbind, gsea_all) else data.frame()

    write_tsv(as.data.frame(ora_df), "${group}_enrichment_ora.tsv")
    write_tsv(as.data.frame(gsea_df), "${group}_enrichment_gsea.tsv")
    """
}