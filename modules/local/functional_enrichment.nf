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
    path "versions.yml", emit: versions

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
        # Same significance definition as the deliverable tables: padj < FDR and
        # |shrunken log2FC| >= LFC (the strict `significant` flag).
        fdr_cut <- ${params.fdr_threshold}
        lfc_cut <- ${params.lfc_threshold}
        de_sig <- ct_data\$gene_id[!is.na(ct_data\$padj) & ct_data\$padj < fdr_cut &
                                   abs(ct_data\$log2FoldChange) >= lfc_cut]
        foreground <- unique(intersect(targets, de_sig))
        foreground_type <- "significant_cis_targets"
        if (length(foreground) < 5) {
            # fallback: all cis targets expressed in this contrast
            foreground <- unique(intersect(targets, universe))
            foreground_type <- "all_cis_targets_expressed"
        }
        cat("Contrast", ct, ": ORA foreground =", length(foreground),
            "cis targets (", foreground_type, ")\\n")

        if (length(foreground) >= 5) {
            # pvalueCutoff/qvalueCutoff = 1: keep every tested term, the researcher
            # filters by p.adjust downstream. Nothing is dropped here.
            ego <- tryCatch({
                enrichGO(gene = foreground,
                         universe = universe,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         minGSSize = 5)
            }, error = function(e) NULL)
            if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                ego_df <- as.data.frame(ego)
                ego_df\$contrast <- ct
                ego_df\$foreground_type <- foreground_type
                ora_all[[ct]] <- ego_df
                cat(sprintf("  ORA: %d GO-BP terms (full table, no cutoff)\n", nrow(ego_df)))
            } else {
                cat(sprintf("  ORA: enrichGO returned 0 terms (foreground=%d, universe=%d)\n",
                            length(foreground), length(universe)))
            }
        } else {
            cat(sprintf("  ORA: insufficient foreground (%d genes, need >= 5)\n",
                        length(foreground)))
        }

        # ---- GSEA: all protein-coding genes ranked by the unshrunken DE stat ----
        gsea <- tryCatch({
            gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1)
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

    # write_tsv() on a 0-column data.frame writes a 0-byte file (no header), which
    # is indistinguishable from a failure. Build a 0-row frame with the expected
    # clusterProfiler columns instead, so an empty result is a valid table.
    empty_frame <- function(cols) {
        d <- data.frame()
        for (col in cols) d[[col]] <- character(0)
        d
    }
    ora_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                  "qvalue", "geneID", "Count", "contrast", "foreground_type")
    gsea_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue",
                   "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment", "contrast")

    ora_df <- if (length(ora_all) > 0) as.data.frame(do.call(rbind, ora_all)) else empty_frame(ora_cols)
    gsea_df <- if (length(gsea_all) > 0) as.data.frame(do.call(rbind, gsea_all)) else empty_frame(gsea_cols)

    write_tsv(ora_df, "${group}_enrichment_ora.tsv")
    write_tsv(gsea_df, "${group}_enrichment_gsea.tsv")

    ver <- paste0("R ", paste(R.version\$major, R.version\$minor, sep = "."))
    for (pkg in c("clusterProfiler", "org.Hs.eg.db")) {
        pv <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA)
        if (!is.na(pv)) ver <- paste0(ver, "; ", pkg, " ", pv)
    }
    writeLines(c(sprintf('"%s":', "${task.process}"),
                 '    tool: "clusterProfiler"',
                 sprintf('    version: "%s"', ver),
                 sprintf('    container: "%s"', "${task.container}")),
               "versions.yml")
    """
}