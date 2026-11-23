// B6: Functional enrichment (ORA via clusterProfiler).
// One over-representation analysis on coding genes (GO-BP): the coding neighbours
// of the lncRNAs that are DE-significant in the contrast, using every cis pair
// (significant or not). MSTRG genes have no functional annotation, so the lncRNA
// signal is read through the coding neighbours they are associated with.
// Background = coding genes in the window of any tested lncRNA (the cis family),
// not the whole transcriptome. No direction is assumed, so it answers "which
// processes do the neighbours of the changing lncRNAs fall in?".
// Fewer than 5 usable genes: the analysis is not run (the report says so).

process FUNCTIONAL_ENRICHMENT {

    input:
    tuple val(group), path(cis_pairs_all), path(deseq2_summary)
    path gene_catalog

    output:
    tuple val(group), path("${group}_enrichment_ora.tsv"), emit: ora_results
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

    dir.create("${group}_enrichment_plots", showWarnings = FALSE)

    catalog <- read_tsv("${gene_catalog}", show_col_types = FALSE)
    pc_ids     <- unique(catalog\$gene_id[catalog\$gene_biotype == "protein_coding"])
    lncrna_ids <- unique(catalog\$gene_id[catalog\$origin %in% c("novel", "annotated_lncrna")])

    # The cis family: every coding gene in the window of a tested lncRNA. Feeds both
    # the neighbourhood foreground and its background.
    cis_all <- read_tsv("${cis_pairs_all}", show_col_types = FALSE)
    neigh_universe <- unique(cis_all\$target_gene_id[!is.na(cis_all\$target_gene_id)])

    # An antisense-overlapping lncRNA can be DE only because its host leaks reads, so
    # its host would enter the foreground and the ORA would just echo the host's DE.
    # Drop leak-suspect pairs from the foreground (the background keeps the full
    # family). leak_suspect is present in the full cis_pairs table.
    cis_fg <- if ("leak_suspect" %in% colnames(cis_all)) {
        cis_all[!(cis_all\$leak_suspect %in% TRUE), , drop = FALSE]
    } else {
        cis_all
    }

    de_results <- read_tsv("${deseq2_summary}", show_col_types = FALSE)
    contrasts <- if ("contrast" %in% colnames(de_results)) unique(de_results\$contrast) else "all"
    fdr_cut <- ${params.fdr_threshold}

    label <- "de_lncrna_neighbours"
    ora_all <- list()

    for (ct in contrasts) {
        ct_data <- if (ct == "all") de_results else subset(de_results, contrast == ct)
        if (nrow(ct_data) == 0) next

        # Universe = protein-coding genes tested in this contrast.
        pc_data  <- ct_data[ct_data\$gene_id %in% pc_ids, ]
        universe <- unique(pc_data\$gene_id)
        if (length(universe) < 10) {
            cat("Contrast", ct, ": too few protein-coding genes tested, skipping\\n")
            next
        }

        # Foreground: coding neighbours of the DE-significant lncRNAs (any cis pair,
        # leak-suspect pairs dropped).
        de_lnc <- ct_data\$gene_id[!is.na(ct_data\$padj) & ct_data\$padj < fdr_cut &
                                   ct_data\$gene_id %in% lncrna_ids]
        fg <- intersect(unique(cis_fg\$target_gene_id[cis_fg\$lncrna_id %in% de_lnc]),
                        universe)
        bg <- intersect(neigh_universe, universe)
        if (length(fg) < 5) {
            cat(sprintf("Contrast %s: ORA skipped, %d usable neighbours (need >= 5)\\n",
                        ct, length(fg)))
            next
        }

        # pvalueCutoff/qvalueCutoff = 1: keep every tested term, the researcher
        # filters by p.adjust downstream. Nothing is dropped here.
        ego <- tryCatch({
            enrichGO(gene = fg, universe = bg, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                     ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                     minGSSize = 5)
        }, error = function(e) NULL)
        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
            cat(sprintf("  ORA: enrichGO returned 0 terms (fg=%d, bg=%d)\\n", length(fg), length(bg)))
            next
        }
        d <- as.data.frame(ego)
        d\$contrast <- ct
        d\$collection <- "GO_BP"
        # input_set/foreground_type: the exact gene list that entered the test, so
        # the table is self-describing in the deliverables and the report.
        d\$input_set <- label
        d\$foreground_type <- label
        png(file.path("${group}_enrichment_plots", paste0(ct, "_ora_dotplot.png")),
            width = 900, height = 700)
        print(dotplot(ego, showCategory = 20))
        dev.off()
        pdf(file.path("${group}_enrichment_plots", paste0(ct, "_ora_dotplot.pdf")),
            width = 9, height = 7)
        print(dotplot(ego, showCategory = 20))
        dev.off()
        cat(sprintf("  ORA: %d GO-BP terms (full table, no cutoff)\\n", nrow(d)))
        ora_all[[ct]] <- d
    }

    # write_tsv() on a 0-column data.frame writes a 0-byte file (no header), which
    # is indistinguishable from a failure. Build a 0-row frame with the expected
    # clusterProfiler columns instead, so an empty result is a valid table.
    empty_frame <- function(cols) {
        d <- data.frame()
        for (col in cols) d[[col]] <- character(0)
        d
    }
    ora_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment",
                  "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count",
                  "contrast", "collection", "input_set", "foreground_type")
    ora_df <- if (length(ora_all) > 0) as.data.frame(do.call(rbind, ora_all)) else empty_frame(ora_cols)

    write_tsv(ora_df, "${group}_enrichment_ora.tsv")

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
