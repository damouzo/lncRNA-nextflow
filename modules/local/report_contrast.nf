// B7: Per-contrast HTML reporting (Quarto, self-contained)
// Renders one self-contained HTML report per contrast using the frozen catalog,
// DESeq2 results, cis associations and functional enrichment of a norm_group.
// Reports are produced for every contrast declared in the input contrasts file
// that applies to this norm_group (numerator/denominator conditions present in
// the group's coldata), even when no results were obtained — empty analyses are
// reported as such so absence of results can be told apart from pipeline errors.
// Output is a flat set of HTML files, one per contrast, landing directly under
// results/reports (no per-group subdirectory).

process REPORT_CONTRAST {

    // Physical copies instead of symlinks: quarto resolves the real path of its
    // input document, so a symlinked template would render next to the source
    // file in assets/ instead of in the task dir.
    stageInMode 'copy'

    input:
    tuple val(group), path(de_summary), path(cis_pairs),
          path(cis_pairs_sig), path(ora), path(gsea),
          path(heterogeneity), path(heterogeneity_cv),
          path(coldata), path(contrasts_tsv), path(tpm_rds)
    path template         // report_template.qmd
    path gene_catalog     // gene_catalog.tsv (origin/class_code/biotype annotations)

    output:
    path "*.html", emit: reports

    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)

    de <- read_tsv("${de_summary}", show_col_types = FALSE)

    # Condition levels present in this norm_group's coldata
    group_conditions <- tryCatch({
        unique(na.omit(read.csv("${coldata}", stringsAsFactors = FALSE)[["condition"]]))
    }, error = function(e) character())

    # Intended contrasts: the user declared them upstream; keep only those that
    # apply to this group (both conditions present in the coldata)
    ct_intended <- tryCatch({
        read_tsv("${contrasts_tsv}", show_col_types = FALSE,
                 col_names = c("contrast_name", "numerator", "denominator", "batch")) %>%
            filter(numerator %in% group_conditions,
                   denominator %in% group_conditions) %>%
            pull(contrast_name) %>%
            unique(.)
    }, error = function(e) character())

    # Fall back to whatever made it into the DE summary if the input was unusable
    contrasts <- if (length(ct_intended) > 0) {
        ct_intended
    } else if ("contrast" %in% colnames(de)) {
        unique(na.omit(de[["contrast"]]))
    } else {
        "all"
    }

    for (ct in contrasts) {
        fname <- paste0(gsub("[^A-Za-z0-9_.-]", "_", ct), "_report.html")
        if (file.exists(fname)) next
        cmd <- paste(
            "quarto render ${template}",
            "-M", shQuote(paste0("title=lncRNA report - ", ct)),
            "-M", shQuote(paste0("subtitle=norm_group: ${group}")),
            "-P", shQuote(paste0("contrast=", ct)),
            "-P", shQuote(paste0("group=${group}")),
            "-P", shQuote(paste0("de_summary=${de_summary}")),
            "-P", shQuote(paste0("cis_pairs=${cis_pairs}")),
            "-P", shQuote(paste0("cis_pairs_sig=${cis_pairs_sig}")),
            "-P", shQuote(paste0("cis_window=${params.cis_window}")),
            "-P", shQuote(paste0("ora=${ora}")),
            "-P", shQuote(paste0("gsea=${gsea}")),
            "-P", shQuote(paste0("heterogeneity=${heterogeneity}")),
            "-P", shQuote(paste0("heterogeneity_cv=${heterogeneity_cv}")),
            "-P", shQuote(paste0("heterogeneity_concordant_frac=${params.heterogeneity_concordant_frac}")),
            "-P", shQuote(paste0("heterogeneity_discordant_frac=${params.heterogeneity_discordant_frac}")),
            "-P", shQuote(paste0("gene_catalog=${gene_catalog}")),
            "-P", shQuote(paste0("species=${params.species}")),
            "-P", shQuote(paste0("genome_build=${params.genome_build ?: ''}")),
            "-P", shQuote(paste0("conservation_metric=${params.conservation_metric}")),
            "-P", shQuote(paste0("conservation_configured=${params.conservation_bigwig ? true : false}")),
            "-P", shQuote(paste0("synteny_configured=${(params.synteny_chain_file && params.synteny_target_gtf) ? true : false}")),
            "-P", shQuote(paste0("synteny_target_species=${params.synteny_target_species ?: ''}")),
            "-P", shQuote(paste0("tpm=${tpm_rds}")),
            "--output", shQuote(fname)
        )
        cat("Rendering contrast:", ct, "\n")
        status <- system(cmd)
        if (status != 0) stop(paste("quarto render failed for contrast:", ct))
    }
    """
}
