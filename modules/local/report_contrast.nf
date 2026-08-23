// B7: Per-contrast HTML reporting (Quarto, self-contained)
// Renders one self-contained HTML report per contrast using the frozen catalog,
// DESeq2 results, cis associations and functional enrichment of a norm_group.

process REPORT_CONTRAST {

    input:
    tuple val(group), path(de_summary), path(cis_pairs),
          path(cis_pairs_sig), path(ora), path(gsea)
    path template         // report_template.qmd

    output:
    path "${group}_reports/*.html", emit: reports

    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)

    de <- read_tsv("${de_summary}", show_col_types = FALSE)
    contrasts <- if ("contrast" %in% colnames(de)) {
        unique(na.omit(de[["contrast"]]))
    } else {
        "all"
    }

    dir.create("${group}_reports", showWarnings = FALSE)

    for (ct in contrasts) {
        fname <- paste0(gsub("[^A-Za-z0-9_.-]", "_", ct), "_report.html")
        out <- file.path("${group}_reports", fname)
        if (file.exists(out)) next
        cmd <- paste(
            "quarto render ${template}",
            "-M", shQuote(paste0("title=lncRNA report - ", ct)),
            "-M", shQuote(paste0("subtitle=norm_group: ${group}")),
            "-P", shQuote(paste0("contrast=", ct)),
            "-P", shQuote(paste0("group=${group}")),
            "-P", shQuote(paste0("de_summary=${de_summary}")),
            "-P", shQuote(paste0("cis_pairs=${cis_pairs}")),
            "-P", shQuote(paste0("cis_pairs_sig=${cis_pairs_sig}")),
            "-P", shQuote(paste0("ora=${ora}")),
            "-P", shQuote(paste0("gsea=${gsea}")),
            "--output-dir", shQuote("${group}_reports"),
            "--output", shQuote(fname)
        )
        cat("Rendering contrast:", ct, "\n")
        status <- system(cmd)
        if (status != 0) stop(paste("quarto render failed for contrast:", ct))
    }
    """
}
