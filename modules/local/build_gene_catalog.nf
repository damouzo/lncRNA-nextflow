// Gene catalog of the full analysis universe
// One row per gene_id, tagging each as:
//   novel            — frozen discovery lncRNAs (with gffcompare class_code)
//   annotated_lncrna — reference genes whose gene_biotype is in the configurable list
//   reference_other  — everything else in the reference annotation
// Used downstream (reporting, cis, enrichment) without touching quantification.

process BUILD_GENE_CATALOG {

    input:
    path frozen_gtf         // frozen novel lncRNA GTF (class_code in attributes)
    path reference_gtf      // reference annotation

    output:
    path "gene_catalog.tsv", emit: gene_catalog
    path "catalog_stats.tsv", emit: catalog_stats

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)

    ANNOT_BIOTYPES <- c(${params.annotated_lncrna_biotypes.collect { "\"$it\"" } .join(", ")})

    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }
    strip_version <- function(x) sub("\\\\.[0-9]+\$", "", x)

    fill_table <- function(lines) {
        if (length(lines) == 0) {
            return(data.frame(
                chromosome = character(),
                start      = integer(),
                end        = integer(),
                strand     = character(),
                attrs      = character(),
                stringsAsFactors = FALSE
            ))
        }
        parts <- strsplit(lines, "\\t", fixed = TRUE)
        idx <- lengths(parts) >= 9L
        if (!any(idx)) stop("No parseable lines in this GTF subset")
        parts <- parts[idx]
        data.frame(
            chromosome = vapply(parts, `[`, character(1), 1L),
            start      = as.integer(vapply(parts, `[`, character(1), 4L)),
            end        = as.integer(vapply(parts, `[`, character(1), 5L)),
            strand     = vapply(parts, `[`, character(1), 7L),
            attrs      = vapply(parts, `[`, character(1), 9L),
            stringsAsFactors = FALSE
        )
    }

    # ---- Reference genes (one line per gene) ----
    system(paste("grep -P '\\tgene\\t'", shQuote("${reference_gtf}"),
                 "> ref_gene_lines.txt"))
    ref <- fill_table(readLines("ref_gene_lines.txt"))
    ref\$gene_id      <- strip_version(extract_attr(ref\$attrs, "gene_id"))
    ref\$gene_name    <- extract_attr(ref\$attrs, "gene_name")
    ref\$gene_biotype <- extract_attr(ref\$attrs, "gene_biotype")

    ref_cat <- ref %>%
        filter(!is.na(gene_id)) %>%
        mutate(
            origin     = ifelse(gene_biotype %in% ANNOT_BIOTYPES, "annotated_lncrna", "reference_other"),
            class_code = NA_character_
        ) %>%
        distinct(gene_id, .keep_all = TRUE) %>%
        select(gene_id, gene_name, chromosome, start, end, strand,
               gene_biotype, origin, class_code)

    # ---- Novel lncRNAs (transcript-level, aggregated to gene) ----
    # MSTRG IDs are not versioned like Ensembl — keep them verbatim
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${frozen_gtf}"),
                 "> novel_transcript_lines.txt"))
    nov_lines <- readLines("novel_transcript_lines.txt")
    if (length(nov_lines) == 0) {
        cat("WARNING: no novel transcripts parsed from frozen GTF, novel origin will be empty\\n")
    }
    nov <- fill_table(nov_lines)
    nov\$gene_id <- extract_attr(nov\$attrs, "gene_id")
    nov\$class_code <- extract_attr(nov\$attrs, "class_code")

    nov_cat <- nov %>%
        filter(!is.na(gene_id)) %>%
        group_by(gene_id) %>%
        summarise(
            chromosome = dplyr::first(chromosome),
            start      = min(start),
            end        = max(end),
            strand     = dplyr::first(strand),
            class_code = paste(sort(unique(na.omit(class_code))), collapse = ","),
            .groups    = "drop"
        ) %>%
        mutate(
            origin       = "novel",
            gene_biotype = NA_character_,
            gene_name    = NA_character_
        ) %>%
        select(gene_id, gene_name, chromosome, start, end, strand,
               gene_biotype, origin, class_code)

    # ---- Merge: novel wins if a gene_id ever collides with reference ----
    catalog <- bind_rows(nov_cat, ref_cat) %>%
        distinct(gene_id, .keep_all = TRUE) %>%
        arrange(origin, gene_id)

    if (nrow(catalog) == 0) stop("Gene catalog is empty — check reference GTF parsing")

    write_tsv(catalog, "gene_catalog.tsv")

    stats <- catalog %>%
        count(origin, name = "n_genes") %>%
        bind_rows(
            catalog %>% filter(!is.na(gene_biotype)) %>%
                count(gene_biotype, name = "n_genes") %>%
                rename(origin = gene_biotype)
        )
    write_tsv(stats, "catalog_stats.tsv")

    cat("Catalog genes:", nrow(catalog), "\\n")
    print(table(catalog\$origin, useNA = "ifany"))
    """
}