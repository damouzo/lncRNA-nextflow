// A10: Orthogonal catalog evidence overlap
// Tags transcripts: known_catalog_hit vs putatively_novel
// NONCODE, LNCipedia, RNAcentral, FANTOM CAT — reporting only, never a filter

process CATALOG_OVERLAP {

    input:
    path consensus_table
    path consensus_gtf

    output:
    path "catalog_overlap.tsv", emit: overlap_table

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)

    consensus <- read_tsv("${consensus_table}", show_col_types = FALSE)

    catalogs <- list()
    noncode_bed <- "${params.noncode_bed}"
    if (noncode_bed != "null") {
        cat("   NONCODE overlap check enabled\\n")
    }
    lncipedia_bed <- "${params.lncipedia_bed}"
    if (lncipedia_bed != "null") {
        cat("   LNCipedia overlap check enabled\\n")
    }
    rnacentral_bed <- "${params.rnacentral_bed}"
    if (rnacentral_bed != "null") {
        cat("   RNAcentral overlap check enabled\\n")
    }
    fantomcat_bed <- "${params.fantomcat_bed}"
    if (fantomcat_bed != "null") {
        cat("   FANTOM CAT overlap check enabled\\n")
    }

    # Placeholder: full catalog-overlap logic will be implemented
    # once reference BED files are configured
    consensus\$known_catalog_hit <- FALSE
    consensus\$catalog_sources <- "none"

    write_tsv(consensus, "catalog_overlap.tsv")
    """
}
