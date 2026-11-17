// A9: Coding-potential consensus — CPAT ∩ CPC2 agreement
// Requires both tools to agree for "high-confidence non-coding" status

process CODING_CONSENSUS {

    input:
    path cpat_tsv
    path cpc2_txt
    path filtered_gtf

    output:
    path "coding_consensus.tsv", emit: consensus_table
    path "consensus_noncoding.gtf", emit: consensus_gtf

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)

    cpat <- read_tsv("${cpat_tsv}", show_col_types = FALSE)
    cpc2 <- read.table("${cpc2_txt}", header = TRUE, sep = "\\t",
                       stringsAsFactors = FALSE, comment.char = "")

    colnames(cpc2)[1] <- "ID"

    # Keep every CPC2 call (not only the non-coding subset) so the CPAT-vs-CPC2
    # classifier figure has a numeric probability per transcript. The non-coding
    # call is still the threshold comparison, as before.
    cpc2 <- cpc2[, intersect(c("ID", "coding_probability"), colnames(cpc2)), drop = FALSE]
    cpc2\$coding_probability <- suppressWarnings(as.numeric(cpc2\$coding_probability))

    cpat_cols <- intersect(c("ID", "mRNA", "cpat_coding_prob"), colnames(cpat))
    merged <- merge(cpat[, cpat_cols, drop = FALSE], cpc2, by = "ID", all = TRUE)
    merged\$CPAT_noncoding <- ifelse(is.na(merged\$mRNA), FALSE, merged\$mRNA == "no")
    merged\$CPC2_noncoding <- !is.na(merged\$coding_probability) &
        merged\$coding_probability < ${params.cpc2_coding_threshold}
    merged\$consensus_noncoding <- merged\$CPAT_noncoding & merged\$CPC2_noncoding

    write_tsv(merged, "coding_consensus.tsv")

    # Subset GTF to transcripts with consensus non-coding status
    ids_keep <- merged\$ID[merged\$consensus_noncoding]
    gtf_lines <- read_lines("${filtered_gtf}")
    has_tx <- grepl('transcript_id "', gtf_lines)
    tx_id <- sub('.*transcript_id "([^"]+)".*', '\\\\1', gtf_lines[has_tx])
    write_lines(gtf_lines[has_tx][tx_id %in% ids_keep], "consensus_noncoding.gtf")
    """
}
