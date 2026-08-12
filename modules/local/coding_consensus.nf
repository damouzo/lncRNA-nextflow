// A9: Coding-potential consensus — CPAT ∩ CPC2 agreement
// Requires both tools to agree for "high-confidence non-coding" status

process CODING_CONSENSUS {
    label 'small_task'

    input:
    path cpat_tsv
    path cpc2_txt

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
                       stringsAsFactors = FALSE)

    colnames(cpc2)[1] <- "ID"
    cpc2_noncode <- subset(cpc2, coding_probability < 0.5)

    merged <- merge(
        cpat[, c("ID", "mRNA")],
        cpc2_noncode[, c("ID", "coding_probability")],
        by = "ID",
        all = TRUE
    )
    merged\$CPAT_noncoding <- ifelse(is.na(merged\$mRNA), FALSE, merged\$mRNA == "no")
    merged\$CPC2_noncoding <- !is.na(merged\$coding_probability)
    merged\$consensus_noncoding <- merged\$CPAT_noncoding & merged\$CPC2_noncoding

    write_tsv(merged, "coding_consensus.tsv")
    """
}
