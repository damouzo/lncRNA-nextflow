// A1b: Deterministic BED12 for RSeQC infer_experiment.py
//
// Built once from the reference GTF instead of per-sample: genes whose span does
// not overlap any other gene (any strand) and that are protein_coding, one
// representative (longest) transcript per gene. Non-overlapping CDS-like regions
// keep the `failed` fraction of infer_experiment meaningful, which the whole-GTF
// BED12 (overlapping / antisense exons) inflates.
//
// Deterministic: output is sorted by locus and one transcript per gene is chosen
// by (longest exon sum, then transcript_id), so the file is byte-stable across runs.

process BUILD_INFER_BED {

    input:
    path reference_gtf

    output:
    path "ref_exon.bed12", emit: bed12

    script:
    """
    #!/usr/bin/env Rscript
    suppressWarnings({
        library(rtracklayer)
        library(GenomicRanges)
    })

    gtf <- import("${reference_gtf}")
    gtf <- gtf[width(gtf) > 0]
    seqlevels(gtf) <- seqlevelsInUse(gtf)

    feats <- as.character(mcols(gtf)\$type)

    g <- gtf[feats == "gene"]
    if (length(g) == 0) stop("BUILD_INFER_BED: no 'gene' features in reference GTF")

    genes <- data.frame(
        gene_id  = as.character(mcols(g)\$gene_id),
        biotype  = as.character(mcols(g)\$gene_biotype),
        stringsAsFactors = FALSE
    )

    # Drop genes that overlap any other gene on either strand. findOverlaps on the
    # GRanges respects seqnames (an IRanges-only call would merge chromosomes) and,
    # with ignore.strand = TRUE, also antisense overlaps — required by F0.8.
    hits <- findOverlaps(g, g, type = "any", ignore.strand = TRUE)
    self <- queryHits(hits) == subjectHits(hits)
    ov   <- hits[!self]
    overlapped <- unique(c(queryHits(ov), subjectHits(ov)))

    coding_rows <- which(genes\$biotype == "protein_coding")
    keep_rows   <- setdiff(coding_rows, overlapped)
    keep_genes  <- genes\$gene_id[keep_rows]
    cat(sprintf("Protein-coding genes: %d, non-overlapping: %d\\n",
                length(coding_rows), length(keep_rows)))

    ex  <- gtf[feats == "exon"]
    if (length(ex) == 0) stop("BUILD_INFER_BED: no 'exon' features in reference GTF")

    ex_id  <- as.character(mcols(ex)\$transcript_id)
    ex_gen <- as.character(mcols(ex)\$gene_id)
    ex_ok  <- !is.na(ex_id) & ex_gen %in% keep_genes
    ex     <- ex[ex_ok]
    ex_id  <- ex_id[ex_ok]
    ex_gen <- ex_gen[ex_ok]

    # Representative transcript per gene: longest summed exon length, tie by id.
    tx_len <- tapply(width(ex), ex_id, sum)
    tx_ids <- names(tx_len)
    tx_gen <- ex_gen[match(tx_ids, ex_id)]
    rep <- data.frame(gene = tx_gen, tx = tx_ids, len = as.integer(tx_len),
                      stringsAsFactors = FALSE)
    rep <- rep[!is.na(rep\$gene), ]
    rep <- rep[order(rep\$gene, -rep\$len, rep\$tx), ]
    rep <- rep[!duplicated(rep\$gene), ]
    if (nrow(rep) == 0) stop("BUILD_INFER_BED: no non-overlapping protein_coding genes to emit")

    ex_by_tx <- split(seq_along(ex_id), ex_id)

    rows <- vector("list", nrow(rep))
    for (i in seq_len(nrow(rep))) {
        idx <- ex_by_tx[[rep\$tx[i]]]
        if (is.null(idx) || length(idx) == 0) next
        e <- ex[idx]
        o <- order(start(e), end(e))
        e <- e[o]
        st <- start(e); en <- end(e)
        cs <- min(st) - 1L
        ce <- max(en)
        rows[[i]] <- data.frame(
            seqnames   = as.character(seqnames(e)[1]),
            start      = cs,
            end        = ce,
            name       = rep\$gene[i],
            score      = 0,
            strand     = as.character(strand(e)[1]),
            thickStart = cs,
            thickEnd   = ce,
            itemRgb    = "0",
            blockCount = length(idx),
            blockSizes = paste(en - st + 1L, collapse = ","),
            blockStarts = paste(st - 1L - cs, collapse = ","),
            stringsAsFactors = FALSE
        )
    }
    bed <- do.call(rbind, rows)
    bed <- bed[order(bed\$seqnames, bed\$start, bed\$name), ]

    write.table(bed, "ref_exon.bed12", sep = "\\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    cat(sprintf("Wrote %d BED12 records\\n", nrow(bed)))
    """
}
