// BUILD_TX2GENE — the single producer of tx2gene_detailed.tsv for the whole run.
//
// The quantification universe (quant_universe.gtf) holds the catalog transcripts
// AND the reference-derived novel isoforms (isoform_classes, e.g. `j`). This
// module builds the canonical transcript->gene map the DE universe is defined
// from:
//   * reference transcripts      -> their reference gene;
//   * novel catalog transcripts  -> their own MSTRG gene;
//   * novel isoform transcripts  -> the REFERENCE PARENT gene (ref_gene_id).
//
// Remapping the isoforms to the parent is what keeps the parent gene's counts
// whole: transcript-level EM still splits reads between the reference isoforms
// and the novel isoform, and tximport sums them back under one gene_id. Without
// it the parent is diluted by a second, near-identical MSTRG gene.
//
// No other process may construct tx2gene_detailed.tsv (enforced by
// bin/tests/check_wiring.py). Isoforms are reported in novel_isoforms.tsv, never
// as novel lncRNAs, and never registered in the catalog.

process BUILD_TX2GENE {

    label 'process_single'

    input:
    path quant_universe_gtf    // ANNOTATION_FREEZE quant_universe.gtf (catalog + isoforms)
    path reference_gtf

    output:
    path "tx2gene_detailed.tsv", emit: tx2gene_detailed
    path "novel_isoforms.tsv",   emit: novel_isoforms
    path "tx2gene_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({ library(readr); library(dplyr) })

    # Keep in sync with ANNOTATION_FREEZE / BUILD_GENE_CATALOG: strip Ensembl
    # versions but leave MSTRG ids verbatim.
    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }
    strip_version <- function(x) sub("\\\\.[0-9]+\$", "", x)

    ISO_CLASSES <- c(${params.isoform_classes.collect { cc -> "\"${cc}\"" } .join(", ")})

    parse_lines <- function(path) {
        lines <- readLines(path)
        parts <- strsplit(lines, "\\t", fixed = TRUE)
        ok    <- lengths(parts) >= 9L
        list(feature = vapply(parts[ok], function(p) p[[3L]], character(1)),
             chrom   = vapply(parts[ok], function(p) p[[1L]], character(1)),
             start   = vapply(parts[ok], function(p) p[[4L]], character(1)),
             end     = vapply(parts[ok], function(p) p[[5L]], character(1)),
             strand  = vapply(parts[ok], function(p) p[[7L]], character(1)),
             attrs   = vapply(parts[ok], function(p) p[[9L]], character(1)))
    }

    # ---- Reference: transcripts -> gene, and gene -> biotype ----
    # A full reference GTF is huge; never load it whole. Only transcript and gene
    # lines are needed (matches ANNOTATION_FREEZE / BUILD_GENE_CATALOG).
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${reference_gtf}"),
                 "> ref_tx_lines.txt"))
    system(paste("grep -P '\\tgene\\t'", shQuote("${reference_gtf}"),
                 "> ref_gene_lines.txt"))
    ref_tx_gtf   <- parse_lines("ref_tx_lines.txt")
    ref_gene_gtf <- parse_lines("ref_gene_lines.txt")

    ref_tx <- data.frame(
        transcript_id = strip_version(extract_attr(ref_tx_gtf\$attrs, "transcript_id")),
        gene_id       = strip_version(extract_attr(ref_tx_gtf\$attrs, "gene_id")),
        stringsAsFactors = FALSE)
    ref_tx <- ref_tx[!is.na(ref_tx\$transcript_id) & !is.na(ref_tx\$gene_id), ]
    ref_tx\$category <- "reference"
    ref_genes <- unique(ref_tx\$gene_id)
    ref_biotype <- setNames(
        extract_attr(ref_gene_gtf\$attrs, "gene_biotype"),
        strip_version(extract_attr(ref_gene_gtf\$attrs, "gene_id")))

    # ---- Novel transcripts from the quantification universe (small) ----
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${quant_universe_gtf}"),
                 "> nov_tx_lines.txt"))
    q <- parse_lines("nov_tx_lines.txt")
    if (length(q\$attrs) == 0) {
        stop("No novel transcripts in quant_universe.gtf — cannot build tx2gene")
    }
    nov <- data.frame(
        transcript_id = extract_attr(q\$attrs, "transcript_id"),
        mstrg_gene_id = extract_attr(q\$attrs, "gene_id"),
        class_code    = extract_attr(q\$attrs, "class_code"),
        ref_gene_id   = strip_version(extract_attr(q\$attrs, "ref_gene_id")),
        ref_gene_name = extract_attr(q\$attrs, "cmp_ref_gene"),
        chrom = q\$chrom, start = q\$start, end = q\$end, strand = q\$strand,
        stringsAsFactors = FALSE)
    nov <- nov[!is.na(nov\$transcript_id) & !is.na(nov\$mstrg_gene_id), ]
    if (nrow(nov) == 0) {
        stop("No novel transcripts in quant_universe.gtf — cannot build tx2gene")
    }

    is_iso <- !is.na(nov\$class_code) & nov\$class_code %in% ISO_CLASSES

    # ---- Isoform -> parent gene (the remap) ----
    iso <- nov[is_iso, ]
    if (any(!nzchar(iso\$ref_gene_id) | is.na(iso\$ref_gene_id))) {
        bad <- head(iso\$transcript_id[is.na(iso\$ref_gene_id) | !nzchar(iso\$ref_gene_id)], 5)
        stop(sprintf("isoform transcript(s) without a reference parent (ref_gene_id): %s",
                     paste(bad, collapse = ", ")))
    }
    missing_parent <- setdiff(unique(iso\$ref_gene_id), ref_genes)
    if (length(missing_parent) > 0) {
        stop(sprintf("%d isoform parent gene(s) not in the reference tx2gene (e.g. %s)",
                     length(missing_parent), paste(head(missing_parent, 5), collapse = ", ")))
    }
    iso_map <- data.frame(transcript_id = iso\$transcript_id,
                          gene_id = iso\$ref_gene_id,
                          category = "novel_isoform",
                          stringsAsFactors = FALSE)

    # ---- Novel catalog transcripts keep their own gene ----
    nvl <- nov[!is_iso, c("transcript_id", "mstrg_gene_id")]
    colnames(nvl) <- c("transcript_id", "gene_id")
    nvl\$category <- "novel_lncRNA"

    tx2gene <- unique(rbind(ref_tx, nvl, iso_map))

    # Hard invariant: no isoform transcript may resolve to an MSTRG gene_id.
    leftover <- tx2gene\$category == "novel_isoform" & grepl("^MSTRG", tx2gene\$gene_id)
    if (any(leftover)) {
        stop(sprintf("%d isoform transcript(s) still map to an MSTRG gene_id — remap failed (e.g. %s)",
                     sum(leftover), paste(head(tx2gene\$transcript_id[leftover], 5), collapse = ", ")))
    }
    write_tsv(tx2gene, "tx2gene_detailed.tsv")

    # ---- novel_isoforms.tsv: reported, never catalogued ----
    iso_out <- data.frame(
        transcript_id = iso\$transcript_id,
        gene_id       = iso\$ref_gene_id,
        class_code    = iso\$class_code,
        host_gene_id  = iso\$ref_gene_id,
        host_gene_name = iso\$ref_gene_name,
        host_biotype  = unname(ref_biotype[iso\$ref_gene_id]),
        mstrg_gene_id = iso\$mstrg_gene_id,
        chrom = iso\$chrom, start = iso\$start, end = iso\$end, strand = iso\$strand,
        stringsAsFactors = FALSE)
    iso_out <- iso_out[order(iso_out\$host_gene_id, iso_out\$transcript_id), ]
    write_tsv(iso_out, "novel_isoforms.tsv")

    cat(sprintf("tx2gene: %d reference + %d novel catalog + %d remapped isoform(s) = %d rows\\n",
                nrow(ref_tx), nrow(nvl), nrow(iso_map), nrow(tx2gene)))

    sink("tx2gene_manifest.json")
    cat(sprintf('{
  "n_reference_transcripts": %d,
  "n_novel_catalog_transcripts": %d,
  "n_isoform_transcripts": %d,
  "n_total_transcripts": %d,
  "n_isoform_parent_genes": %d
}
', nrow(ref_tx), nrow(nvl), nrow(iso_map), nrow(tx2gene), length(unique(iso_map\$gene_id))))
    sink()

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
