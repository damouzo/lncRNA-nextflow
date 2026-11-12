#!/usr/bin/env Rscript
# BUILD_TX2GENE — single producer of tx2gene_detailed.tsv.
#
# The quantification universe (quant_universe.gtf) holds the catalog transcripts
# AND the reference-derived novel isoforms (isoform_classes, e.g. `j`):
#   * reference transcripts      -> their reference gene;
#   * novel catalog transcripts  -> their own MSTRG gene;
#   * novel isoform transcripts  -> the REFERENCE PARENT gene (ref_gene_id).
#
# Remapping the isoforms to the parent keeps the parent gene's counts whole:
# transcript-level EM splits reads between the reference isoforms and the novel
# isoform, and tximport sums them back under one gene_id.
#
# usage:
#   build_tx2gene.R --reference ref.gtf --universe quant_universe.gtf \
#       --isoform-classes j --outdir .

suppressPackageStartupMessages({ library(readr); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

reference_gtf     <- get_arg("--reference")
quant_universe_gtf <- get_arg("--universe")
isoform_classes   <- get_arg("--isoform-classes", "j")
outdir            <- get_arg("--outdir", ".")
if (is.null(reference_gtf) || is.null(quant_universe_gtf)) {
  stop("usage: build_tx2gene.R --reference ref.gtf --universe quant_universe.gtf ",
       "--isoform-classes j --outdir .")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Keep in sync with ANNOTATION_FREEZE / BUILD_GENE_CATALOG: strip Ensembl
# versions but leave MSTRG ids verbatim.
extract_attr <- function(attrs, key) {
    pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
    m <- regexpr(pat, attrs, perl = TRUE)
    ans <- rep(NA_character_, length(attrs))
    ok <- m > 0L
    if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*$'), '\\1', attrs[ok], perl = TRUE)
    ans
}
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

ISO_CLASSES <- trimws(strsplit(isoform_classes, ",", fixed = TRUE)[[1]])
ISO_CLASSES <- ISO_CLASSES[nzchar(ISO_CLASSES)]

parse_lines <- function(path) {
    lines <- readLines(path)
    parts <- strsplit(lines, "\t", fixed = TRUE)
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
tmp <- tempfile("tx2gene_")
ref_tx_lines   <- paste0(tmp, "_ref_tx.txt")
ref_gene_lines <- paste0(tmp, "_ref_gene.txt")
system(paste("grep -P '\ttranscript\t'", shQuote(reference_gtf), ">", shQuote(ref_tx_lines)))
system(paste("grep -P '\tgene\t'", shQuote(reference_gtf), ">", shQuote(ref_gene_lines)))
ref_tx_gtf   <- parse_lines(ref_tx_lines)
ref_gene_gtf <- parse_lines(ref_gene_lines)

ref_tx <- data.frame(
    transcript_id = strip_version(extract_attr(ref_tx_gtf$attrs, "transcript_id")),
    gene_id       = strip_version(extract_attr(ref_tx_gtf$attrs, "gene_id")),
    stringsAsFactors = FALSE)
ref_tx <- ref_tx[!is.na(ref_tx$transcript_id) & !is.na(ref_tx$gene_id), ]
ref_tx$category <- "reference"
ref_genes <- unique(ref_tx$gene_id)
ref_biotype <- setNames(
    extract_attr(ref_gene_gtf$attrs, "gene_biotype"),
    strip_version(extract_attr(ref_gene_gtf$attrs, "gene_id")))

# ---- Novel transcripts from the quantification universe (small) ----
nov_tx_lines <- paste0(tmp, "_nov_tx.txt")
system(paste("grep -P '\ttranscript\t'", shQuote(quant_universe_gtf), ">", shQuote(nov_tx_lines)))
q <- parse_lines(nov_tx_lines)
if (length(q$attrs) == 0) {
    stop("No novel transcripts in quant_universe.gtf — cannot build tx2gene")
}
nov <- data.frame(
    transcript_id = extract_attr(q$attrs, "transcript_id"),
    mstrg_gene_id = extract_attr(q$attrs, "gene_id"),
    class_code    = extract_attr(q$attrs, "class_code"),
    ref_gene_id   = strip_version(extract_attr(q$attrs, "ref_gene_id")),
    ref_gene_name = extract_attr(q$attrs, "cmp_ref_gene"),
    chrom = q$chrom, start = q$start, end = q$end, strand = q$strand,
    stringsAsFactors = FALSE)
nov <- nov[!is.na(nov$transcript_id) & !is.na(nov$mstrg_gene_id), ]
if (nrow(nov) == 0) {
    stop("No novel transcripts in quant_universe.gtf — cannot build tx2gene")
}

is_iso <- !is.na(nov$class_code) & nov$class_code %in% ISO_CLASSES

# ---- Isoform -> parent gene (the remap) ----
iso <- nov[is_iso, ]
if (any(!nzchar(iso$ref_gene_id) | is.na(iso$ref_gene_id))) {
    bad <- head(iso$transcript_id[is.na(iso$ref_gene_id) | !nzchar(iso$ref_gene_id)], 5)
    stop(sprintf("isoform transcript(s) without a reference parent (ref_gene_id): %s",
                 paste(bad, collapse = ", ")))
}
missing_parent <- setdiff(unique(iso$ref_gene_id), ref_genes)
if (length(missing_parent) > 0) {
    stop(sprintf("%d isoform parent gene(s) not in the reference tx2gene (e.g. %s)",
                 length(missing_parent), paste(head(missing_parent, 5), collapse = ", ")))
}
iso_map <- data.frame(transcript_id = iso$transcript_id,
                      gene_id = iso$ref_gene_id,
                      category = "novel_isoform",
                      stringsAsFactors = FALSE)

# ---- Novel catalog transcripts keep their own gene ----
nvl <- nov[!is_iso, c("transcript_id", "mstrg_gene_id")]
colnames(nvl) <- c("transcript_id", "gene_id")
nvl$category <- "novel_lncRNA"

tx2gene <- unique(rbind(ref_tx, nvl, iso_map))

# Hard invariant: no isoform transcript may resolve to an MSTRG gene_id.
leftover <- tx2gene$category == "novel_isoform" & grepl("^MSTRG", tx2gene$gene_id)
if (any(leftover)) {
    stop(sprintf("%d isoform transcript(s) still map to an MSTRG gene_id — remap failed (e.g. %s)",
                 sum(leftover), paste(head(tx2gene$transcript_id[leftover], 5), collapse = ", ")))
}
write_tsv(tx2gene, file.path(outdir, "tx2gene_detailed.tsv"))

# ---- novel_isoforms.tsv: reported, never catalogued ----
iso_out <- data.frame(
    transcript_id = iso$transcript_id,
    gene_id       = iso$ref_gene_id,
    class_code    = iso$class_code,
    host_gene_id  = iso$ref_gene_id,
    host_gene_name = iso$ref_gene_name,
    host_biotype  = unname(ref_biotype[iso$ref_gene_id]),
    mstrg_gene_id = iso$mstrg_gene_id,
    chrom = iso$chrom, start = iso$start, end = iso$end, strand = iso$strand,
    stringsAsFactors = FALSE)
iso_out <- iso_out[order(iso_out$host_gene_id, iso_out$transcript_id), ]
write_tsv(iso_out, file.path(outdir, "novel_isoforms.tsv"))

cat(sprintf("tx2gene: %d reference + %d novel catalog + %d remapped isoform(s) = %d rows\n",
            nrow(ref_tx), nrow(nvl), nrow(iso_map), nrow(tx2gene)))

sink(file.path(outdir, "tx2gene_manifest.json"))
cat(sprintf('{
  "n_reference_transcripts": %d,
  "n_novel_catalog_transcripts": %d,
  "n_isoform_transcripts": %d,
  "n_total_transcripts": %d,
  "n_isoform_parent_genes": %d
}
', nrow(ref_tx), nrow(nvl), nrow(iso_map), nrow(tx2gene), length(unique(iso_map$gene_id))))
sink()
