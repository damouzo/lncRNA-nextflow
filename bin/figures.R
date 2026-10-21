#!/usr/bin/env Rscript
# F5.1: catalog figures. Non-blocking: the module uses errorStrategy 'ignore',
# and every panel is wrapped in tryCatch so a missing package or an odd table
# produces a manifest row instead of failing the run.
#
# Writes PNGs + figures_manifest.tsv (figure, path, description, status).
#
# gffcompare class codes (source: gffcompare manual) — only `x` is antisense:
#   c contained in a reference transcript | k contains a reference transcript
#   m retained introns, all matching | n retained introns, not all matching
#   j at least one junction shared | e single exon overlapping an exon and >=10 bp
#   of an intron (pre-mRNA fragment) | o other same-strand exonic overlap
#   s intron matching on the opposite strand (probable mapping error)
#   x exonic overlap on the opposite strand (antisense) | i entirely within an intron
#   y contains a reference transcript within its introns | p possible run-on
#   r repeat | u intergenic (lincRNA-like)
#
# usage:
#   figures.R --catalog catalog_genes.tsv --by-class tier_final_by_class.tsv --outdir .

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
catalog_path <- get_arg("--catalog")
by_class_path <- get_arg("--by-class")
outdir <- get_arg("--outdir", ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

manifest <- list()
add <- function(figure, path, description, status) {
  manifest[[length(manifest) + 1L]] <<- data.frame(
    figure = figure, path = path, description = description, status = status,
    stringsAsFactors = FALSE)
}
land <- function(figure, description, expr) {
  tryCatch({
    p <- expr
    f <- file.path(outdir, paste0(figure, ".png"))
    ggsave(f, p, width = 8, height = 6, dpi = 150)
    add(figure, basename(f), description, "ok")
    cat("figure:", figure, "ok\n")
  }, error = function(e) {
    add(figure, NA_character_, description, paste0("skipped: ", conditionMessage(e)))
    cat("figure:", figure, "skipped:", conditionMessage(e), "\n")
  })
}

catalog <- tryCatch(read_tsv(catalog_path, show_col_types = FALSE), error = function(e) NULL)
if (is.null(catalog) || nrow(catalog) == 0) {
  add("catalog", NA_character_, "catalog empty or unreadable", "skipped")
} else {
  land("fig_catalog_origin_class", "genes by origin and class_code", {
    d <- as.data.frame(table(origin = catalog$origin, class_code = catalog$class_code))
    ggplot(d, aes(x = origin, y = Freq, fill = class_code)) +
      geom_col() + theme_minimal() +
      labs(title = "Catalog composition", x = "origin", y = "genes")
  })

  if ("n_tier1_transcripts" %in% colnames(catalog)) {
    land("fig_tier1_transcripts", "Tier 1 transcripts per gene", {
      d <- catalog
      d$n_tier1_transcripts <- suppressWarnings(as.numeric(d$n_tier1_transcripts))
      d <- d[!is.na(d$n_tier1_transcripts), , drop = FALSE]
      ggplot(d, aes(x = n_tier1_transcripts)) +
        geom_histogram(binwidth = 1) + theme_minimal() +
        labs(title = "Tier 1 transcripts per gene", x = "n tier1 transcripts")
    })
  }
}

by_class <- tryCatch(read_tsv(by_class_path, show_col_types = FALSE), error = function(e) NULL)
land("fig_tier_by_class", "tier_final_by_class counts", {
  if (is.null(by_class) || nrow(by_class) == 0) stop("tier_final_by_class empty")
  d <- by_class
  d$n <- suppressWarnings(as.numeric(d$n))
  ggplot(d, aes(x = class_set, y = n, fill = unit)) +
    geom_col(position = "dodge") + theme_minimal() +
    labs(title = "Tier counts by class set", x = "class set", y = "n")
})

out <- if (length(manifest) > 0) do.call(rbind, manifest) else
  data.frame(figure = character(), path = character(), description = character(), status = character())
write_tsv(out, file.path(outdir, "figures_manifest.tsv"))
cat(sprintf("figures: %d panel(s), %d ok\n", nrow(out), sum(out$status == "ok")))
