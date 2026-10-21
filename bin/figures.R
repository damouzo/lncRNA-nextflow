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
#   figures.R --catalog catalog_genes.tsv --by-class tier_final_by_class.tsv \
#             [--tpm <group>_tpm.tsv ...] [--sample-metadata sample_metadata.tsv] --outdir .

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
get_args_all <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(character(0))
  vapply(idx, function(i) if (i < length(args)) args[i + 1] else NA_character_,
         character(1))
}
catalog_path <- get_arg("--catalog")
by_class_path <- get_arg("--by-class")
tpm_paths <- get_args_all("--tpm")
sample_meta_path <- get_arg("--sample-metadata")
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

  land("fig_catalog_chromosome", "genes per chromosome by origin", {
    d <- as.data.frame(table(chromosome = catalog$chromosome, origin = catalog$origin))
    ggplot(d, aes(x = chromosome, y = Freq, fill = origin)) +
      geom_col() + theme_minimal() +
      labs(title = "Catalog genes per chromosome", x = "chromosome", y = "genes")
  })

  # F5.4: lncRNAs vs reference protein-coding genes (gene span). Only lncRNAs
  # carry a class; reference_other rows are NA and are dropped here.
  if ("lnc_position_class" %in% colnames(catalog)) {
    land("fig_lnc_position_class", "lncRNAs by position class", {
      keep <- c("antisense", "intronic", "overlapping", "intergenic")
      d <- catalog[catalog$lnc_position_class %in% keep, , drop = FALSE]
      if (nrow(d) == 0) stop("no classified lncRNAs")
      d <- as.data.frame(table(cls = d$lnc_position_class, origin = d$origin))
      ggplot(d, aes(x = cls, y = Freq, fill = origin)) +
        geom_col(position = "dodge") + theme_minimal() +
        labs(title = "lncRNA position class", x = "lnc_position_class", y = "genes")
    })

    land("fig_lnc_position_expr", "expressed lncRNA fraction by position class", {
      if (length(tpm_paths) == 0) stop("no TPM tables provided")
      keep <- c("antisense", "intronic", "overlapping", "intergenic")
      cls_rows <- catalog[catalog$lnc_position_class %in% keep, , drop = FALSE]
      cmap <- setNames(cls_rows$lnc_position_class, cls_rows$gene_id)
      meta <- if (!is.null(sample_meta_path) && file.exists(sample_meta_path))
        read_tsv(sample_meta_path, show_col_types = FALSE) else NULL
      acc <- list()
      for (p in tpm_paths) {
        tp <- tryCatch(read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
        if (is.null(tp) || !"gene_id" %in% colnames(tp)) next
        tp <- tp[tp$gene_id %in% names(cmap), , drop = FALSE]
        if (nrow(tp) == 0) next
        tp$cls <- unname(cmap[tp$gene_id])
        samples <- setdiff(colnames(tp), c("gene_id", "cls"))
        for (s in samples) {
          v <- suppressWarnings(as.numeric(tp[[s]]))
          for (cl in unique(tp$cls)) {
            idx <- tp$cls == cl
            acc[[length(acc) + 1L]] <- data.frame(
              sample = s, cls = cl,
              frac = mean(v[idx] >= 1, na.rm = TRUE), stringsAsFactors = FALSE)
          }
        }
      }
      if (length(acc) == 0) stop("no lncRNAs of the catalog in the TPM tables")
      d <- do.call(rbind, acc)
      if (!is.null(meta) && all(c("sample", "condition") %in% colnames(meta))) {
        d$condition <- meta$condition[match(d$sample, meta$sample)]
      } else {
        d$condition <- d$sample
      }
      d$condition[is.na(d$condition)] <- "NA"
      ggplot(d, aes(x = cls, y = frac, fill = condition)) +
        geom_boxplot() + theme_minimal() +
        labs(title = "Expressed lncRNA fraction by position class",
             x = "lnc_position_class", y = "fraction TPM >= 1")
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
