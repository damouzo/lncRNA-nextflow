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
#             [--tpm <group>_tpm.tsv ...] [--sample-metadata sample_index.tsv] \
#             [--comparisons comparisons.csv] --outdir .

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
comparisons_path <- get_arg("--comparisons")
outdir <- get_arg("--outdir", ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

manifest <- list()
add <- function(figure, path, description, status) {
  manifest[[length(manifest) + 1L]] <<- data.frame(
    figure = figure, path = path, description = description, status = status,
    stringsAsFactors = FALSE)
}
land <- function(figure, description, expr, width = 8, height = 6) {
  tryCatch({
    p <- expr
    f <- file.path(outdir, paste0(figure, ".png"))
    ggsave(f, p, width = width, height = height, dpi = 150)
    add(figure, basename(f), description, "ok")
    cat("figure:", figure, "ok\n")
  }, error = function(e) {
    add(figure, NA_character_, description, paste0("skipped: ", conditionMessage(e)))
    cat("figure:", figure, "skipped:", conditionMessage(e), "\n")
  })
}

# Canonical GRCh38 chromosomes in biological order. Anything else (patches,
# scaffolds, alt contigs) is dropped so the axis stays readable.
CANON_CHR <- c(as.character(1:22), "X", "Y", "MT")

catalog <- tryCatch(read_tsv(catalog_path, show_col_types = FALSE), error = function(e) NULL)
if (is.null(catalog) || nrow(catalog) == 0) {
  add("catalog", NA_character_, "catalog empty or unreadable", "skipped")
} else {
  # Catalog composition by gene_group (coding vs annotated lncRNA vs novel vs other).
  land("fig_catalog_gene_group", "genes by gene_group", {
    if (!"gene_group" %in% colnames(catalog)) stop("no gene_group column")
    d <- as.data.frame(table(gene_group = catalog$gene_group))
    ggplot(d, aes(x = gene_group, y = Freq)) +
      geom_col(fill = "#4c72b0") + theme_minimal() +
      labs(title = "Catalog genes by gene_group", x = "gene_group", y = "genes")
  })

  # What makes up `gene_group == other`: reference genes that are neither
  # protein-coding nor annotated lncRNAs (pseudogenes, TEC, snRNA, miRNA, ...).
  # Horizontal bars so a long biotype list stays readable.
  land("fig_catalog_gene_group_other", "other reference genes by biotype", {
    if (!all(c("gene_group", "gene_biotype") %in% colnames(catalog)))
      stop("no gene_group/gene_biotype column")
    d <- as.data.frame(table(biotype = catalog$gene_biotype[catalog$gene_group == "other"]))
    d <- d[d$Freq > 0, , drop = FALSE]
    if (nrow(d) == 0) stop("no 'other' genes")
    d$biotype <- factor(d$biotype, levels = d$biotype[order(d$Freq)])
    ggplot(d, aes(x = biotype, y = Freq)) +
      geom_col(fill = "#4c72b0") + coord_flip() + theme_minimal() +
      labs(title = "gene_group = other, by biotype", x = "gene_biotype", y = "genes")
  }, width = 9, height = 7)

  # Novel genes only: class_code on y, tier on x, count in each cell. A gene whose
  # transcripts carry several gffcompare codes appears once, under its full code
  # set, so the same gene is never counted twice across single-code rows.
  land("fig_catalog_novel_class_tier", "novel genes by class_code and tier (cell counts)", {
    need <- c("origin", "class_code", "tier")
    if (!all(need %in% colnames(catalog))) stop("missing origin/class_code/tier")
    d <- catalog[catalog$origin == "novel", , drop = FALSE]
    if (nrow(d) == 0) stop("no novel genes")
    d <- as.data.frame(table(class_code = d$class_code, tier = d$tier))
    d <- d[d$Freq > 0, , drop = FALSE]
    tot <- tapply(d$Freq, d$class_code, sum)
    d$class_code <- factor(d$class_code, levels = names(sort(tot, decreasing = TRUE)))
    ggplot(d, aes(x = tier, y = class_code, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), size = 2.8, color = "grey15") +
      scale_fill_gradient(low = "#e8f0fb", high = "#08306b") +
      theme_minimal() +
      labs(title = "Novel lncRNAs by class_code and tier",
           subtitle = "Cell = number of genes. Each gene is counted once, under its full gffcompare class_code set.",
           x = "tier", y = "class_code", fill = "genes")
  }, width = 7, height = 14)

  # Genes per canonical chromosome by origin: one panel with counts, a second
  # with the within-chromosome origin fraction (same figure, free y per panel).
  land("fig_catalog_chromosome", "genes per chromosome by origin (count and fraction)", {
    if (!"chromosome" %in% colnames(catalog)) stop("no chromosome column")
    d <- catalog[catalog$chromosome %in% CANON_CHR, , drop = FALSE]
    if (nrow(d) == 0) stop("no canonical chromosomes")
    d$chromosome <- factor(d$chromosome, levels = CANON_CHR)
    cnt <- as.data.frame(table(chromosome = d$chromosome, origin = d$origin))
    tot <- tapply(cnt$Freq, cnt$chromosome, sum)
    cnt$panel <- "genes"
    # Only chromosomes with genes get a fraction (avoids NaN rows for empty ones).
    frac <- cnt[as.numeric(tot[as.character(cnt$chromosome)]) > 0, , drop = FALSE]
    frac$Freq <- frac$Freq / as.numeric(tot[as.character(frac$chromosome)])
    frac$panel <- "origin_fraction"
    dd <- rbind(cnt, frac)
    dd$panel <- factor(dd$panel, levels = c("genes", "origin_fraction"))
    ggplot(dd, aes(x = chromosome, y = Freq, fill = origin)) +
      geom_col() +
      facet_wrap(~panel, ncol = 1, scales = "free_y") +
      scale_x_discrete(limits = CANON_CHR) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = "Catalog genes per chromosome", x = "chromosome",
           y = "genes / origin fraction")
  }, width = 8, height = 8)

  # F5.4: lncRNAs vs reference protein-coding genes (gene span). Only lncRNAs
  # carry a class; reference_other rows are NA and are dropped here.
  if ("lnc_position_class" %in% colnames(catalog)) {
    land("fig_lnc_position_class", "lncRNAs by position class", {
      keep <- c("antisense", "intronic", "overlapping", "no_coding_overlap")
      d <- catalog[catalog$lnc_position_class %in% keep, , drop = FALSE]
      if (nrow(d) == 0) stop("no classified lncRNAs")
      d <- as.data.frame(table(cls = d$lnc_position_class, origin = d$origin))
      ggplot(d, aes(x = cls, y = Freq, fill = origin)) +
        geom_col(position = "dodge") + theme_minimal() +
        labs(title = "lncRNA position class", x = "lnc_position_class", y = "genes")
    })

    land("fig_lnc_position_expr", "expressed lncRNA fraction by position class", {
      if (length(tpm_paths) == 0) stop("no TPM tables provided")
      keep <- c("antisense", "intronic", "overlapping", "no_coding_overlap")
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
        # Exclude the annotation columns the TPM tables now carry (gene_group,
        # catalog_id, origin, class_code, tier, lnc_position_class) so they are
        # never mistaken for sample columns.
        samples <- setdiff(colnames(tp), c("gene_id", "cls", "gene_group", "catalog_id",
                                           "origin", "class_code", "tier", "lnc_position_class"))
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
        d$norm_group <- if ("norm_group" %in% colnames(meta))
          meta$norm_group[match(d$sample, meta$sample)] else "all"
      } else {
        d$condition <- d$sample
        d$norm_group <- "all"
      }
      d$condition[is.na(d$condition)] <- "NA"
      d$norm_group[is.na(d$norm_group)] <- "unknown"

      # Control conditions = the denominators declared in the comparisons file;
      # fall back to a name heuristic (Scramble / Healthy / Control / WT) so the
      # panel still separates groups when the file is unavailable.
      ctrl_conds <- character(0)
      if (!is.null(comparisons_path) && file.exists(comparisons_path)) {
        cmp <- tryCatch(read_tsv(comparisons_path, show_col_types = FALSE),
                        error = function(e) NULL)
        if (!is.null(cmp) && "denominator" %in% colnames(cmp))
          ctrl_conds <- unique(cmp$denominator)
      }
      is_ctrl <- function(c)
        c %in% ctrl_conds | grepl("scramble|healthy|control|wildtype", c, ignore.case = TRUE)
      conds <- sort(unique(d$condition))
      ctrl <- sort(conds[vapply(conds, is_ctrl, logical(1))])
      case <- sort(setdiff(conds, ctrl))
      lev <- c(ctrl, case)
      d$condition <- factor(d$condition, levels = lev)
      pal <- c(
        if (length(ctrl)) colorRampPalette(c("#9ecae1", "#08306b"))(length(ctrl)),
        if (length(case)) colorRampPalette(c("#fcae91", "#67000d"))(length(case)))
      names(pal) <- lev
      ggplot(d, aes(x = cls, y = frac, fill = condition)) +
        geom_boxplot() +
        facet_wrap(~norm_group, ncol = 2) +
        scale_fill_manual(values = pal, drop = FALSE) +
        theme_minimal() +
        labs(title = "Expressed lncRNA fraction by position class",
             x = "lnc_position_class", y = "fraction TPM >= 1")
    }, width = 10, height = 7)
  }
}

by_class <- tryCatch(read_tsv(by_class_path, show_col_types = FALSE), error = function(e) NULL)
land("fig_tier_by_class", "Tier 1 composition by gffcompare class set", {
  if (is.null(by_class) || nrow(by_class) == 0) stop("tier_final_by_class empty")
  d <- by_class
  d$n <- suppressWarnings(as.numeric(d$n))
  # `tier1_allow_single_exon` changes the Tier 1 set; both settings are in the
  # table, so facet on it. Without the facet the two rows per (class_set, unit)
  # would overlap and only the last drawn setting would be visible.
  if (!"tier1_allow_single_exon" %in% colnames(d)) stop("no tier1_allow_single_exon column")
  ggplot(d, aes(x = class_set, y = n, fill = unit)) +
    geom_col(position = "dodge") + theme_minimal() +
    facet_wrap(~tier1_allow_single_exon, nrow = 1) +
    labs(title = "Tier 1 counts by class set",
         subtitle = "Composition of the Tier 1 catalog by gffcompare class set of its transcripts",
         x = "class set", y = "n", fill = "unit")
}, width = 12, height = 6)

out <- if (length(manifest) > 0) do.call(rbind, manifest) else
  data.frame(figure = character(), path = character(), description = character(), status = character())
write_tsv(out, file.path(outdir, "figures_manifest.tsv"))
cat(sprintf("figures: %d panel(s), %d ok\n", nrow(out), sum(out$status == "ok")))
