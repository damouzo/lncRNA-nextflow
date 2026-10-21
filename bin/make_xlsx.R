#!/usr/bin/env Rscript
# Slice D: TSV -> xlsx wrapper. The module uses errorStrategy 'ignore', so a missing
# openxlsx or a bad sheet never fails the run; the TSV stays the source of truth.
# Gene symbols and ids are written as text so Excel does not mangle them.

suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

inp <- get_arg("--in")
out <- get_arg("--out")
sheet <- get_arg("--sheet", "table")
if (is.null(inp) || is.null(out)) {
  stop("usage: make_xlsx.R --in table.tsv --out table.xlsx [--sheet name]")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("openxlsx is not installed in this image")
}

d <- as.data.frame(read_tsv(inp, show_col_types = FALSE))
d[] <- lapply(d, as.character)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheet)
openxlsx::writeData(wb, sheet, d)
openxlsx::saveWorkbook(wb, out, overwrite = TRUE)
cat(sprintf("wrote %s (%d rows)\n", out, nrow(d)))
