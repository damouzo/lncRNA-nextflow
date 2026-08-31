// DESeq2 differential expression
// Flexible design formulas, batch-subsettable contrasts, apeglm shrinkage, lfcThreshold
//
// Contrasts CSV format (tab-separated):
//   contrast_name<TAB>numerator<TAB>denominator<TAB>batch
// If batch is non-empty, DESeq2 runs on the batch-restricted subset of samples.

process DESEQ2_DE {

    input:
    tuple val(group), path(counts_rds), path(coldata_csv), path(contrasts_tsv)
    val design_formula

    output:
    path "${group}_deseq2_results/*",              emit: all_results
    path "${group}_deseq2_results/dds.rds",        emit: dds_rds
    path "${group}_deseq2_results/design_diagnostics.txt", emit: design_diagnostics
    tuple val(group), path("${group}_deseq2_results/contrasts_summary.tsv"), emit: contrasts_summary

    script:
    def formula = design_formula.toString()
    """
    #!/usr/bin/env Rscript
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(apeglm)

    counts <- readRDS("${counts_rds}")
    coldata <- read.csv("${coldata_csv}", stringsAsFactors = FALSE)
    rownames(coldata) <- colnames(counts)

    # Detect if batch is constant — if so, simplify design to ~ condition
    formula_used <- "${formula}"
    if (grepl("batch", formula_used) && "batch" %in% colnames(coldata)) {
        batch_levels <- unique(coldata[["batch"]])
        batch_levels <- batch_levels[nchar(batch_levels) > 0]
        if (length(batch_levels) <= 1) {
            formula_used <- "~ condition"
            cat(sprintf("Batch has %d unique non-empty level(s) — using ~ condition\\n",
                        length(batch_levels)))
        }
    }

    # Pre-flight: check design matrix rank. A rank-deficient design (e.g. batch fully
    # confounded with condition) will silently produce unstable or uninterpretable
    # DESeq2 results. Warn early with a clear diagnostic.
    dir.create(paste0("${group}_deseq2_results"), showWarnings = FALSE, recursive = TRUE)
    diag_file <- file.path("${group}_deseq2_results", "design_diagnostics.txt")
    sink(diag_file)
    n_coef <- ncol(model.matrix(as.formula(formula_used), data = coldata))
    rank_mm <- qr(model.matrix(as.formula(formula_used), data = coldata))\$rank
    cat(sprintf("formula: %s\\n", formula_used))
    cat(sprintf("rank: %d / %d coefficients\\n", rank_mm, n_coef))
    if (rank_mm < n_coef) {
        cat("WARNING: design matrix is rank-deficient.\\n")
        cat("This usually means a covariate is confounded with condition or batch.\\n")
        if ("batch" %in% colnames(coldata) && nchar(paste(unique(coldata\$batch), collapse = "")) > 0) {
            cat("\\nCondition x batch contingency table:\\n")
            print(table(coldata\$condition, coldata\$batch))
        }
        cat("\\nThe pipeline will continue, but downstream results may be unreliable.\\n")
    } else {
        cat("OK\\n")
    }
    sink()
    # Also echo to stdout so it's visible in the Nextflow log
    writeLines(readLines(diag_file))

    contrasts <- read_tsv("${contrasts_tsv}", show_col_types = FALSE,
                          col_names = c("contrast_name", "numerator", "denominator", "batch"))

    dir.create("${group}_deseq2_results", showWarnings = FALSE)
    all_results <- list()

    for (i in seq_len(nrow(contrasts))) {
        ct_name  <- contrasts\$contrast_name[i]
        num      <- contrasts\$numerator[i]
        den      <- contrasts\$denominator[i]
        ct_batch <- contrasts\$batch[i]
        ct_batch <- if (is.na(ct_batch) || nchar(ct_batch) == 0) NULL else ct_batch

        cat(sprintf("\\n=== Contrast: %s (%s vs %s)", ct_name, num, den))
        if (!is.null(ct_batch)) {
            cat(sprintf(" [batch: %s]", ct_batch))
        }
        cat("\\n")

        # Subset to batch if specified
        if (!is.null(ct_batch)) {
            batch_col <- "batch"
            if (!(batch_col %in% colnames(coldata))) {
                cat(sprintf("  SKIP: 'batch' column not found in coldata\\n"))
                next
            }
            idx <- coldata[[batch_col]] == ct_batch
            if (sum(idx) < 4) {
                cat(sprintf("  SKIP: only %d samples in batch '%s'\\n", sum(idx), ct_batch))
                next
            }
            sub_counts <- counts[, idx, drop = FALSE]
            sub_coldata <- coldata[idx, , drop = FALSE]
            sub_coldata\$condition <- factor(sub_coldata\$condition)
            design <- as.formula("~ condition")
        } else {
            sub_counts <- counts
            sub_coldata <- coldata
            design <- as.formula(formula_used)
        }

        # Ensure condition is a factor with correct levels
        sub_coldata\$condition <- factor(sub_coldata\$condition,
                                          levels = unique(sub_coldata\$condition))

        # Check required levels exist
        if (!(num %in% levels(sub_coldata\$condition))) {
            cat(sprintf("  SKIP: numerator '%s' not in condition levels\\n", num))
            next
        }
        if (!(den %in% levels(sub_coldata\$condition))) {
            cat(sprintf("  SKIP: denominator '%s' not in condition levels\\n", den))
            next
        }

        # Relevel so 'den' is the reference: apeglm shrinkage below picks its
        # coefficient by name (condition_<num>_vs_<reference>), which only
        # matches the intended num-vs-den contrast when den is the reference
        # level. Without this, groups with >2 condition levels would silently
        # shrink toward the wrong comparison whenever the arbitrary default
        # reference (first level encountered) wasn't den.
        sub_coldata\$condition <- relevel(sub_coldata\$condition, ref = den)

        res <- tryCatch({
            dds <- DESeqDataSetFromMatrix(countData = sub_counts,
                                          colData = sub_coldata,
                                          design = design)
            # Drop unused factor levels
            dds\$condition <- droplevels(dds\$condition)
            dds <- DESeq(dds)
            results(dds, contrast = c("condition", num, den),
                    lfcThreshold = ${params.lfc_threshold},
                    alpha = ${params.fdr_threshold})
        }, error = function(e) {
            cat(sprintf("  ERROR: %s\\n", e\$message))
            NULL
        })

        if (is.null(res)) next

        # LFC shrinkage — pass res = res so the shrunk table keeps the same
        # test statistics/padj computed with the configured alpha and
        # lfcThreshold, instead of lfcShrink silently recomputing both from
        # scratch with DESeq2 defaults (alpha = 0.1, lfcThreshold = 0).
        res_shrunk <- tryCatch({
            lfcShrink(dds, res = res, coef = resultsNames(dds)[grep(
                paste0("condition_", make.names(num)), resultsNames(dds))],
                type = "apeglm", quiet = TRUE)
        }, error = function(e) {
            cat(sprintf("  WARNING: apeglm failed (%s), using ashr\\n", e\$message))
            lfcShrink(dds, res = res, coef = resultsNames(dds)[grep(
                paste0("condition_", make.names(num)), resultsNames(dds))],
                type = "ashr", quiet = TRUE)
        })

        res_df <- as.data.frame(res_shrunk)
        res_df\$contrast <- ct_name
        res_df\$gene_id <- rownames(res_df)
        res_df\$batch_subset <- if (is.null(ct_batch)) "all" else ct_batch

        fname <- file.path("${group}_deseq2_results", paste0(ct_name, "_DE.tsv"))
        write_tsv(res_df, fname)
        cat(sprintf("  Wrote %s (%d genes)\\n", fname, nrow(res_df)))

        all_results[[ct_name]] <- res_df

        # Save first valid dds
        if (i == 1 || !exists("saved_dds")) {
            saved_dds <- dds
        }
    }

    if (length(all_results) > 0) {
        summary <- do.call(rbind, all_results)
        write_tsv(as.data.frame(summary), "${group}_deseq2_results/contrasts_summary.tsv")
    } else {
        write_tsv(data.frame(), "${group}_deseq2_results/contrasts_summary.tsv")
    }

    if (exists("saved_dds")) {
        saveRDS(saved_dds, "${group}_deseq2_results/dds.rds")
    }
    """
}
