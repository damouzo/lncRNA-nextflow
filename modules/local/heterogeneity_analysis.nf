// B5.5: Expression-heterogeneity diagnostic (optional, opt-in per norm_group)
// Complementary to DESeq2 (B4), not a replacement — see Carrasco-Leon et al. 2021
// (Leukemia): a standard group-mean DE test under-detects lncRNAs because their
// expression is far more heterogeneous across patients than protein-coding genes.
// Runs per norm_group, looping over contrasts internally (same convention as
// DESEQ2_DE/REPORT_CONTRAST); case/control assignment and batch subsetting are
// copied verbatim from DESEQ2_DE's contrasts_tsv handling.
//
// Two analyses per contrast:
//   a. CV of log2(TPM+1) expression within the case group: lncRNA (novel + annotated)
//      vs protein-coding, compared with a Welch two-sample t-test. The
//      "*_heterogeneity_cv_summary.tsv" output is long format — one row per gene with
//      its own CV plus the contrast-level test result repeated on every row — so the
//      report can plot both raw CV distributions and read the t-test from one file.
//   b. Per-gene (lncRNA universe only), per-case-sample up/down/no_change call against
//      the control group's distribution, then a summary "altered" call.
//      Control distribution = mean +/- 1 SD of log2(TPM+1) in the control group
//      (documented choice — the paper's methods section doesn't specify one explicitly).
//      A gene is "altered" if concordant (same direction) in
//      >= params.heterogeneity_concordant_frac of case samples AND discordant (opposite
//      direction) in < params.heterogeneity_discordant_frac; otherwise "not_altered".

process HETEROGENEITY_ANALYSIS {

    input:
    tuple val(group), path(counts_rds), path(tpm_rds), path(coldata_csv), path(contrasts_tsv)
    path gene_catalog       // gene_catalog.tsv: origin/gene_biotype annotations

    output:
    path "${group}_heterogeneity_results/*",                                              emit: all_results
    tuple val(group), path("${group}_heterogeneity_results/heterogeneity_summary.tsv"),    emit: heterogeneity_summary
    tuple val(group), path("${group}_heterogeneity_results/heterogeneity_cv_summary.tsv"), emit: heterogeneity_cv_summary

    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)

    CONCORDANT_FRAC <- ${params.heterogeneity_concordant_frac}
    DISCORDANT_FRAC <- ${params.heterogeneity_discordant_frac}

    counts <- readRDS("${counts_rds}")
    tpm    <- readRDS("${tpm_rds}")
    coldata <- read.csv("${coldata_csv}", stringsAsFactors = FALSE)
    # Same convention as DESEQ2_DE: coldata rows are assumed pre-ordered to match the
    # expression matrix columns (both derive from the same tximport_filter.nf sample order).
    rownames(coldata) <- colnames(counts)

    # Restrict TPM to the DESeq2-tested gene universe (rows of counts_rds), same sample columns
    shared_genes <- intersect(rownames(counts), rownames(tpm))
    tpm <- tpm[shared_genes, colnames(counts), drop = FALSE]
    log2tpm <- log2(tpm + 1)

    catalog <- tryCatch(read_tsv("${gene_catalog}", show_col_types = FALSE), error = function(e) NULL)
    lnc_ids <- character()
    pc_ids  <- character()
    if (!is.null(catalog)) {
        lnc_ids <- catalog\$gene_id[catalog\$origin %in% c("novel", "annotated_lncrna")]
        pc_ids  <- catalog\$gene_id[catalog\$gene_biotype == "protein_coding"]
        if (length(pc_ids) == 0) {
            pc_ids <- catalog\$gene_id[catalog\$origin == "reference_other"]
            cat("WARNING: no protein_coding biotype in catalog, using all non-lncRNA genes for the CV comparison\\n")
        }
    } else {
        cat("WARNING: gene catalog unreadable, lncRNA/coding CV comparison will be empty\\n")
    }

    contrasts <- read_tsv("${contrasts_tsv}", show_col_types = FALSE,
                          col_names = c("contrast_name", "numerator", "denominator", "batch"))

    dir.create("${group}_heterogeneity_results", showWarnings = FALSE)

    safe_mean   <- function(x) if (length(x) == 0) NA_real_ else mean(x)
    safe_median <- function(x) if (length(x) == 0) NA_real_ else median(x)

    # Per-gene CV (sd/mean) on log2(TPM+1) across the given samples. Genes with a
    # non-positive or non-finite mean are dropped — CV is undefined/unstable there.
    gene_cv <- function(mat) {
        if (nrow(mat) == 0) return(numeric())
        m <- rowMeans(mat)
        s <- apply(mat, 1, sd)
        keep <- is.finite(m) & m > 0 & is.finite(s)
        s[keep] / m[keep]
    }

    gene_results <- list()
    cv_results   <- list()

    for (i in seq_len(nrow(contrasts))) {
        ct_name  <- contrasts\$contrast_name[i]
        num      <- contrasts\$numerator[i]
        den      <- contrasts\$denominator[i]
        ct_batch <- contrasts\$batch[i]
        ct_batch <- if (is.na(ct_batch) || nchar(ct_batch) == 0) NULL else ct_batch

        cat(sprintf("\\n=== Contrast: %s (%s vs %s)", ct_name, num, den))
        if (!is.null(ct_batch)) cat(sprintf(" [batch: %s]", ct_batch))
        cat("\\n")

        sub_coldata <- coldata
        if (!is.null(ct_batch)) {
            if (!("batch" %in% colnames(coldata))) {
                cat("  SKIP: 'batch' column not found in coldata\\n")
                next
            }
            idx <- coldata[["batch"]] == ct_batch
            if (sum(idx) < 4) {
                cat(sprintf("  SKIP: only %d samples in batch '%s'\\n", sum(idx), ct_batch))
                next
            }
            sub_coldata <- coldata[idx, , drop = FALSE]
        }

        if (!(num %in% sub_coldata\$condition) || !(den %in% sub_coldata\$condition)) {
            cat("  SKIP: numerator/denominator not present in this (sub)coldata\\n")
            next
        }

        case_samples    <- rownames(sub_coldata)[sub_coldata\$condition == num]
        control_samples <- rownames(sub_coldata)[sub_coldata\$condition == den]

        if (length(control_samples) < 2) {
            cat(sprintf("  SKIP: only %d control sample(s), SD undefined for classification\\n",
                        length(control_samples)))
            next
        }
        if (length(case_samples) < 1) {
            cat("  SKIP: no case samples\\n")
            next
        }

        case_expr    <- log2tpm[, case_samples, drop = FALSE]
        control_expr <- log2tpm[, control_samples, drop = FALSE]

        # ---- (a) CV comparison: lncRNA vs protein-coding, within the case group ----
        cv_lnc <- gene_cv(case_expr[rownames(case_expr) %in% lnc_ids, , drop = FALSE])
        cv_pc  <- gene_cv(case_expr[rownames(case_expr) %in% pc_ids,  , drop = FALSE])

        cv_test <- if (length(cv_lnc) >= 2 && length(cv_pc) >= 2) {
            tryCatch(t.test(cv_lnc, cv_pc), error = function(e) NULL)
        } else NULL

        # Long format: one row per gene with its own CV, plus the contrast-level
        # summary/t-test columns repeated on every row (denormalized on purpose, so the
        # report can both plot the two raw CV distributions and read the test result
        # from the same file without a second join).
        cv_long <- bind_rows(
            if (length(cv_lnc) > 0) {
                data.frame(gene_id = names(cv_lnc), origin_group = "lncrna", cv = unname(cv_lnc), stringsAsFactors = FALSE)
            } else data.frame(gene_id = character(), origin_group = character(), cv = numeric()),
            if (length(cv_pc) > 0) {
                data.frame(gene_id = names(cv_pc), origin_group = "coding", cv = unname(cv_pc), stringsAsFactors = FALSE)
            } else data.frame(gene_id = character(), origin_group = character(), cv = numeric())
        )
        cv_long\$contrast         <- ct_name
        cv_long\$n_lncrna_genes   <- length(cv_lnc)
        cv_long\$n_coding_genes   <- length(cv_pc)
        cv_long\$mean_cv_lncrna   <- safe_mean(cv_lnc)
        cv_long\$mean_cv_coding   <- safe_mean(cv_pc)
        cv_long\$median_cv_lncrna <- safe_median(cv_lnc)
        cv_long\$median_cv_coding <- safe_median(cv_pc)
        cv_long\$t_statistic      <- if (!is.null(cv_test)) unname(cv_test\$statistic) else NA_real_
        cv_long\$p_value          <- if (!is.null(cv_test)) cv_test\$p.value else NA_real_
        cv_results[[ct_name]] <- cv_long

        # ---- (b) per-gene, per-case-sample concordance vs the control distribution ----
        gene_ids_b <- intersect(rownames(case_expr), lnc_ids)
        if (length(gene_ids_b) == 0) {
            cat("  No lncRNA genes available for per-sample classification\\n")
        } else {
            ctrl_sub  <- control_expr[gene_ids_b, , drop = FALSE]
            ctrl_mean <- rowMeans(ctrl_sub)
            ctrl_sd   <- apply(ctrl_sub, 1, sd)

            case_sub <- case_expr[gene_ids_b, , drop = FALSE]
            n_case   <- ncol(case_sub)

            # sd == 0 (identical control replicates) means any non-zero difference
            # already counts as up/down — no extra special-casing needed here.
            up_mat <- case_sub > (ctrl_mean + ctrl_sd)
            dn_mat <- case_sub < (ctrl_mean - ctrl_sd)

            n_up <- rowSums(up_mat)
            n_dn <- rowSums(dn_mat)
            n_no <- n_case - n_up - n_dn

            frac_up <- n_up / n_case
            frac_dn <- n_dn / n_case
            concordant_frac <- pmax(frac_up, frac_dn)
            discordant_frac <- pmin(frac_up, frac_dn)
            direction <- ifelse(frac_up >= frac_dn, "up", "down")

            altered_call <- ifelse(
                concordant_frac >= CONCORDANT_FRAC & discordant_frac < DISCORDANT_FRAC,
                "altered", "not_altered"
            )

            gene_results[[ct_name]] <- data.frame(
                contrast         = ct_name,
                gene_id          = gene_ids_b,
                n_case           = n_case,
                n_control        = length(control_samples),
                n_up             = n_up,
                n_down           = n_dn,
                n_no_change      = n_no,
                concordant_frac  = concordant_frac,
                discordant_frac  = discordant_frac,
                direction        = direction,
                altered_call     = altered_call,
                stringsAsFactors = FALSE,
                row.names        = NULL
            )
        }

        # Per-contrast files, mirroring DESEQ2_DE's "\${ct_name}_DE.tsv" convention
        fname_gene <- file.path("${group}_heterogeneity_results", paste0(ct_name, "_heterogeneity.tsv"))
        fname_cv   <- file.path("${group}_heterogeneity_results", paste0(ct_name, "_heterogeneity_cv_summary.tsv"))
        write_tsv(if (!is.null(gene_results[[ct_name]])) gene_results[[ct_name]] else data.frame(), fname_gene)
        write_tsv(cv_results[[ct_name]], fname_cv)
        cat(sprintf("  Wrote %s, %s\\n", fname_gene, fname_cv))
    }

    gene_summary <- if (length(gene_results) > 0) bind_rows(gene_results) else data.frame()
    cv_summary   <- if (length(cv_results)   > 0) bind_rows(cv_results)   else data.frame()

    write_tsv(gene_summary, "${group}_heterogeneity_results/heterogeneity_summary.tsv")
    write_tsv(cv_summary,   "${group}_heterogeneity_results/heterogeneity_cv_summary.tsv")

    cat(sprintf("\\nDONE: %d contrast(s) with gene-level results, %d with CV summary\\n",
                length(gene_results), length(cv_results)))
    """
}
