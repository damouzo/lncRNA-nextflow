## Feature: Heterogeneity Module

### Status: IN PROGRESS

New, optional, opt-in-per-norm_group feature: per-contrast expression-heterogeneity
diagnostic, based on Carrasco-Leon et al. 2021 (Leukemia) — lncRNA CV comparison
(lncRNA vs coding) + per-sample concordance-based "altered gene" calling against the
control distribution. Complementary to DESeq2, not a replacement. Already anticipated
in [pipeline_plan.md](pipeline_plan.md) item 8 ("Optional heterogeneity/recurrence
diagnostic, explicitly secondary to DESeq2 (§B5)").

Hard constraints respected: no `nextflow run`, no git commits, all comments/logs in
English. Not mixed with the F1-F4 audit findings above — this is purely new feature work.

### Read before writing code
[modules/local/deseq2_de.nf](modules/local/deseq2_de.nf),
[workflows/quantification.nf](workflows/quantification.nf),
[conf/modules.config](conf/modules.config),
[assets/report_template.qmd](assets/report_template.qmd),
[modules/local/report_contrast.nf](modules/local/report_contrast.nf),
[modules/local/tximport_filter.nf](modules/local/tximport_filter.nf),
[modules/local/cis_associations.nf](modules/local/cis_associations.nf),
[modules/local/functional_enrichment.nf](modules/local/functional_enrichment.nf),
[modules/local/build_gene_catalog.nf](modules/local/build_gene_catalog.nf),
[modules/local/split_normgroup.nf](modules/local/split_normgroup.nf),
[main.nf](main.nf), [nextflow.config](nextflow.config),
[nextflow_schema.json](nextflow_schema.json), samplesheet/comparisons templates and the
real coldata/comparisons.csv under the CRUK_DDX41_DHX34 test dataset (read-only).

### Key findings that shape the design
- **Granularity clarification:** DESEQ2_DE and REPORT_CONTRAST are Nextflow processes
  invoked **once per norm_group**, not once per contrast — they loop over all rows of
  the group's `contrasts_tsv` *inside the R script* and write one file per contrast
  (`${ct_name}_DE.tsv`) plus one stacked `contrasts_summary.tsv` (with a `contrast`
  column) that downstream steps (CIS_ASSOCIATIONS, FUNCTIONAL_ENRICHMENT, REPORT_CONTRAST)
  consume, keyed by `group` only. "Per-contrast granularity" in this codebase means "one
  row/file per contrast inside a per-group process", not "one Nextflow task per contrast".
  The new module follows the exact same convention (one process call per group, internal
  loop over contrasts) — no restructuring of the process graph needed.
- **Case/control assignment logic already exists and is reused as-is**, copied verbatim
  from `deseq2_de.nf`: read `contrasts_tsv` (`contrast_name, numerator, denominator,
  batch`), subset `coldata` to the batch if given, skip if <4 samples in batch subset,
  skip if numerator/denominator not present as condition levels. Coldata row order is
  assumed to already match `counts`/`tpm` column order (`rownames(coldata) <-
  colnames(counts)`), same pre-existing implicit assumption as `deseq2_de.nf` — not
  re-audited here, out of scope for this feature.
- **`tpm_rds` (from `TXIMPORT_AND_FILTER`) has zero downstream consumers today** — grepped
  the whole repo to confirm. Its process output is declared as a bare `path`, not a
  `tuple val(group), path(...)`, so it cannot safely be `.join()`-ed by group for >1
  norm_group (same class of hazard as F1). Since nothing currently depends on it, keying
  it by group is a strictly additive, zero-blast-radius change and is required to
  correctly compute CV on **normalized** (TPM) expression rather than raw filtered counts
  (raw counts would confound between-sample CV with library-size differences — especially
  relevant for heterogeneous primary-sample groups, exactly where this feature is meant to
  be turned on). This is the one pre-existing file touched outside the new module; it does
  not touch DESEQ2_DE / CIS_ASSOCIATIONS / FUNCTIONAL_ENRICHMENT / REPORT_CONTRAST wiring.
- Gene universe for the per-gene classification step (b): restricted to `origin %in%
  c("novel", "annotated_lncrna")` from `gene_catalog.tsv` — same origin categories
  `cis_associations.nf` already uses for "lncRNA side" — since that's the paper's actual
  focus and keeps runtime/output bounded. CV comparison (step a) contrasts this lncRNA set
  against `gene_biotype == "protein_coding"` (same fallback to `reference_other` as
  `cis_associations.nf` if `protein_coding` biotype is absent from the catalog).
  Both steps restrict to genes present in `counts_rds` (the DESeq2-tested universe), so the
  report's "altered but not already significant in DESeq2" comparison is apples-to-apples.

### Design decisions requiring documentation (paper is not fully explicit)
- **Control distribution rule:** mean ± 1 SD of the control group's `log2(TPM+1)` per
  gene, computed on the same batch-subsetted sample set DESeq2 would use for that
  contrast. A case sample is "up" if `expr > control_mean + control_sd`, "down" if
  `expr < control_mean - control_sd`, else "no_change". If the control subset has <2
  samples (SD undefined), the whole contrast is skipped for classification (documented
  in module comments), consistent with DESeq2's own `>=4 samples in batch subset` guard
  being a similar practical floor.
- **CV statistic:** per-gene CV computed on `log2(TPM+1)` across case-group samples only
  (`sd/mean` of the log2 values, matching the brief's "CV of expression (log2 scale)"
  instruction literally); genes with `mean(log2TPM) <= 0` in the case group are excluded
  from the CV comparison (count logged) since CV is undefined/unstable there. Welch two-
  sample t-test compares the per-gene CV of the lncRNA set vs the protein-coding set.
- **Altered call:** `concordant_frac` = fraction of case samples matching the gene's
  majority direction (up or down, whichever is larger); `discordant_frac` = fraction in
  the opposite direction; denominator is total case-group `n` for that contrast (samples
  classified "no_change" count toward neither numerator but do count toward the
  denominator). "altered" iff `concordant_frac >= params.heterogeneity_concordant_frac`
  (default 0.50) AND `discordant_frac < params.heterogeneity_discordant_frac` (default
  0.25); else "not_altered".

### Toggle mechanism
- `params.heterogeneity_analysis` = per-`norm_group` boolean map in `nextflow.config`,
  default `false` for any group not listed (fail-safe, opt-in). NOT in comparisons.csv/
  samplesheet.csv, matching existing convention that pipeline-behavior switches live in
  `nextflow.config`/`nextflow_schema.json` (like `cpat_coding_threshold`).
- `params.heterogeneity_concordant_frac` (default 0.50) / `params.heterogeneity_discordant_frac`
  (default 0.25): pipeline-wide, documented in `nextflow_schema.json` under a new
  `phase_b_thresholds` entry (existing section, since these are Phase B statistics
  thresholds like `fdr_threshold`/`lfc_threshold`).
- Wiring: `HETEROGENEITY_ANALYSIS` is only ever scheduled for norm_groups where the map
  says `true` — implemented via **channel `.join()` filtering** (join the per-group input
  channel against the enabled-groups channel; Nextflow `.join()` drops non-matching keys
  by default), not a process-level `when:` — this is more robust than `when:` (task is
  never submitted at all, not submitted-then-skipped) and there is no existing `when:`
  convention in this codebase to follow instead (grepped `workflows/**` for `when:` /
  `.filter` / `.branch` — none found).
- Disabled groups still need a channel entry so the `REPORT_CONTRAST` join (keyed by
  group) doesn't silently drop them. Solved with `Channel....collectFile()` generating a
  static, header-only placeholder TSV (no process needed — `collectFile` can synthesize a
  file from a string) that's `.combine()`-broadcast to every disabled group, mirroring the
  existing empty-`data.frame()`-tsv convention already used for "no results" cases in
  `cis_associations.nf`. `REPORT_CONTRAST` treats a 0-row heterogeneity table as "not
  present" and skips the new report section entirely — so disabled-group / non-existent-
  contrast reports are byte-for-byte identical to the pre-feature pipeline.

### Files touched (planned)
- NEW `modules/local/heterogeneity_analysis.nf`
- `modules/local/tximport_filter.nf` — key `tpm_rds` emit by group (additive only)
- `nextflow.config` — new params + container routing (added to the existing `r_bioc`
  container group, same as DESEQ2_DE/CIS_ASSOCIATIONS)
- `nextflow_schema.json` — document the 2 pipeline-wide fraction params + the per-group map
- `conf/modules.config` — publishDir/label for `HETEROGENEITY_ANALYSIS`
- `workflows/quantification.nf` — new channel wiring (additive; see below for the "when
  disabled, identical to before" verification)
- `modules/local/report_contrast.nf` — 2 new optional path inputs + `-P` params passed to
  quarto (always non-null; placeholder file when disabled/absent)
- `assets/report_template.qmd` — new conditional "Heterogeneity Analysis" section

### Verification of "no change when disabled" (requirement 3)
With `params.heterogeneity_analysis` empty/all-false (the default), `ch_het_flagged`
filters every group into `ch_het_disabled_groups`; `HETEROGENEITY_ANALYSIS` channel input
(`ch_het_in`, built via `.join()` against the enabled-groups channel) receives zero items
for every group, so the process is **never invoked** — confirmed by inspection of the
`.join()` semantics (Nextflow docs: default `remainder: false` drops unmatched keys from
both sides, i.e. an inner join). `ch_rpt`'s new two `.join()`s against
`ch_het_summary_all`/`ch_het_cv_summary_all` therefore always resolve to the placeholder
files for every group in the all-disabled case — identical arity/shape to before, and
`REPORT_CONTRAST`'s script only changes behavior (renders the new section) when it reads
a non-empty heterogeneity TSV, which never happens in the all-disabled case. DESEQ2_DE /
CIS_ASSOCIATIONS / FUNCTIONAL_ENRICHMENT channels are entirely untouched (no edits to
those modules or their call sites).

### Progress log
- [x] `heterogeneity_analysis.nf` written — per-group process, internal per-contrast loop,
  same skip/batch-subsetting rules as `deseq2_de.nf`. Writes per-contrast
  `${ct_name}_heterogeneity.tsv` (per-gene concordance/altered call, lncRNA universe only)
  and `${ct_name}_heterogeneity_cv_summary.tsv` (long format: one row per gene with its
  own CV plus the contrast's Welch t-test result repeated on every row — needed so the
  report can plot the raw CV distributions, not just the aggregate stats), plus group-level
  stacked `heterogeneity_summary.tsv` / `heterogeneity_cv_summary.tsv` (same `contrasts_summary.tsv`
  pattern as `deseq2_de.nf`).
- [x] `tximport_filter.nf` — `tpm_rds` emit changed to `tuple val(group), path(...)`.
  Confirmed via repo-wide grep that this output had zero downstream consumers before this
  change, so it's additive-only; does not touch DESEQ2_DE/CIS_ASSOCIATIONS/FUNCTIONAL_ENRICHMENT.
- [x] `nextflow.config`: added `params.heterogeneity_analysis` (per-group map, default
  `[:]` → false for unlisted groups), `params.heterogeneity_concordant_frac` (0.50),
  `params.heterogeneity_discordant_frac` (0.25); added `HETEROGENEITY_ANALYSIS` to the
  existing `r_bioc` container-routing regex (same container as DESEQ2_DE/CIS_ASSOCIATIONS).
  `nextflow_schema.json`: documented all three under the existing `phase_b_thresholds`
  definition.
  `conf/modules.config`: added a `withName: HETEROGENEITY_ANALYSIS` block
  (`quantification/heterogeneity` publishDir, `process_medium`/16GB, matching
  CIS_ASSOCIATIONS/FUNCTIONAL_ENRICHMENT sizing).
- [x] `workflows/quantification.nf` wiring added: `ch_het_enabled_groups` /
  `ch_het_disabled_groups` derived from `params.heterogeneity_analysis` via `.filter()`;
  `HETEROGENEITY_ANALYSIS`'s input channel is built by `.join()`-ing the per-group
  TXIMPORT outputs against `ch_het_enabled_groups`, so it is never scheduled for a
  disabled group. Disabled groups get a header-only placeholder TSV pair (built with
  `Channel.of(...).collectFile()` — no process needed) broadcast via `.combine()`, mixed
  with the real outputs so every group has exactly one entry before the `REPORT_CONTRAST`
  join. DESEQ2_DE/CIS_ASSOCIATIONS/FUNCTIONAL_ENRICHMENT channel wiring is untouched.
- [x] `report_contrast.nf` + `report_template.qmd` updated: `REPORT_CONTRAST` gained 2 new
  path inputs (`heterogeneity`, `heterogeneity_cv`) and passes them plus the 2 fraction
  params to quarto via `-P`, always (never conditionally — the file is just empty for
  disabled groups). The qmd adds a `## Heterogeneity Analysis` section, conditionally
  emitted via an inline-R heading (`` `r if (show_heterogeneity) "## Heterogeneity Analysis"` ``)
  and `eval=show_heterogeneity` chunk options — genuinely absent from the rendered HTML
  when there's no data, not just an empty/placeholder section, per requirement 4. Content:
  a violin+boxplot of the two raw CV distributions (lncRNA vs coding) with the Welch
  t-test result printed underneath, and a table of genes with `altered_call == "altered"`
  whose `gene_id` is NOT in the current contrast's DESeq2 `padj < 0.05` set (joined against
  `gene_catalog` for `gene_name`/`class_code`).
- [x] Final read-through done: `get_errors` run on all touched files. Two pre-existing,
  unrelated lint findings in `nextflow.config` (`ghcr_base` config-in-process-block style,
  `memory = '48 GB'` string-vs-MemoryUnit) predate this feature and were not introduced by
  it — left alone, out of scope. The "Parameter was not used" / "`Channel` deprecated"
  notices on the new `quantification.nf` code match the file's own pre-existing style
  (unprefixed destructured closure params, capital-`Channel` factories used throughout
  `main.nf`/this file already) — consistent with the codebase, not new problems.

### Verified NOT touched
`modules/local/deseq2_de.nf`, `modules/local/cis_associations.nf`,
`modules/local/functional_enrichment.nf` — no edits. Their channel wiring in
`quantification.nf` (`ch_deseq`, `ch_cis`, `ch_enrich`, and all their `.join()`s) is
byte-identical to before this feature. The only pre-existing files touched are
`tximport_filter.nf` (additive `tpm_rds` keying, justified above) and `report_contrast.nf`
(2 new optional inputs, discussed above and explicitly permitted by the brief as long as
the disabled-group behavior is unchanged).

### Status: DONE (implementation complete; not executed — hard constraint: no `nextflow run`)
Recommend a human/next-agent do a `-preview`/`-stub-run` or a real `-resume` test run with
`params.heterogeneity_analysis = [MNC: true]` (or similar) against the real
CRUK_DDX41_DHX34 dataset to confirm: (1) `HETEROGENEITY_ANALYSIS` shows `1 of 4` (only MNC),
(2) disabled-group reports (NB4/MSCline/MSC) are byte-identical to a pre-feature run, (3)
the MNC report's new "Heterogeneity Analysis" section renders without error and the
"altered but not DE-significant" table is non-trivial for at least one contrast.

---

# Internal Audit Log
## Status: IN PROGRESS (core reporting bug root-caused and fixed; broader audit incomplete — see Handoff)

Audit of the lncRNA-nextflow downstream pipeline (runs after nf-core/rnaseq).
Scope: DISCOVERY phase, QUANTIFICATION phase, per-contrast reporting bug, resource sizing.
Test data reference (read-only, no new runs): /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/lncRNA

Hard constraints respected: no `nextflow run`, no git commits, all comments/logs in English.

## Findings

### F1 — CRITICAL: only 1 of N norm_groups ever reaches Phase B (QUANTIFICATION)
- File: [main.nf](main.nf), [workflows/quantification.nf](workflows/quantification.nf)
- What is wrong: `tx2gene_detailed` and `gene_catalog` are DISCOVERY-phase process outputs
  (single execution, plain queue channel with exactly one emission). They are passed
  directly — without `.first()` / `Channel.value()` — into `QUANTIFICATION()`, which is
  invoked with `ch_per_group`, a queue channel carrying one tuple per norm_group (N>1).
  Per Nextflow semantics ("Multiple inputs" docs), when a process receives two or more
  **queue channels**, execution is limited to the channel with the *fewest* items — extra
  items in the larger channel are silently dropped (equivalent to `merge()`, not a
  cross-join). Concretely:
  - `TXIMPORT_AND_FILTER(ch_txi_in, tx2gene_detailed)` in [modules/local/tximport_filter.nf](modules/local/tximport_filter.nf)
  - `CIS_ASSOCIATIONS(ch_cis, gene_catalog)` in [modules/local/cis_associations.nf](modules/local/cis_associations.nf)
  - `FUNCTIONAL_ENRICHMENT(ch_enrich, gene_catalog)` in [modules/local/functional_enrichment.nf](modules/local/functional_enrichment.nf)
  - `REPORT_CONTRAST(ch_rpt, ..., gene_catalog)` in [modules/local/report_contrast.nf](modules/local/report_contrast.nf)
  all zip against a 1-item channel and therefore run exactly once, no matter how many
  norm_groups exist.
  - Note the pipeline author was aware of this exact hazard elsewhere: `ch_salmon_index =
    BUILD_DECOY_INDEX.out.index_dir.first()` in [main.nf](main.nf) explicitly guards
    against it. The same guard is missing for `tx2gene_detailed` / `gene_catalog`.
- Evidence: production run log
  `/data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/lncRNA/logs/nextflow_24146169.out`
  shows `NORM_GROUPS_FOR_PHASE_B: [MNC, MSC, MSCline, NB4]` (4 groups reach the split
  step correctly) but the final, completed run shows
  `TXIMPORT_AND_FILTER (1) | 1 of 1`, `DESEQ2_DE (1) | 1 of 1`, `REPORT_CONTRAST (1) | 1 of 1`.
  Only `results/reports/MNC_DDX41_vs_Healthy_report.html` exists on disk — the 8 other
  declared contrasts belonging to NB4/MSCline/MSC (see `comparisons.csv`, 10 contrasts
  total) never ran, and the pipeline still exited 0 ("lncRNA-discovery pipeline
  finished"). This is a bigger and more silent problem than the reported "no per-contrast
  report" symptom: entire norm_groups (3 of 4) are dropped, not just extra contrasts
  within one group.
- Status: fix applied — see Actions Taken.

### F2 — Per-contrast report generation: current code no longer has the literal
  `quarto --output <path>` bug, but risk remains
- File: [modules/local/report_contrast.nf](modules/local/report_contrast.nf)
- The `--output` value (`fname`) is already a bare basename (`gsub(...) + "_report.html"`,
  no directory component), so the specific historical failure mode described ("--output
  only accepts a bare filename") is not reproducible against the current script as written.
- However, `path template` is staged into the task work dir as a **symlink** by default
  (no `stageInMode` override found in any `conf/*.config`). Quarto is known to resolve
  the real path of its input document for some project/output-directory bookkeeping,
  which is a latent risk: if it ever resolves the symlink back to `assets/`, the render
  would write next to the source template instead of into the task directory, and the
  declared `path "*.html"` output would then find nothing.
- Correction: [assets/ct1_report.html](assets/ct1_report.html) / [assets/ct2_report.html](assets/ct2_report.html)
  are intentionally committed example reports (git-tracked, added in the same commit that
  overhauled this module: `20309ca "fix: contrasts bug"`), not accidental leftovers from a
  symlink-escape failure — this is not evidence of the risk actually firing. No cleanup
  needed here; treat the symlink risk as theoretical/defensive only.
- Status: hardened defensively anyway (copy the template into the task dir before
  rendering, so quarto's working directory and the physical file location always
  coincide) — cheap, no behavior change for the already-working case.

### F3 — WARNING: DESeq2 `params.lfc_threshold` / `params.fdr_threshold` silently ignored
  in the final (shrunk) results
- File: [modules/local/deseq2_de.nf](modules/local/deseq2_de.nf)
- What is wrong: `res <- results(dds, contrast=..., lfcThreshold=params.lfc_threshold,
  alpha=params.fdr_threshold)` performs the hypothesis test and independent filtering
  with the user-configured thresholds — but this object was then discarded.
  `res_shrunk <- lfcShrink(dds, coef=..., type="apeglm")` was called **without** `res =
  res`, which means `lfcShrink()` silently recomputed its own internal `results()` call
  using DESeq2 defaults (`alpha = 0.1`, `lfcThreshold = 0`) to get the test statistics
  before shrinking log2FoldChange. `res_df` (the file that actually gets written and
  consumed downstream) came entirely from `res_shrunk`, i.e. from the default-parameter
  statistics, not from the ones actually requested via `params.lfc_threshold` /
  `params.fdr_threshold`. Net effect: those two pipeline parameters exist in
  `nextflow_schema.json` / `params.yaml` and look like they control significance calling,
  but had zero effect on the exported `pvalue`/`padj`/independent-filtering behaviour.
- Evidence: read the `apeglm`/`lfcShrink` documentation contract — `res` is the
  documented way to keep `lfcShrink()`'s output consistent with a prior `results()` call
  that used non-default `alpha`/`lfcThreshold`. Without it, `lfcShrink()` re-derives its
  own results object.
- Status: fix applied — added `res = res` to both the apeglm and ashr `lfcShrink()`
  calls in `deseq2_de.nf`.

### F4 — WARNING: cis-association DE annotation is contrast-blind in multi-contrast groups
- File: [modules/local/cis_associations.nf](modules/local/cis_associations.nf) (~line 178)
- What is wrong: `CIS_ASSOCIATIONS` runs once per norm_group (not per contrast) and
  receives the group's full, multi-contrast `deseq2_summary.tsv` (one stacked table with
  a `contrast` column, one block of rows per contrast — see `deseq2_de.nf`). To annotate
  each significant cis pair with a DE log2FC/padj, the code does:
  `de_idx <- match(cis_sig$lncrna_id, de$gene_id)` — `match()` only returns the *first*
  matching row. The code comment even says so: "DE is per-contrast; annotate the lncRNA
  side on the first match". For any norm_group with more than one contrast (e.g. NB4 has
  4: DDX41sh1/DDX41sh2/DHX34sh1/DHX34sh2 vs Scramble), `lncrna_DE_log2FC`/`lncrna_DE_padj`
  in `${group}_cis_pairs_significant.tsv` come from an arbitrary single contrast (whichever
  appears first in the stacked table) — yet this same file is reused, unchanged, by
  `REPORT_CONTRAST` for **every** contrast's report of that group. So 3 of NB4's 4
  per-contrast reports would show cis-association DE numbers that belong to a different
  contrast than the one the report is about, with no indication that this happened.
- Evidence: read `cis_associations.nf` lines ~176-182 and the calling convention in
  `workflows/quantification.nf` (`CIS_ASSOCIATIONS` is invoked once per group, joined
  only on `group`, never on `contrast`).
- Status: found, NOT fixed. This needs a design decision rather than a mechanical patch:
  either (a) make `CIS_ASSOCIATIONS`'s DE annotation contrast-aware (would require
  restructuring the module and its channel wiring to run per-contrast, a real scope
  change), or (b) drop/rename the ambiguous `lncrna_DE_*` columns and let
  `report_template.qmd` re-join against the contrast-specific `de_summary` it already
  receives (cheaper, recommended — the report already loads `de_ann` filtered to
  `contrast == params$contrast`, so it could look up the lncRNA's own DE stats from
  there instead of trusting the pre-computed, ambiguous columns). Flagging for another
  pass rather than guessing; left untouched pending a decision.

## Actions Taken

- [internal_audition.md](internal_audition.md) created and kept up to date incrementally.
- [workflows/quantification.nf](workflows/quantification.nf): added `.first()` on
  `tx2gene_detailed` and `gene_catalog` at the top of the `main:` block (F1 fix), and
  updated `TXIMPORT_AND_FILTER`, `CIS_ASSOCIATIONS`, `FUNCTIONAL_ENRICHMENT`, and
  `REPORT_CONTRAST` call sites to use the guarded channels. Not yet re-verified against
  a real run (hard constraint: no `nextflow run`); reasoning is based on documented
  Nextflow "Multiple inputs" semantics + the observed 1-of-4-groups log evidence.
- [modules/local/report_contrast.nf](modules/local/report_contrast.nf): copy the staged
  `template` to a local `report_template.qmd` in the task dir before invoking `quarto
  render`, so quarto's output-directory resolution can never depend on symlink behaviour
  (F2 hardening).
- [modules/local/deseq2_de.nf](modules/local/deseq2_de.nf): pass `res = res` to both
  `lfcShrink()` calls (apeglm + ashr fallback) so shrinkage reuses the `alpha`/
  `lfcThreshold`-aware test statistics from `results()` instead of silently recomputing
  with DESeq2 defaults (F3 fix).
- [assets/report_template.qmd](assets/report_template.qmd): fixed the `session` chunk —
  it computed `sessioninfo::session_info()`/`sessionInfo()` into an unused variable
  (`devtools_or`) and never printed it; now it's actually printed at the end of the
  report.
- F4 (cis-association DE annotation is contrast-blind) documented but **not** fixed —
  needs a design decision, see Handoff.

## Context: this is a same-day, in-progress fix
`git log -1` shows the working tree's parent commit is `20309ca "fix: contrasts bug"`
(today), which already reworked `workflows/quantification.nf` from a `multiMap`-based
per-group split (which risked scrambling groups via non-keyed pairing) to the current
`.join()`-on-group-key design, and introduced `gene_catalog` (new module
`build_gene_catalog.nf`) to replace passing `frozen_gtf` + `reference_gtf` separately
into `CIS_ASSOCIATIONS`. That rework correctly fixed the intra-group key-matching
problem, but it (re)introduced F1: `tx2gene_detailed` and `gene_catalog` are
single-execution Discovery outputs consumed directly by per-group processes without
`.first()`. So the "only one contrast/report gets produced" symptom the user is chasing
is best explained by F1, not by a quarto `--output` path issue (see F2).

## Handoff — Next Actions For Another Agent

Read so far (this session): main.nf, workflows/validate_input.nf, workflows/discovery.nf,
workflows/quantification.nf, modules/local/split_normgroup.nf, modules/local/report_contrast.nf,
modules/local/tximport_filter.nf, modules/local/cis_associations.nf (full),
modules/local/deseq2_de.nf (full), assets/report_template.qmd (full), conf/modules.config
(full), conf/apocrita.config (partial). Cross-referenced against the real run at
CRUK_DDX41_DHX34/.../lncRNA (logs, results, comparisons.csv, samplesheet.csv, git history).

Not yet reviewed: modules/local/{functional_enrichment,gffcompare,cpat,cpat_build_model,
cpc2,coding_consensus,catalog_overlap,expression_recurrence,annotation_freeze,
build_gene_catalog,stringtie_per_sample,stringtie_merge,featurecounts_expression,
salmon_quantify,build_decoy_index,infer_strandedness,bam_quickcheck,length_filter}.nf,
conf/base.config, conf/slurm.config, conf/test.config, nextflow.config, nextflow_schema.json.

Fixes applied this session (F1, F2, F3) are believed correct based on Nextflow/DESeq2/
quarto documented semantics, but could NOT be verified with an actual `nextflow run`
(hard constraint). Recommend the next reviewer/human re-run once with `-resume` against
the same real dataset and confirm:
  - `TXIMPORT_AND_FILTER`/`DESEQ2_DE`/`CIS_ASSOCIATIONS`/`FUNCTIONAL_ENRICHMENT`/
    `REPORT_CONTRAST` now show `4 of 4` (one per norm_group: MNC, MSC, MSCline, NB4),
    and `results/reports/` ends up with 10 files (one per row of `comparisons.csv`
    whose numerator/denominator both exist in that group's coldata).
  - DESeq2 output columns (`pvalue`/`padj`) shift slightly vs the old cached results now
    that `lfcShrink(res = res)` is wired in (expected, this is the fix).

Next steps in priority order:
1. Decide and implement a fix for F4 (cis-association DE annotation is contrast-blind
   in multi-contrast norm_groups) — recommended approach: drop `lncrna_DE_log2FC`/
   `lncrna_DE_padj` from `cis_associations.nf`'s output entirely (it can't be correct at
   the group level) and instead have `report_template.qmd` look up
   `de_ann$log2FoldChange[match(cis_sig$lncrna_id, de_ann$gene_id)]` itself, since
   `de_ann` there is already filtered to the current `contrast`.
2. `functional_enrichment.nf` was not reviewed — check for the same "is this per-contrast
   or per-group, and does everything downstream know which one it is" issue that F4
   surfaced in `cis_associations.nf` (functional enrichment is also fed the full,
   multi-contrast `deseq2_summary`, per `quantification.nf`).
3. Audit gffcompare.nf / coding_consensus.nf / catalog_overlap.nf / annotation_freeze.nf /
   build_gene_catalog.nf for the "silent correctness bugs" the original brief called out
   specifically (class-code handling, ID consistency between tx2gene and gene_catalog —
   `report_template.qmd`'s `catalog-summary` chunk already contains a defensive check for
   >5% of tested genes missing from the catalog, worth confirming it never actually fires
   on real data, i.e. that gene IDs are consistent end-to-end).
4. Resource sizing (conf/modules.config): spot-checked only; SALMON_QUANTIFY/BUILD_DECOY_INDEX/
   STRINGTIE2_PER_SAMPLE look generously sized for 50 real samples + human genome, no
   obvious undersizing found, but conf/base.config and conf/slurm.config (queue/time
   limits) were not opened this session — worth a pass in case wall-time limits are too
   tight for the real per-sample BAM sizes in the CRUK_DDX41_DHX34 dataset.
5. conf/test.config and nextflow_schema.json not reviewed — worth confirming the
   documented `params.lfc_threshold`/`params.fdr_threshold` defaults are sane now that F3
   makes them actually take effect.
