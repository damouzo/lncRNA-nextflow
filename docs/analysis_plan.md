# Pre-specified analysis plan — lncRNA catalog tiers and flags

**Status:** F1.1a frozen 2026-10-07. F1.1b (numeric thresholds) is **pending** the
final control run and must be signed before any Phase B DE result is opened.
Post-pilot revision 2026-10-14: F1.3 evaluated and **not adopted** (§5); mono-exonic
demotion is now the categorical rule (§4); F1.4 host resolution and strand convention
fixed (§6).
**Scope:** documents the rules that define which transcripts/genes are reported
as lncRNA candidates, and the procedures that fix the flag thresholds. It is a
contract, not a results document.

## 0. Honesty statement

DE results from the legacy Salmon run (library type `IU`, i.e. unstranded because
of the `-l SR` single-end bug) have already been inspected. **No threshold in this
plan is derived from those results.** Every threshold is fixed from the controls
described in §5 and §6 (annotated mono-exonic lncRNAs, introns of highly expressed
genes, annotated host–antisense pairs) and from the strandedness fidelity measured
in `run_info/strandedness_summary.tsv`. The legacy DE is only used later, in F6.4,
as a pre/post comparison.

The rules below are categorical wherever the property is not library-specific
(transcript class, exon count, parent biotype). The only data-derived number is the
F1.4 antisense threshold, which comes from the empirical leakage null of each
library and carries a fallback and a recorded origin. Categorical rules were fixed
without reference to any DE result.

## 1. Universe

Rules apply to the frozen discovery catalog (`frozen_lncrna.gtf`, `gene_catalog.tsv`)
produced by Phase A. Class codes come from gffcompare (`class_code`). The reference
parent gene for tier 1b (`j`/`o`) is the `ref_gene_id` attribute; the F1 antisense/
host resolution uses the reference transcript in `cmp_ref` mapped through the
reference tx2gene, because `cmp_ref_gene` holds a gene *name* (or an Ensembl id only
when the reference gene has no name). Nothing here changes the frozen GTF or the
assigned IDs.

## 2. Transcript-level tier (`novel_type`)

Pure function of two catalog columns: `class_code` and the biotype of the
reference parent gene (`parent_biotype`, resolved through `ref_gene_id`; present
for the classes that need it, `j`/`o`).

| `novel_type` | Rule |
|---|---|
| `1a` | `class_code ∈ {u, x, i}` — new loci |
| `1b` | `class_code ∈ {j, o}` **and** the reference parent gene biotype is in the whitelist `tier1_parent_biotypes` (default `[lncRNA]`) |
| Tier 2 | everything else: `j`/`o` whose parent is absent or not in the whitelist (pseudogenes, TEC, small-RNA hosts, …), and `c, k, m, n` |

`j`/`o` are isoforms of known genes and are mostly coding; keeping them in Tier 1
without the parent-biotype whitelist would let "novel lncRNA" mean "isoform of a
coding gene". `i` (intronic) is kept in `1a` per the researcher decision.

Parameters: `tier1_classes` (default `[u,x,i,j,o]`), `tier1a_classes` (default `[u,x,i]`),
`tier1_parent_biotypes` (default `[lncRNA]`). The gene catalog carries `parent_biotype`.

## 3. Gene-level aggregation

One row per `gene_id`:

- `n_transcripts` — novel transcripts of the gene in the frozen GTF.
- `n_tier1_transcripts` — those that are Tier 1 by class **and** pass the
  mono-exonic rule (§4).
- `frac_tier1_transcripts` — `n_tier1_transcripts / n_transcripts`.
- `tier_static` — `n_tier1_transcripts ≥ 1` (best-evidence rule: one clear
  transcript is enough; "majority" would hide genes whose best isoform is clear).
- `single_exon_low_confidence` — any transcript demoted by the §4 mono-exonic rule.
- `novel_type` (gene) — `1a` if any transcript is `1a`, else `1b` if any is `1b`,
  else `NA`. Class-only label; it does not encode the mono-exonic rule, so a `1a`
  gene can be `tier_static = FALSE` when all its `1a` transcripts are mono-exonic.
  Reference-origin rows are `NA` here.

Reproducible from the columns with a one-line rule; no filtering (annotate only).

## 4. Mono-exonic rule

**Categorical rule (no derived threshold).** A mono-exonic transcript of any class
is demoted to Tier 2 and marked `single_exon_low_confidence`, unless
`tier1_allow_single_exon` is set (default `false`). The rule is applied at the
transcript level before gene aggregation: `tier_static(tx) = class tier AND
(multi-exonic OR tier1_allow_single_exon)`. A gene survives only if at least one
of its Tier 1 transcripts is multi-exonic (or the parameter is set).

Rationale: mono-exonic transcripts reconstructed in ribodepleted total RNA are the
most pre-mRNA/gDNA-prone calls, and a categorical structural rule is dataset-
independent and needs no tuning. Nothing is discarded — every gene stays in the
catalog with the `single_exon_low_confidence` mark; a sensitivity run
(`tier1_allow_single_exon=true`) reports how many genes would return to Tier 1.

`u_single_exon_low_confidence` (class `u`, gene level; pre-existing) stays
annotation-only and is **not** a filter.

## 5. F1.3 — flank-coverage flag (gDNA / pre-mRNA)

**Outcome: evaluated, not adopted.** On the 2026-10 pilot controls the derived
threshold failed its pre-specified validation: the `neg_flank` rejection rate at
`t_flank` was **0.135 < 0.90** (`run_info/thresholds_used.tsv`), i.e. the metric does
not separate positives from negatives. The controls were **not** redesigned to force
a pass. The check is now **exploratory and off by default** (`run_flank_check=false`)
and does **not** gate `tier_final`; the categorical mono-exonic rule (§4) replaces
it. The method below is documented so the code can be re-evaluated on other data;
anything in `results/f1/` from F1.3 is exploratory.

The rest of this section describes the (retained but disabled) F1.3 method.
Module `FLANK_COVERAGE_RATIO` (Phase A, post-freeze). For each mono-exonic frozen
Tier 1 transcript, strand-specific paired counts (featureCounts `-p
--countReadPairs`) over SAF windows: the candidate, its left and right flanks
(±`f1_flank_bp`), on the same strand (`-s 2`) and on the opposite strand (`-s 1`).
Flanks are **masked** against reference-gene exons (a flank that overlaps another
gene's exon only contributes its unmasked pieces); a flank with less than
`f1_min_flank_len` unmasked bases is `NA`. For class `i` the flank is the **host
intron** (coding host; a non-coding host falls back to genomic flanks).
Boundary reads are counted per feature by featureCounts, so a read crossing the
candidate/flank boundary can contribute to both; this is documented, not corrected.

Output per element and norm_group: `flank_ratio_median`, `opp_strand_ratio_median`,
`n_samples_used`, `flank_len_used`. A sample contributes only with ≥ `f1_min_reads`
same-strand candidate reads.

- `flank_ratio = density(candidate) / mean(density(flanks), same strand)`
- `opp_strand_ratio = opposite-strand candidate reads / same-strand candidate reads`

gDNA gives `flank_ratio ≈ 1`; a real mono-exonic transcript gives `flank_ratio ≫ 1`.

**Threshold rule (fixed, derived in-run by `DERIVE_THRESHOLDS`):**
`t_flank = percentile(f1_flank_percentile=10) of pos_flank.flank_ratio_median`,
i.e. the value that keeps ≥ 90% of the positive controls. The neg_flank rejection
rate at `t_flank` is reported. `flank_ratio_threshold` (param) overrides. If a
control set has fewer than `f1_min_controls` elements, or neg_flank rejection
< `f1_neg_reject_min` (0.90), the run falls back to `flank_ratio_default` and marks
`thresholds_unreliable` (propagated to the tier tables).

Controls (deterministic: seed `f1_seed`, cap `f1_control_cap`), processed by the
same module with a `set` column:
- `pos_flank` — annotated mono-exonic expressed lncRNAs (expression proxy from the
  Phase A recurrence table).
- `neg_flank` — intronic windows of highly expressed coding genes, length-matched
  to the candidates, avoiding reference-gene exons.

**"Highly expressed" definition and bias.** The only condition-blind expression
signal available in Phase A is the recurrence table, which is built from the
candidate set (gffcompare `=` / novel transcripts). "Highly expressed coding gene"
therefore means a protein-coding gene with at least one transcript of the
candidate universe passing recurrence (`n_reps_expressed ≥
expression_recurrence_min_reps`), ranked by the best transcript's
`n_reps_expressed`. Consequences:
- Coding genes absent from the candidate universe cannot be selected; the control
  set is biased towards genes with a recovered reference transcript.
- Expression is a replicate-count proxy, not a magnitude, so ranking is coarse.
The negative (`neg_flank`) is only used to measure the rejection rate at the
derived threshold and, if unreliable, triggers the default fallback; it never sets
the threshold itself (positives do). The bias is therefore a sensitivity
limitation, not a bias in the threshold value.

## 6. F1.4 — antisense host-leakage flag

Module `FLAG_ANTISENSE_LEAKAGE` (Phase A, post-freeze; redesigned to avoid Salmon
and batch/VST dependencies). For each class-`x` candidate, the host is the
overlapping reference gene (`cmp_ref_gene`; `ref_gene_id` is absent for class `x`)
and the featureCounts region is the **exonic overlap** of the candidate with the
host's coding exons. The SAF feature is placed on the **host strand**, so `-s 2`
counts host (sense) reads and `-s 1` counts candidate (antisense) reads — the same
convention as the `null_antisense` control.

Main metric: `ratio = antisense_reads / sense_reads = candidate / host` in that
region. Secondary: Spearman `rho` between sense and antisense counts within the
norm_group (reported with `n`; not used for the flag).

Output per element and norm_group: `region_len`, `sense_reads`, `antisense_reads`,
`ratio_median`, `rho`, `n`.

**Threshold rule:** `r_thr = percentile(f1_antisense_percentile=95) of
null_antisense.ratio_median` (the empirical leakage floor). A candidate is
`x_host_leak_suspect` when `ratio_median ≤ r_thr` (its antisense signal does not
clear the floor). `antisense_ratio_threshold` (param) overrides; unreliable controls
fall back to `antisense_ratio_default`.

**Specificity diagnostic (does not change the value).** The fraction of
`pos_antisense` flagged at `r_thr` is reported as the antisense `separation` column.
If the positives are not above the null (`pos_median ≤ null_median`) the status is
`thresholds_inverted`; if more than `f1_antisense_pos_flag_max` of positives are
flagged it is `thresholds_not_specific`. Both are surfaced rather than silently
tolerated.

Controls:
- `pos_antisense` — annotated lncRNA / coding-gene opposite-strand exonic overlaps,
  SAF on the coding (host) strand.
- `null_antisense` — exons of expressed coding genes with no annotated element on
  the opposite strand (the leakage null), SAF on the coding gene's own strand (the
  host). "Expressed" uses the same recurrence-based definition and the same bias as
  §5: only coding genes present in the candidate universe are eligible, and only
  their clean exons (no annotated opposite-strand overlap in the reference) are
  used. The null sets the leakage floor `r_thr`, so a systematically "quiet" null
  would make the floor too permissive; this is checked through the control-set size,
  the `pos_antisense` separation and the recorded distribution in
  `run_info/f1_controls/`.

## 7. `tier_final`

`tier_final` is computed at the **transcript** level and then aggregated to gene
(Phase A, module `COMPUTE_TIER_FINAL`, terminal — nothing downstream depends on it):

```
tier_static(tx) = class tier (1a/1b) AND (multi-exonic OR tier1_allow_single_exon)
tier_final(tx)  = tier_static(tx) AND not antisense-leakage flag
gene tier_final = (n_tier1_transcripts_final >= 1)
```

`tier_static` already contains the mono-exonic rule (§4). `tier_final` adds only the
F1.4 antisense flag; F1.3 (flank) is exploratory and does **not** participate. A
flag marks the element if it triggers in at least one norm_group. The outputs
(`tier_final_transcripts.tsv`, `tier_final_genes.tsv`) reproduce `tier_final` from
the columns with the rule above and are consumed by F4.2/the report. They also carry
`single_exon_low_confidence`, `flank_*` (exploratory) and `antisense_*` columns so
every demotion is traceable. A post-run check
(`bin/tests/check_tier_consistency.py`) asserts `tier_static` and
`n_tier1_transcripts` agree between `gene_catalog.tsv` (R) and `tier_final_genes.tsv`
(Python), so the duplicated rule cannot drift.

## 8. Sensitivity analysis (F1.6)

Report the number of Tier 1 genes under 2–3 alternative antisense thresholds
(`tier_final_sensitivity.tsv`, percentiles p90/p95/p99 of the null) and, separately,
the number recovered with `tier1_allow_single_exon=true` (the mono-exonic rule is
categorical, so this is a switch, not a threshold sweep). Table shipped in
`deliverables/catalog/`.

## 9. Ancillary decisions

- **Coding-potential cutoffs:** `CODING_CONSENSUS` intersects CPAT and CPC2. CPAT
  uses `cpat_coding_threshold = 0.364` (CPAT's own human default) and CPC2 the tool
  default `< 0.5`; both are the tools' defaults, not tuned here. Requiring both to
  agree is the stricter of the two choices and is documented as such.
- **DE FDR:** `padj` comes from the full universe of tested genes (DESeq2
  independent filtering); it is **not** re-adjusted within Tier 1. Tier 1 genes are
  reported as "genes with `padj < fdr_threshold`" from that universe. Selecting a
  subset and recomputing FDR would define a different test universe and would require
  its own pre-specification.
- **FEELnc:** opt-in, annotation-only, run on the same set as CPAT/CPC2, **outside**
  Tier 1 and outside `CODING_CONSENSUS` (frozen GTF and IDs unchanged). Time cap 2 days.
- **External catalogs:** LNCipedia + NONCODE minimum; RNAcentral if an hg38 BED is
  easily available; FANTOM CAT only if hg38. BEDs failing the assembly check are
  dropped or lifted over; version/date/md5 recorded in `run_info/external_catalogs.tsv`.
- **xlsx:** TSV is the source of truth; xlsx is generated alongside for the wet lab,
  with symbol columns forced to text (avoids SEPT2/MARCH1 being read as dates).

## 10. F1.1b — thresholds (TO BE FILLED)

Filled **after** the F1 controls are run (Phase A) and **before** opening any DE.
Values come from `results/run_info/thresholds_used.tsv`; record the origin
(`derived`/`override`/`default`), the control-set sizes and the separation metric.
If `thresholds_unreliable`, record the fallback and the reason. With F1.3 disabled
by default, only the antisense row is present; the flank row is filled only if
`run_flank_check=true`.

| Param | Value | Origin | n controls | separation | notes |
|---|---|---|---|---|---|
| `antisense_ratio_threshold` (derived) | | | null_antisense = | pos flag rate = | r = p95 of null; status from §6 |
| `flank_ratio_threshold` (exploratory) | n/a | disabled | pos_flank = | neg rejection = | filled only if `run_flank_check=true` |
| `tier_final` Tier 1 genes (baseline) | | | | | from `tier_final_genes.tsv` |
| sensitivity p90 / p95 / p99 (antisense) | | | | | `tier_final_sensitivity.tsv` |
| recovery with `tier1_allow_single_exon=true` | | | | | count of genes demoted by §4 |

## 11. Sign-off

F1.1b is considered signed when this table is filled and the document is committed
with a date, before the final run (F7.0).
