# lncRNA Downstream Analysis Pipeline — Design Document 

**Scope:** disease-agnostic Nextflow pipeline for downstream lncRNA discovery, curation, decoy-aware quantification, and differential expression, run after `nf-core/rnaseq`. The current application is a total-RNA (ribodepleted, Novogene) AML cohort comparison for publication, but nothing below is hard-coded to AML — see the callouts marking where disease-specific inputs plug in.

**How to read this document:** This document outlines an evolutionary implementation plan. The immediate objective is to fully build and validate the **CORE** path, which constitutes a complete, self-contained, and publishable pipeline on its own. **OPTIONAL** modules are explicitly deferred to future iterations; they represent potential extensions for when additional multi-omics data or external cohorts become available, and they are intentionally excluded from the current implementation scope.

## Contents

0. Scope & Design Principles
1. What Changed vs. the v1 Draft
2. Repository Layout
3. Inputs & Validation
4. Containerization & CI/CD
5. Reporting Engine
6. Phase A — Discovery, Curation & Freezing (A1–A7)
7. Phase B — Quantification & Statistics (B1–B8)
8. Phase C — Orthogonal & External Validation (C1–C2)
9. QC / Sign-off Checklist
10. Explicitly Out of Scope
11. Design Rationale & Sources Consulted

---

## 0. Scope and Design Principles

- **Core principle (unchanged from v1):** strict separation of **Phase A — Discovery, Curation, Freezing** from **Phase B — Quantification & Statistics**. This **prevents condition-driven selection and substantially reduces selection bias; discovery still utilizes the same biological data unless an independent cohort is provided.** This is the single most important methodological guard in the design; nothing below compromises it.
- **Generalizability:** every disease-specific input (clinical covariates, a normal-tissue reference panel, an external validation cohort) is a *configurable slot*, not pipeline code. The AML details throughout are worked examples to make the design concrete, not requirements.
- **Baseline:** conceptually informed by Carrasco-Leon et al. 2021 (*Leukemia* 35:1438–1450), a ssRNA-seq lncRNA landscape study in multiple myeloma, as a proof-of-concept that this class of analysis is publishable and reviewer-tested. Several of that paper's 2019–2020-era tool choices are superseded below — see §11.

## 1. What Changed vs. the v1 Draft

| # | Gap identified | Fix in this version (v6.1) |
|---|---|---|
| 1 | No control for intronic/gDNA contamination in monoexonic transcripts | Intron-coverage ratio converted to a configurable risk flag rather than hard exclusion (§A3.4) |
| 2 | "x" (antisense) class accepted with no artifact control | Highly-expressed-neighbor check converted to an artifact risk flag + empirical strandedness via RSeQC (`infer_experiment.py`) mapping to Salmon explicit parameters (§A2, §A3.3) |
| 3 | Coding-potential filter had no conservation or known-catalog layer | PhyloCSF (species-permitting) + catalog overlap, reported not filtered (§A4–A5) |
| 4 | Disease covariates unaddressed | Config-driven covariate columns, AML given as the worked example (§3, §B4) |
| 5 | "LFC shrinkage" left unspecified | Explicit apeglm/ashr rule by contrast type, clearly separating statistical significance (FDR) from magnitude of effect (§B4) |
| 6 | No multiple-testing correction on cis-association p-values | BH-FDR mandated on the full family of tests per contrast (§B7) |
| 7 | Silently assumed nf-core/rnaseq produced genome BAMs | Explicit early validation gate via content-checking (`samtools quickcheck`) (§A1) |
| 8 *(new)* | — | Optional heterogeneity/recurrence diagnostic, explicitly secondary to DESeq2 (§B5) |
| 9 *(new)* | — | Optional disease-vs-normal-lineage-continuum module, generalized from the MM paper (§B6) |
| 10 *(new)* | — | Optional epigenomic corroboration module (§C1) |
| 11 *(new)* | — | Optional external-cohort validation with explicit multiplicity control (§C2) |

## 2. Repository Layout
.
├── main.nf
├── nextflow.config
├── conf/
│   ├── modules.config
│   └── base.config
├── modules/local/        # one .nf file per analytical step
├── workflows/            # phase-level sub-workflows (phase_a.nf, phase_b.nf)
├── bin/                  # R / Python / Bash scripts, chmod +x
├── assets/               # Quarto/RMarkdown templates, CSS, reference manifests
├── tests/                # nf-test on a minimal/subsampled dataset
└── .github/workflows/    # lint, test, container build → GHCR

## 3. Inputs & Validation (CORE)

- **Samplesheet:** `sample, condition, batch, subject, fastq_1, fastq_2`, plus any number of user-defined covariate columns. The pipeline treats these as opaque columns usable in the design formula — it never hardcodes a disease-specific field name.
  *AML worked example:* `cytogenetic_risk, blast_pct, fab_subtype`. Substitute freely for other diseases.
- Validate at startup: required columns present, no duplicated sample IDs, absolute paths resolve, FASTQ pairing correct.
- **Confounding check:** compute `rank(model.matrix(<user formula>))` vs. the number of coefficients; hard-fail on rank deficiency (complete confounding), warn on near-singularity (e.g. a batch level with <2 samples in some condition).
- **Contrast matrix:** CSV of named contrasts (e.g. `KO_vs_WT`); every referenced level is validated against the samplesheet before the run starts, not at the DESeq2 step.

## 4. Containerization & CI/CD (CORE)

- Bioconductor/Python base images; every tool version-pinned, no `latest` tags.
- GitHub Actions: `nf-core lint`, unit tests on a minimal/subsampled dataset, container build + push to GHCR on tag.

## 5. Reporting Engine (CORE)

- Parameterized Quarto/RMarkdown → self-contained HTML (base64 assets, `plotly`/`DT`).
- Per-contrast fault isolation: one failed contrast yields a flagged, empty report section, not a broken build.
- Run manifest: pipeline version, git commit, container digests, checksums of the Phase A frozen files.

---

## 6. PHASE A — Discovery, Curation & Freezing (CORE)

### A1. Upstream Compatibility Gate (CORE)

- Confirm the upstream `nf-core/rnaseq` run used `--aligner star_salmon` (pipeline default) or `hisat2` — both produce genome-coordinate BAMs suitable for StringTie. Guard explicitly against `--aligner bowtie2_salmon` (transcriptome-only pseudo-alignment), which does not produce genome BAMs.
- **Content-based BAM validation:** Execute `samtools quickcheck` on every input BAM to verify structural integrity, valid headers, and EOF blocks, avoiding silent failures from corrupted files. **Hard-fail with a clear message** if validation fails.
- Verify `genome.fa`, the reference GTF, and the upstream BAM headers share identical chromosome naming and coordinate build.

### A2. Transcriptome Assembly & Strandedness Detection (CORE)

- **Empirical strandedness check:** Run RSeQC's `infer_experiment.py` on a subsample of each BAM to determine library strand specificity. Feed this empirical strand orientation directly into downstream tools rather than relying on guesswork or default `-l A` parameters in Salmon.
- Per-sample **StringTie2** → `stringtie --merge` against the reference GTF → single discovery GTF.

### A3. Structural Filtering & `gffcompare` Classification (CORE)

1. Length filter: summed exon length > 200 nt.
2. `gffcompare` vs. the reference GTF. High-confidence novel-locus classes: **u** (intergenic), **i** (intronic), **x** (antisense). Novel-isoform-of-known-gene classes (**j, o, c, k, m, n**) tracked separately. Artifact-prone classes (**e, p, s, r**) excluded from the high-confidence set, kept only for sensitivity analysis. Class **y** (rare) reviewed manually rather than auto-included or auto-excluded.
3. **Antisense ("x") artifact guard:** Flag any "x" transcript whose opposite-strand overlapping gene falls in the top expression percentile of the cohort (default top 5%, configurable) as an *antisense artifact risk* rather than performing a blind automatic exclusion.
4. **Monoexonic handling:** Compute an intron-coverage ratio per candidate (`bedtools coverage` or `featureCounts` over the transcript span vs. any overlapping intron of a neighboring gene). Candidates above a configurable ratio (default 0.3) are flagged with high gDNA/pre-mRNA contamination risk for sensitivity evaluation rather than silent deletion.
5. Exclude pseudogenes, rRNA, mitochondrial transcripts, and high-repeat-content loci via strict overlap checks.

### A4. Coding-Potential Consensus (CORE)

- **CPAT** (species-specific pre-built model with validated cutoff) **AND CPC2** (< 0.5). Require agreement between both for the "high-confidence non-coding" tier; report UpSet/Venn concordance in the QC report.
- *(Note: PhyloCSF is deferred to future optional extensions).*

### A5. Orthogonal Catalog Evidence (CORE, reporting-only — never a filter)

- Overlap candidates against NONCODE, LNCipedia, RNAcentral, and FANTOM CAT. Tags `known_catalog_hit` vs. `putatively_novel` in the metadata table.

### A6. Reproducibility / Recurrence Filter (CORE)

- Require expression support in a configurable minimum number of biological replicates **across the whole experiment**, condition-blind.

### A7. Annotation Freezing (CORE)

- Compile the final custom GTF and `tx2gene_detailed.tsv`.
- Checksum and tag `FROZEN`. **Strict technical barrier:** Phase B receives *only* the frozen artifacts and their checksums from Phase A. Any change to the GTF automatically increments the catalog version hash.

---

## 7. PHASE B — Quantification & Statistics (CORE)

### B1–B2. Unified Reference, Whole-Genome Decoys, Decoy-Aware Salmon Index & Requantification (CORE)

- Extract frozen novel lncRNA sequences (`gffread`) → concatenate with full GENCODE transcript set → append whole-genome decoys + `decoys.txt` → build decoy-aware Salmon index → requantify samples using the empirically determined strand parameter (e.g., `-l SF` or `-l SR` derived from A2, replacing loose `-l A`).

```bash
salmon quant -i decoy_aware_index -l <empirical_strand> \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --gcBias --seqBias --validateMappings \
  -p 8 -o sample_quant
```

### B3. Import & Independent Filtering (CORE)

* **Quantification Import:** `tximport` (`type="salmon"`) utilizing the frozen `tx2gene_detailed.tsv` mapping file.
* **Independent Filtering:** Condition-blind low-count filtering executed via `DESeqDataSetFromTximport` or `edgeR::filterByExpr()`.
* **Sample-Level QC:** Execution of PCA and sample-sample correlation checks prior to differential expression testing to identify and isolate potential outliers or technical artifacts.

### B4. Differential Expression — DESeq2 (CORE)

* **Flexible Design Formulas:** Supports user-defined formulas (e.g., `~ batch + condition`, paired designs like `~ subject + condition`, or custom clinical covariates). Batch effects are modeled directly as GLM covariates.
* **Shrinkage Decision Rule & Effect Size Separation:**
  * Strict separation between statistical significance (FDR < 0.05) and the magnitude of biological effect sizes.
  * Utilization of `lfcThreshold` inside DESeq2 when testing against a minimum biological log2 fold-change requirement.
  * Application of `apeglm` for direct model coefficients or `ashr` for complex contrasts/interactions during LFC shrinkage.

### B5. Cis-Regulatory Candidate Associations (CORE)

* **Spatial Search:** `GenomicRanges`-based window search (\pm 100 kb default) to identify protein-coding genes located near lncRNA loci.
* **Correlation Metrics:** Computation of global, within-condition, and partial correlation coefficients.
* **Rigorous Multiplicity Correction:** Mandatory Benjamini-Hochberg (BH) FDR correction applied across the **entire family of tested lncRNA–gene pairs per contrast**. Outputs are labeled strictly as **candidate cis-regulatory associations**.

### B6. Functional Enrichment (CORE)

* **Target Enrichment:** Execution of Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) via `clusterProfiler` on associated protein-coding targets.
* **Background Universe:** Utilizes the experiment's own expressed-gene universe as the proper statistical background rather than a generic genome-wide default.

---

## 8. PHASE C — Orthogonal & External Validation (DEFERRED / OUT OF SCOPE FOR CURRENT PLAN)

*All modules in this phase (Epigenomic Corroboration and External Cohort Validation) are deferred to future roadmap phases and are not part of the current implementation scope.*

---

## 9. QC / Sign-off Checklist

- [ ] Samplesheet valid; confounding/design-singularity check passed or explicitly justified
- [ ] Contrast matrix conditions verified against the samplesheet
- [ ] Upstream BAMs validated for content integrity via `samtools quickcheck`
- [ ] `genome.fa` / GTF / upstream BAMs coordinate-system-matched
- [ ] Empirical strandedness determined via RSeQC and passed explicitly to Salmon (no unguided `-l A`)
- [ ] Monoexonic and antisense ("x") features evaluated as risk-flagged metadata rather than rigid deletions
- [ ] Coding-potential consensus (CPAT+CPC2) computed and summarized via UpSet/Venn visualization
- [ ] Known-catalog overlap tagged (not filtered)
- [ ] Novel-lncRNA catalog frozen, checksum-verified, and version-locked before Phase B execution
- [ ] Independent filtering (§B3) confirmed condition-blind; sample-level QC reviewed
- [ ] Statistical significance (FDR) decoupled from effect size criteria; shrinkage method matches contrast type
- [ ] Cis-association p-values BH-corrected across the full test family before designation as candidate associations
- [ ] Enrichment background = experiment's expressed-gene universe
- [ ] CI/CD: `nf-core lint` and container build pass on push

## 10. Explicitly Out of Scope

* Wet-lab functional validation (shRNA knockdown, qPCR, phenotypic assays).
* Single-cell RNA-seq (bulk analysis only).
* Optional Phase C extensions (Epigenomic datasets integration, independent clinical cohort prognostic modeling).

## 11. Design Rationale & Sources Consulted

* **Carrasco-Leon et al. 2021 (*Leukemia*):** Serves as the conceptual workflow framework (discovery -> downstream analysis), upgraded with modern tools (StringTie2, CPAT+CPC2 consensus, decoy-aware Salmon quantification, strict multiple-testing controls).
* **github.com/GudaLab/Lncrna-framework:** Reviewed and rejected due to outdated single-sample scripting architecture, flawed reference indexing methods, and lack of comprehensive decoy protection.