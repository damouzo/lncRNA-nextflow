# lncRNA-nextflow

Discovery, curation, and differential expression pipeline for lncRNAs.

Two phases, kept strictly separate so contrast information never influences which transcripts get kept:

- **Phase A — Discovery.** Transcript assembly (StringTie2), classification (gffcompare), coding-potential consensus (CPAT + CPC2), catalog overlap tagging, expression recurrence filtering. Runs on all samples, condition-blind.
- **Phase B — Quantification.** Decoy-aware Salmon requantification against the frozen catalog, tximport, DESeq2 per norm_group, cis-regulatory associations, functional enrichment, per-contrast reports.

## Quick start

### 1. Build containers

```bash
docker build -t ghcr.io/damouzo/lncrna-nextflow/discovery:latest containers/discovery/
docker build -t ghcr.io/damouzo/lncrna-nextflow/coding_potential:latest containers/coding_potential/
docker build -t ghcr.io/damouzo/lncrna-nextflow/r_bioc:latest containers/r_bioc/
docker build -t ghcr.io/damouzo/lncrna-nextflow/reporting:latest containers/reporting/
```

GitHub Actions builds and pushes these to GHCR on every push to `main`.

### 2. Input files

**samplesheet.csv** — one row per lane, same format as nf-core/rnaseq. Multi-lane samples repeat the sample ID; lanes are grouped automatically.

```csv
sample,condition,batch,norm_group,fastq_1,fastq_2,bam
sample1,control,batch1,NB4,/path/lane1_R1.fq.gz,/path/lane1_R2.fq.gz,/path/sample1.bam
sample1,control,batch1,NB4,/path/lane2_R1.fq.gz,/path/lane2_R2.fq.gz,/path/sample1.bam
sample2,treated,batch1,NB4,/path/lane1_R1.fq.gz,/path/lane1_R2.fq.gz,/path/sample2.bam
```

**comparisons.csv** — one row per contrast.

```csv
contrast_name,numerator,denominator
treated_vs_control,treated,control
```

### 3. params.yaml

```yaml
input:       samplesheet.csv
contrasts:   comparisons.csv
outdir:      results
genome:      /path/to/genome.fa
gtf:         /path/to/annotation.gtf
design_formula: '~ batch + condition'

# Single strandedness contract for Salmon, featureCounts and StringTie.
# reverse = fr-firststrand (dUTP), the usual RNA-seq case.
library_strandedness: reverse

# Optional, per norm_group, default off:
heterogeneity_analysis:
  MNC: true
```

`design_formula` is a base formula — DESeq2 simplifies to `~ condition` within any group where `batch` is constant.

`salmon_libtype` and `featurecounts_strandness` still exist but only as optional
overrides of the derived values; if set they must be coherent with
`library_strandedness` or the pipeline fails. Two strandedness guards run before
and after quantification: `CHECK_STRANDEDNESS` compares the empirical RSeQC
orientation against the contract, and `CHECK_SALMON_LIBTYPE` verifies the actual
Salmon `expected_format` and `compatible_fragment_ratio` per sample.

A filled example lives in `assets/params.example.yaml`.

> Note: `catalog_export/*` is retired and no longer produced. The corrected
> FASTA/BED/GTF and gene/transcript tables will be emitted under `deliverables/catalog/`
> in the F4.3 deliverable layer; until then `results/*/catalog_export/` in older runs
> is stale and must not be used.

### 4. Optional external resources (conservation & synteny)

The pipeline is **species-agnostic** and never downloads anything: you point it at
static reference files in your `params.yaml`, exactly like `genome`/`gtf`. When
`conservation_bigwig` and/or `synteny_chain_file` are set, an early
`CHECK_REFERENCE_COMPATIBILITY` guard verifies they share the reference assembly
(fails fast on a build mismatch — the main silent-bug risk). Both are
**reporting-only**, never filters.

```yaml
# REQUIRED when enabling conservation or synteny:
genome_build: 'GRCh38.p14'

# Sequence conservation (phastCons / phyloP / GERP bigWig, same build as genome)
conservation_bigwig: /path/to/reference/conservation/hg38.phastCons100way.bw

# Synteny (liftOver chain + target-species GTF)
synteny_chain_file:     /path/to/reference/liftover/hg38ToMm39.over.chain.gz
synteny_target_species: 'mouse'
synteny_target_gtf:     /path/to/reference/mouse/Mus_musculus.GRCm39.gtf
```

Where to get each resource (reference only — the pipeline consumes paths):

| Resource | Where to obtain |
|----------|-----------------|
| phastCons/phyloP/GERP bigWig | UCSC Table Browser → *Conservation* track (select your build) |
| liftOver chain file | UCSC (`hgdownload.soe.ucsc.edu` /goldenPath/*sourceBuild*/liftOver/) or Ensembl |
| target-species GTF | Ensembl/Ensembl BioMart (`Mus_musculus.GRCm39.gtf`, etc.) |

### 4b. Stable catalog identity (novel lncRNA registry, optional)

Novel lncRNAs can receive **stable hierarchical IDs** (`LNC000001`, transcripts
`LNC000001.1`, …) assigned by a persistent two-level registry that survives across
runs. Set a registry path *outside* `results/` (must be **absolute** — the registry
lives on shared storage, and a relative path would resolve against the per-task
scratch work dir and break persistence):

```yaml
catalog_registry:         /path/to/outside/catalog/lncrna_catalog_registry.tsv
catalog_registry_project: MY_PROJECT        # provenance label
```

- The registry is read and extended, never regenerated; writes are atomic with a
  `history/` backup, guarded by an atomic `mkdir` lock, and fail fast if the
  build/annotation change. Transcripts are matched by structural signature; the
  gene id is inherited from its matched transcripts (span fallback only when no
  transcript matches). Any ambiguous match writes `registry_review_needed.tsv`
  and stops for manual resolution — never a silent reassignment.
- The registry header records `catalog_registry_project`: a run from a *different*
  project label fails fast unless you explicitly set
  `catalog_registry_allow_cross_project: true` (provenance/merge escape hatch).
- `ASSIGN_CATALOG_IDS` sets `cache false` as a deliberate `-resume` exception; the
  GTF/FASTA/BED of quantification keep the original `MSTRG.*` IDs.
- Omit `catalog_registry` to disable the feature (columns come out `NA`); no
  existing run changes behaviour.

### 4c. External lncRNA catalogs (optional — supports the novelty claim)

Without catalogs the run still completes but records `novelty_status = not_evaluated`
(never "clean"). To back the novelty claim, download **NONCODE** (hg38) yourself and
pass the local BED as a param — **the pipeline never downloads anything**, exactly
like `genome`/`gtf` and the §4 resources. Any subset of the catalogs below works; the
run uses whichever ones you configured.

| Source | Status | Where to obtain |
|--------|--------|-----------------|
| NONCODE v6 (hg38) | recommended (primary) | <http://www.noncode.org/download.php> → BED section, `NONCODEv6_hg38.lncAndGene.bed.gz` |
| RNAcentral (GRCh38) | recommended (second) | <https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/> → `homo_sapiens.GRCh38.bed.gz` |
| LNCipedia | optional / legacy | lncipedia.org is inactive; only from a verifiable mirror. Non-commercial licence |

**NONCODE** ships a ready BED. The host blocks requests without a browser User-Agent,
so pass one (the URL works over plain HTTP):

```bash
http://www.noncode.org/datadownload/NONCODEv6_hg38.lncAndGene.bed.gz
gunzip -k NONCODEv6_hg38.lncAndGene.bed.gz
md5sum NONCODEv6_hg38.lncAndGene.bed.gz  
head -3 NONCODEv6_hg38.lncAndGene.bed 
```

**RNAcentral** ships BED12 with exon blocks, but it is **all ncRNA** (rRNA, tRNA, snoRNA,
miRNA, lncRNA, …), not lncRNA-only — a small ncRNA overlapping ≥50 % of a candidate
still counts as a strict hit, so treat it as broad "known ncRNA" evidence:

```bash
wget -O homo_sapiens.GRCh38.bed.gz \
  https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz
gunzip -k homo_sapiens.GRCh38.bed.gz
```

```yaml
noncode_bed:    /path/to/reference/catalogs/NONCODEv6_hg38.lncAndGene.bed
rnacentral_bed: /path/to/reference/catalogs/homo_sapiens.GRCh38.bed
```

The overlap is **reciprocal and pair-by-pair** (threshold `external_catalog_min_reciprocal`,
default 0.5), strand anchored to the candidate. Only the strict call feeds
`known_catalog_hit` and the novel figures; the per-source coverage (`<src>_pct`) is
informative and never filters.

`CHECK_NOVELTY` fails fast if a configured BED has 0 features, shares no chromosome with
the annotation after normalization (a genuine assembly mismatch — `chr1`/`1` and
`chrM`/`MT` are already handled), yields zero hits, or the annotated lncRNAs have a 0
hit rate (positive control) — a bad download stops the run instead of silently marking
everything novel. It reports a run-wide `novelty_status` over the novel (MSTRG)
candidates (the guard's positive control is the non-MSTRG transcripts that enter the
candidate set, not the whole annotated catalog):

- `novel_for_all_catalogs` — no novel candidate matches any configured catalog;
- `known` — at least one novel candidate matches a catalog (check `known_catalog_hit` /
  `catalog_sources` before claiming novelty);
- `not_evaluated` — no catalog configured (never "clean").

The status lands in `run_info/novelty_status.tsv` (with `n_candidate_*` counts, scoped to
the candidate set, not the frozen catalog) and `run_info/qa_report.tsv`; the
deliverable-level novelty call is per gene in `catalog/catalog_genes.tsv`
(`novelty_class` / `known_catalog_hit` / `catalog_sources`).


### 5. Run

```bash
nextflow run /path/to/lncRNA-nextflow \
    -profile apocrita,apptainer \
    -params-file params.yaml \
    -w /gpfs/scratch/$USER/lncrna_work
```
### 6. Design decisions

- **norm_group separation.** Phase B runs independently per norm_group (NB4, MSCline, MNC, MSC) so cell lines and primary samples never share a size-factor normalization. Phase A uses all samples.
- **tx2gene from the frozen GTF.** `ANNOTATION_FREEZE` parses the frozen GTF directly, so transcript IDs always match what Salmon sees. `gene_id` is anchored to `; ` so `ref_gene_id` is never captured by mistake.
- **Decoy-aware index.** The Salmon gentrome includes the full reference transcriptome, novel lncRNAs, and genome decoys.
- **Batch handling.** DESeq2 falls back to `~ condition` automatically when `batch` is constant within a group.
- **No hard filters on class codes.** Antisense and monoexonic transcripts are flagged as risk, not removed. You decide what to exclude.
- **Heterogeneity module (optional).** Per norm_group, opt-in via `params.heterogeneity_analysis`. Complements DESeq2 with a per-sample concordance check, useful for heterogeneous primary-sample groups — see Carrasco-Leon et al. 2021 for the method. Off by default; does not affect DESeq2/cis-associations/enrichment when disabled.