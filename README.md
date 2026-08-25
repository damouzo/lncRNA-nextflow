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
salmon_libtype: 'SR'

# Optional, per norm_group, default off:
heterogeneity_analysis:
  MNC: true
```

`design_formula` is a base formula — DESeq2 simplifies to `~ condition` within any group where `batch` is constant.

A filled example lives in `assets/params.example.yaml`.

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