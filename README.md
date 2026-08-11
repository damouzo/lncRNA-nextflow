# lncRNA-nextflow

lncRNA discovery, curation, and differential expression pipeline.

Two phases, strictly separated to avoid selection bias:

- **Phase A — Discovery**: transcript assembly (StringTie2), classification (gffcompare), coding-potential consensus (CPAT + CPC2), catalog overlap tagging, expression recurrence filtering. Condition-blind — no contrast influences which transcripts are kept.
- **Phase B — Quantification**: decoy-aware Salmon requantification against the frozen catalog, tximport, DESeq2 per norm_group, cis-regulatory associations, functional enrichment.

## Quick start

### 1. Build containers

```bash
# Push to main — GitHub Actions builds + pushes to GHCR:
#   ghcr.io/<org>/lncrna-nextflow/discovery:latest
#   ghcr.io/<org>/lncrna-nextflow/coding_potential:latest
#   ghcr.io/<org>/lncrna-nextflow/r_bioc:latest

# Or build locally:
docker build -t ghcr.io/damouzo/lncrna-nextflow/discovery:latest containers/discovery/
docker build -t ghcr.io/damouzo/lncrna-nextflow/coding_potential:latest containers/coding_potential/
docker build -t ghcr.io/damouzo/lncrna-nextflow/r_bioc:latest containers/r_bioc/
```

### 2. Prepare input files

**samplesheet.csv** (one row per lane — same format as nf-core/rnaseq)
```csv
sample,condition,batch,norm_group,fastq_1,fastq_2,bam
sample1,control,batch1,NB4,/path/lane1_R1.fq.gz,/path/lane1_R2.fq.gz,/path/sample1.bam
sample1,control,batch1,NB4,/path/lane2_R1.fq.gz,/path/lane2_R2.fq.gz,/path/sample1.bam
sample2,treated,batch1,NB4,/path/lane1_R1.fq.gz,/path/lane1_R2.fq.gz,/path/sample2.bam
```

Multi-lane samples get multiple rows with the same sample ID — lanes are grouped automatically. No need to concatenate FASTQs beforehand; Salmon handles multiple input files natively.

**comparisons.csv**
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
# Used as base formula; DESeq2 simplifies to ~ condition automatically
# within groups where batch is constant (see Design decisions below).
salmon_libtype: 'SR'
```

### 4. Run

```bash
nextflow run /path/to/lncRNA-nextflow \
    -profile apocrita,apptainer \
    -params-file params.yaml \
    -w /gpfs/scratch/$USER/lncrna_work
```


## Output

```
results/
├── discovery/
│   ├── stringtie/           per-sample assemblies
│   ├── stringtie_merge/     merged GTF
│   ├── gffcompare/          class codes
│   ├── cpat/ cpc2/          coding-potential scores
│   ├── coding_consensus/    CPAT+CPC2 agreement table
│   ├── expression_recurrence/  recurrence-filtered candidates
│   └── annotation_freeze/   frozen GTF + tx2gene_detailed.tsv
├── quantification/
│   ├── decoy_index/         Salmon index
│   ├── salmon_quant/        per-sample quantification dirs
│   ├── per_group_coldata/   coldata CSVs split by norm_group
│   ├── tximport_filter/     per-group filtered counts + PCA
│   ├── deseq2/              per-group DE results
│   ├── cis_associations/    lncRNA–gene pairs
│   └── enrichment/          ORA + GSEA tables and plots
└── pipeline_info/           timeline, trace, DAG
```

## Design decisions

- **norm_group separation**: Phase B runs independently per norm_group (e.g., NB4, MSCline, MNC, MSC). This avoids mixing cell lines with primary samples in the same size-factor normalisation. Phase A (discovery) uses all samples — it's condition-blind.
- **Batch handling**: DESeq2 auto-detects whether the batch column has multiple levels. If batch is constant within a group, it falls back to `~ condition`.
- **No hard filters on class codes**: antisense ("x") and monoexonic transcripts are flagged as risk, not deleted. You decide what to exclude.

## What's not done yet

- Catalog overlap module runs but has no real reference BEDs wired up yet
- Monoexonic intron-coverage and antisense-neighbor risk flags are stubs
- No test suite (`nf-test` planned)
- Enrichment module currently runs ORA/GSEA on all DE genes, not just cis targets
