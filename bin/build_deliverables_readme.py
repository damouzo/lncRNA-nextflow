#!/usr/bin/env python3
"""F4.7: generate the deliverables README and its column dictionary from the contract.

The README is derived, never hand-edited: `docs/deliverables_contract.json` is the
single source of truth for the file tree and the required columns. Stdlib only.

Writes two files (names chosen so the Nextflow publishDir can rename them):
  --deliverables-readme   -> deliverables/README.md
  --run-info-readme       -> deliverables/run_info/README.md

CLI:
  build_deliverables_readme.py --contract docs/deliverables_contract.json --outdir . \
      --deliverables-readme __deliverables_README.md \
      --run-info-readme __run_info_README.md [--qa-report qa_report.tsv] \
      [--publish-level core_raw]
"""

import argparse
import csv
import json
import os


def load_contract(path):
    with open(path) as fh:
        return json.load(fh)


def read_kv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        return [row for row in rdr if row]


def _fmt(value):
    if value is None or value == "":
        return "NA"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    return ("%g" % f)


# Meanings for the columns that recur across the deliverable tables. Kept here
# rather than in the contract because the contract only carries name/type.
COLUMN_MEANINGS = {
    "gene_id": "Stable gene id (Ensembl, or catalog `LNC` id for novel genes).",
    "gene_id_stable": "One stable id per gene for downstream joins: the catalog `LNC` id for "
                      "novel genes, the reference `gene_id` otherwise.",
    "gene_name": "Gene symbol; `NA` for novel genes without a symbol.",
    "transcript_id": "Transcript id (Ensembl or catalog transcript id). A `__ENST...` "
                     "suffix (e.g. `MSTRG.351__ENST00000584329`) means StringTie reused a "
                     "reference transcript id for a partial match: it names the closest "
                     "reference transcript, not the lncRNA — use `catalog_transcript_id`.",
    "origin": "`novel` (StringTie/GFFCompare) or `annotated_lncrna`.",
    "gene_group": "Coarse group from `gene_biotype`: annotated_lncRNA / novel_lncRNA / "
                  "protein_coding / other. Use this (not `origin`) to separate coding from "
                  "non-coding reference genes.",
    "chr": "Chromosome / contig of the gene span (from the catalog).",
    "start": "Gene span start (1-based, GRCh38).",
    "end": "Gene span end (1-based, GRCh38).",
    "strand": "Gene span strand (+/-).",
    "genome_build": "Genome build the coordinates refer to (e.g. GRCh38, Ensembl 115).",
    "contig_type": "`canonical` for GRCh38 1-22/X/Y/MT, else `noncanonical` "
                   "(patches, scaffolds, alt contigs). A non-canonical locus with a high "
                   "baseMean is a repeat/rRNA suspect.",
    "baseMean_extreme": "TRUE when baseMean >= the configured `extreme_basemean`. A flag, "
                        "never a filter.",
    "significant_padj": "padj < FDR. This is what 'significant' means by default.",
    "significant_padj_lfc1": "padj < FDR AND |log2FoldChange| >= LFC. Used only where a "
                             "table/panel names it explicitly (e.g. the report Overview tab).",
    "spearman_pvalue": "cis correlation p-value (t-approximation on `df_effective`; `NA` when "
                       "not computable, never floored). `p_method` records how it was obtained.",
    "df_effective": "Degrees of freedom used for the cis p-value (n samples minus one per "
                    "condition mean for the residual design).",
    "p_method": "How the cis p-value was computed (e.g. t_approx_residual_df), or why it is NA.",
    "loo_min_abs_rho": "cis only: weakest |rho| when each sample is dropped in turn "
                       "(leave-one-out). Near 0 means the pair leans on a single sample.",
    "loo_sign_stable": "cis only: TRUE if the correlation sign survives dropping any one "
                       "sample. Prefer TRUE candidates.",
    "loo_max_delta": "cis only: largest |rho| change when dropping one sample.",
    "lncrna_origin": "lncRNA origin: novel (discovered) or annotated_lncrna (reference).",
    "lncrna_class_code": "GFFCOMPARE class code of the lncRNA (novel only; NA for annotated "
                         "lncRNAs and reference genes).",
    "lncrna_position_class": "Position of the lncRNA vs protein-coding genes: antisense / "
                             "intronic / overlapping / no_coding_overlap (defined for every "
                             "lncRNA, novel or annotated).",
    "class_code": "GFFCompare code: u/i/x/j/o = novel, =/c/k/m/n = reference.",
    "tier_static": "Categorical tier from class code + parent biotype + exon count.",
    "tier_final": "Tier after flags; equals tier_static in this version.",
    "single_exon_low_confidence": "TRUE when a mono-exonic candidate is demoted to Tier 2.",
    "antisense_status": "Antisense-leakage flag; annotates, never filters.",
    "baseMean": "Mean normalised count across the samples of the contrast.",
    "log2FoldChange": "Shrunken LFC (lfcShrink/apeglm); the sign is the direction.",
    "lfcSE": "Standard error of the shrunken LFC.",
    "pvalue": "Wald test p-value, unadjusted.",
    "padj": "BH-adjusted p-value; `NA` = independent filtering or baseMean 0.",
    "contrast": "`<numerator>_vs_<denominator>` label of the comparison.",
    "batch_subset": "TRUE when the contrast used a batch-restricted sample subset.",
    "sample": "Sample id from the samplesheet.",
    "norm_group": "Normalisation/comparison group.",
    "condition": "Experimental condition used by the design formula.",
    "batch": "Batch label.",
    "keep": "TRUE when the gene passed the condition-blind expression filter.",
    "filter_reason": "Why a gene was dropped (`low_expression`, ...); empty when kept.",
    "step": "Funnel stage name (catalog, coding consensus, DE universe, ...).",
    "n": "Count for the funnel step (1 on per-gene rows).",
    "n_universe": "Rows in `de_full.tsv` (all quantified genes, tested or not).",
    "n_tested": "Genes DESeq2 actually tested (independent filtering removes the rest).",
    "n_significant": "Rows in `de_full_significant.tsv` (`significant_padj`, padj < FDR).",
    "n_lncrna_full": "Rows in `de_lncrna_full.tsv` (all lncRNAs, significant or not).",
    "n_lncrna_significant": "Rows in `de_lncrna_significant.tsv` (significant lncRNAs).",
    "reason": "Motive attached to a funnel step or a dropped gene.",
    "lnc_position_class": "Position vs protein-coding genes: antisense / intronic / "
                          "overlapping / no_coding_overlap. Only protein-coding genes "
                          "are compared, so an overlap with an annotated lncRNA or "
                          "pseudogene is still no_coding_overlap (not 'intergenic').",
    "gene_biotype": "Reference gene biotype (e.g. protein_coding, lncRNA, pseudogene).",
    "gene_group": "Coarse category used across the deliverables: protein_coding / "
                  "annotated_lncRNA / novel_lncRNA / other (other = reference genes that "
                  "are neither protein-coding nor annotated lncRNAs: pseudogenes, TEC, TR, ...).",
    "catalog_id": "Stable catalog id (gene level) assigned by the registry.",
    "catalog_transcript_id": "Stable catalog transcript id (e.g. `LNC000001.1`), the id "
                             "used in the annotation GTF and FASTA.",
    "mstrg_transcript_id": "Discovery transcript id (`MSTRG.x.y`, or the disambiguated "
                           "`MSTRG.x__ENST...` form) the stable id was assigned to.",
    "closest_reference_transcript": "gffcompare `cmp_ref`: the reference transcript the "
                                    "novel transcript overlaps most (its closest match).",
    "shrink_type": "Which lfcShrink ran for the contrast (apeglm, or ashr fallback).",
    "stat": "Unshrunken Wald test statistic (not used for the significant flag).",
    "lfc_unshrunken": "Unshrunken log2 fold change (the shrunken value is log2FoldChange).",
    "tier": "Categorical tier: 1a / 1b / 2 for novel genes; NA for reference rows. "
            "A label, never a filter.",
    "overlaps_target": "TRUE when the lncRNA and target spans overlap (putative host).",
    "distance_bp": "Gap between the lncRNA and target spans (0 when they overlap).",
    "correlation_input": "What the cis correlation used: residuals_intra_condition "
                         "(within-condition residuals) or expression.",
    "spearman_rho": "Spearman correlation (on intra-condition residuals by default).",
}


def build_column_dictionary(contract):
    """Column name -> type + meaning, from the contract columns plus the glossary."""
    seen = {}
    for spec in contract.get("files", []):
        for col in spec.get("columns", []):
            seen.setdefault(col["name"], col.get("type", "string"))
    lines = ["## Column dictionary", "",
             "Columns that recur across the deliverable tables. `NA` means "
             "*not evaluated*, never *clean*.", "",
             "| column | type | meaning |", "|---|---|---|"]
    for name, typ in seen.items():
        lines.append("| `%s` | %s | %s |" % (name, typ, COLUMN_MEANINGS.get(name, "")))
    lines.append("")
    # Columns documented here but not typed by the contract. The cis and funnel
    # tables can legitimately be empty, so the contract declares no columns for
    # them; without this block their meaning would never render.
    extra = [name for name in COLUMN_MEANINGS if name not in seen]
    if extra:
        lines.append("Columns that recur in the published tables but are not "
                     "declared in the contract:")
        lines.append("")
        for name in sorted(extra):
            lines.append("- `%s` — %s" % (name, COLUMN_MEANINGS[name]))
        lines.append("")
    return lines


def _start_here():
    """Question -> file guide (part a of the deliverables README)."""
    lines = []
    lines.append("## Start here")
    lines.append("")
    lines.append("| I want to… | Open |")
    lines.append("|---|---|")
    lines.append("| see the novel lncRNAs discovered | `catalog/novel_lncrnas.tsv` |")
    lines.append("| see every gene / transcript in the catalog | `catalog/catalog_genes.tsv`, `catalog/catalog_tx.tsv` |")
    lines.append("| see the candidates kept/filtered in discovery | `catalog/candidates_all.tsv` |")
    lines.append("| compare DE counts across contrasts | `contrasts/contrasts_overview.tsv` |")
    lines.append("| see what changed in one contrast | `contrasts/<group>__<contrast>/de_full.tsv` |")
    lines.append("| see only the significant genes | `contrasts/<group>__<contrast>/de_full_significant.tsv` |")
    lines.append("| get a complete lncRNA-only table for Excel | `contrasts/<group>__<contrast>/de_lncrna_full.tsv` |")
    lines.append("| see only the significant lncRNAs | `contrasts/<group>__<contrast>/de_lncrna_significant.tsv` |")
    lines.append("| read the narrative report of a contrast | `contrasts/<group>__<contrast>/<id>_report.html` |")
    lines.append("| see functional enrichment for a contrast | `contrasts/<group>__<contrast>/enrichment_ora.tsv`, `enrichment_gsea.tsv` |")
    lines.append("| check cis candidates of a group | `groups/<group>/cis_pairs_significant.tsv` |")
    lines.append("| get the expression matrix of a group | `groups/<group>/counts.tsv`, `groups/<group>/tpm.tsv` |")
    lines.append("| see run metadata and tool versions | `run_info/` |")
    lines.append("")
    lines.append("Each folder is computed at one level, so it never duplicates work:")
    lines.append("")
    lines.append("| Folder | Computed once per | Per contrast? |")
    lines.append("|---|---|---|")
    lines.append("| `catalog/` | project | no |")
    lines.append("| `groups/<group>/` | group (shared matrix) | no |")
    lines.append("| `contrasts/<id>/` | (group, contrast) | yes |")
    lines.append("")
    return lines


def _tiers_section():
    """How to read the tiers, including the two wet-lab cautions."""
    lines = []
    lines.append("## How to read the tiers")
    lines.append("")
    lines.append("- A tier is a **structural priority label**, not a confidence score: "
                 "`1a` = class `u`/`x`/`i`, `1b` = class `j`/`o` whose reference host is an "
                 "annotated lncRNA, `2` = mono-exonic or other overlaps. `i` (intronic) and "
                 "`x` (antisense) do overlap reference genes; `u`/`x`/`i` are grouped because "
                 "they cannot be a simple isoform of an annotated gene, not because they are "
                 "cleaner.")
    lines.append("- Tier 1 is **not** confirmed. Class `i` on the same strand as its host is "
                 "often nascent pre-mRNA in total RNA without poly-A selection, and class `x` "
                 "depends on strand fidelity (see the strandedness fidelity in `run_info/`). "
                 "Treat Tier 1 as a candidate set to validate.")
    lines.append("- `catalog/tier_final_by_class.tsv`: which gffcompare classes compose Tier 1 "
                 "(transcript and gene level, with and without the mono-exonic rule).")
    lines.append("- `catalog/validation/tier_final_sensitivity.tsv`: audit only. The Tier 1 "
                 "gene count under alternative flank/antisense thresholds; it does not move "
                 "because, after decision D1, those thresholds no longer gate the tier.")
    lines.append("")
    return lines


def _where_to_look():
    """One-line answer for each common question (layers: catalog / groups / contrasts)."""
    lines = []
    lines.append("## Where to look")
    lines.append("")
    lines.append("- Identity, sequence and tier of a novel lncRNA: `catalog/novel_lncrnas.tsv` "
                 "+ `catalog/annotation/lncrna_catalog.fa`.")
    lines.append("- Which candidates were kept or filtered in discovery: "
                 "`catalog/candidates_all.tsv`.")
    lines.append("- Expression matrix and cis candidates for one model/group: "
                 "`groups/<group>/tpm.tsv`, `groups/<group>/cis_pairs_significant.tsv`.")
    lines.append("- What changes between conditions: `contrasts/<id>/de_full_significant.tsv` (any gene), or "
                 "`de_lncrna_full.tsv` / `de_lncrna_significant.tsv` for the lncRNA-only view.")
    lines.append("- Run metadata, thresholds and tool versions: `run_info/`.")
    lines.append("")
    return lines


def build_deliverables_readme(contract, qa_rows, publish_level, fdr=None, lfc=None):
    lines = []
    lines.append("# lncRNA discovery — deliverables")
    lines.append("")
    lines.append("This README and the column dictionary are generated from "
                 "`docs/deliverables_contract.json`; do not edit by hand.")
    lines.append("")
    lines.append("Publish level: **%s**." % publish_level)
    lines.append("")
    lines.extend(_start_here())
    lines.extend(_where_to_look())
    lines.extend(build_column_dictionary(contract))
    lines.append("## Significance")
    lines.append("")
    lines.append("- 'significant' means `significant_padj` = `padj < %s`. "
                 "`significant_padj_lfc1` adds `|log2FoldChange| >= %s` and is used only "
                 "where a table/panel names it explicitly." % (_fmt(fdr), _fmt(lfc)))
    lines.append("- The LFC cutoff is applied only to `significant_padj_lfc1`, on the **shrunken** "
                 "`log2FoldChange` from `lfcShrink`; it is not passed to DESeq2 "
                 "`results()`, so `pvalue`/`padj` keep their standard meaning.")
    lines.append("- Defaults live in `nextflow.config`/`nextflow_schema.json` "
                 "(`fdr_threshold`, `lfc_threshold`) and are recorded in "
                 "`run_info/params_used.tsv`.")
    lines.append("")
    lines.append("## Opening the TSVs in Excel")
    lines.append("")
    lines.append("Do **not** double-click the `.tsv` files. Use **Data → From Text/CSV** "
                 "(or *Importar desde texto*) and set the columns holding gene symbols "
                 "to **Text**: Excel silently turns `SEPT2` into a date and `MARCH1` into "
                 "`1-Mar`. The same applies to ids that look numeric.")
    lines.append("")
    lines.append("## Honesty note")
    lines.append("")
    lines.append("- Tiers are **labels, not filters**. No table drops a gene by tier or "
                 "significance; every gene keeps its marks (`single_exon_low_confidence`, "
                 "`antisense_status`, `x_host_leak_suspect`).")
    lines.append("- `tier_final` equals `tier_static` in this version (the antisense flag "
                 "annotates but does not demote).")
    lines.append("- `NA` means *not evaluated*, never *clean*.")
    lines.append("")
    lines.extend(_tiers_section())
    lines.append("## What \"novel\" means here")
    lines.append("")
    lines.append("- `novel` = a transcript not present in the reference annotation "
                 "(Ensembl 115). It is **not** a claim that the locus was never described "
                 "before: no external catalog was configured for this run, so novelty is "
                 "*not evaluated* against other databases.")
    lines.append("- Coding genes (and other non-lncRNA reference genes) are tested and shipped "
                 "in full in `de_full.tsv` under `gene_group = protein_coding` / `other`. "
                 "`de_lncrna_full.tsv` contains **only** lncRNAs (annotated + novel).")
    lines.append("- `reference_fraction_recovered` (in `catalog/validation/recovery_reference.tsv`) "
                 "is the fraction of reference lncRNA genes detected at >= 1 TPM in >= 3 of the "
                 "samples — a low value is expected for this definition, not a failure.")
    lines.append("- A disambiguated transcript id like `MSTRG.351__ENST00000584329` means StringTie "
                 "reused a reference transcript id for a partial match: the `__ENST` part names the "
                 "reference transcript it overlaps, it is not the id of the lncRNA (`catalog_transcript_id`).")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    for spec in contract.get("files", []):
        cols = ", ".join(c["name"] for c in spec.get("columns", []))
        lines.append("### `%s`" % spec["path"])
        lines.append("")
        lines.append("- kind: %s" % spec.get("kind", "table"))
        lines.append("- min rows: %s" % spec.get("min_rows", 0))
        if cols:
            lines.append("- required columns: %s" % cols)
        lines.append("")
    if qa_rows:
        lines.append("## Run summary")
        lines.append("")
        lines.append("| metric | value |")
        lines.append("|---|---|")
        for row in qa_rows:
            if len(row) >= 2:
                lines.append("| %s | %s |" % (row[0], row[1]))
        lines.append("")
    lines.append("## Provenance")
    lines.append("")
    lines.append("Run metadata, thresholds, strandedness and external catalogs live under "
                 "`../run_info/` (audit) and `run_info/` (summary). The execution trace, "
                 "report, timeline and DAG are under `run_info/pipeline_info/`.")
    lines.append("")
    return "\n".join(lines)


def build_run_info_readme(contract, qa_rows, publish_level):
    lines = []
    lines.append("# run_info — summary")
    lines.append("")
    lines.append("Curated metadata for this run (publish level **%s**). Generated from "
                 "`docs/deliverables_contract.json`." % publish_level)
    lines.append("")
    lines.append("| file | content |")
    lines.append("|---|---|")
    lines.append("| `versions.tsv` | pipeline and tool versions (best effort) |")
    lines.append("| `sample_index.tsv` | sample, norm_group, condition, batch |")
    lines.append("| `qa_report.tsv` | headline QC; `NA` = not evaluated |")
    lines.append("| `pipeline_info/` | Nextflow trace, report, timeline, DAG |")
    lines.append("")
    lines.append("There are two run_info trees by design: this `deliverables/run_info/` holds the "
                 "curated summary (versions, sample index, QA, pipeline_info); the top-level "
                 "`<outdir>/run_info/` holds the raw audit set (params_used, thresholds, "
                 "strandedness, external catalogs, manifest, funnel, qc_controls).")
    lines.append("")
    if qa_rows:
        lines.append("| metric | value |")
        lines.append("|---|---|")
        for row in qa_rows:
            if len(row) >= 2:
                lines.append("| %s | %s |" % (row[0], row[1]))
        lines.append("")
    lines.append("The full audit set (parameters, thresholds, strandedness summary, "
                 "external catalogs, manifest) lives at `<outdir>/run_info/`.")
    lines.append("")
    return "\n".join(lines)


def build_raw_readme(publish_level):
    lines = []
    lines.append("# raw/ — native outputs (audit)")
    lines.append("")
    lines.append("Native, function-named outputs kept for reproducibility (publish level "
                 "**%s**). Not the deliverable: use `../deliverables/` for the shippable "
                 "contract and `../run_info/` for run metadata." % publish_level)
    lines.append("")
    lines.append("| Folder | Content |")
    lines.append("|---|---|")
    lines.append("| `deseq2/` | DESeq2 results per group (`*_DE.tsv`, `dds.rds`) |")
    lines.append("| `tximport_filter/` | tximport counts/TPM, filter stats, pre-filter universe |")
    lines.append("| `enrichment/` | group-level ORA/GSEA tables (repartitioned per contrast in deliverables) |")
    lines.append("| `cis_associations/` | cis pair tables per group |")
    lines.append("| `stringtie/`, `stringtie_merge/`, `gffcompare/` | assembly and merge |")
    lines.append("| `salmon_quant/`, `decoy_index/` | quantification |")
    lines.append("| `heterogeneity/` | opt-in heterogeneity diagnostic (also embedded in the reports) |")
    lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--deliverables-readme", default="__deliverables_README.md")
    ap.add_argument("--run-info-readme", default="__run_info_README.md")
    ap.add_argument("--raw-readme", default="__raw_README.md")
    ap.add_argument("--qa-report", default=None)
    ap.add_argument("--publish-level", default="core_raw")
    ap.add_argument("--fdr", default=None)
    ap.add_argument("--lfc", default=None)
    args = ap.parse_args()

    contract = load_contract(args.contract)
    qa_rows = read_kv(args.qa_report)
    os.makedirs(args.outdir, exist_ok=True)

    with open(os.path.join(args.outdir, args.deliverables_readme), "w") as fh:
        fh.write(build_deliverables_readme(contract, qa_rows, args.publish_level,
                                           fdr=args.fdr, lfc=args.lfc))
    with open(os.path.join(args.outdir, args.run_info_readme), "w") as fh:
        fh.write(build_run_info_readme(contract, qa_rows, args.publish_level))
    with open(os.path.join(args.outdir, args.raw_readme), "w") as fh:
        fh.write(build_raw_readme(args.publish_level))
    print("build_deliverables_readme: wrote %s, %s and %s"
          % (args.deliverables_readme, args.run_info_readme, args.raw_readme))


if __name__ == "__main__":
    main()
