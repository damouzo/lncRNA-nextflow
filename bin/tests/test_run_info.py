#!/usr/bin/env python3
"""F4.6 / F4.7 (run_info + generated README) tests. Stdlib only.

Run: python3 bin/tests/test_run_info.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import build_run_info as R          # noqa: E402
import build_deliverables_readme as D  # noqa: E402
import validate_deliverables as V   # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def make_fixture(d):
    write(os.path.join(d, "samples.tsv"),
          "sample\tnorm_group\tcondition\tbatch\n"
          "S1\tMNC\tDDX41\tb1\nS2\tMNC\tHealthy\tb1\n")
    write(os.path.join(d, "strandedness_summary.tsv"),
          "sample\tdata_type\tdetected_orientation\tdeclared_strandedness\tfidelity\n"
          "S1\tpaired\treverse\treverse\t0.81\n"
          "S2\tpaired\treverse\treverse\t0.88\n")
    write(os.path.join(d, "thresholds_used.tsv"),
          "metric\tvalue\torigin\tn_controls\tseparation\tstatus\n"
          "antisense_ratio\t0.32\tderived\t897\t0.78\tthresholds_not_specific\n")
    write(os.path.join(d, "gene_catalog.tsv"),
          "gene_id\torigin\ttier_static\nG1\tnovel\tTRUE\nG2\tannotated_lncrna\tNA\n")
    write(os.path.join(d, "tier_final_sensitivity.tsv"),
          "changed\tpercentile\tn_tier1_final_genes\nbaseline\tderived\t2507\n")
    write(os.path.join(d, "novelty_status.tsv"),
          "metric\tvalue\nnovelty_status\tevaluated\nn_catalogs\t2\n")
    write(os.path.join(d, "versions_a.yml"),
          '"STRINGTIE2_MERGE":\n'
          '    tool: "stringtie"\n'
          '    version: "stringtie 3.0.0"\n'
          '    container: "quay.io/x:stringtie-3.0.0"\n'
          '"DESEQ2_DE":\n'
          '    tool: "DESeq2/apeglm"\n'
          '    version: "R 4.3.2; DESeq2 1.42.0"\n'
          '    container: "img/r_bioc:1"\n')


def test_run_info():
    d = tempfile.mkdtemp()
    make_fixture(d)
    out = os.path.join(d, "out")
    old = sys.argv
    sys.argv = ["build_run_info.py", "--outdir", out,
                "--samples", os.path.join(d, "samples.tsv"),
                "--strandedness", os.path.join(d, "strandedness_summary.tsv"),
                "--thresholds", os.path.join(d, "thresholds_used.tsv"),
                "--gene-catalog", os.path.join(d, "gene_catalog.tsv"),
                "--sensitivity", os.path.join(d, "tier_final_sensitivity.tsv"),
                "--versions-yml", os.path.join(d, "versions_a.yml"),
                "--novelty-status", os.path.join(d, "novelty_status.tsv"),
                "--pipeline-version", "0.1.0", "--publish-level", "core_raw",
                "--outdir-param", "results",
                "--catalog", "lncipedia=/does/not/exist.bed"]
    try:
        R.main()
        ok = True
    finally:
        sys.argv = old
    check(ok, "build_run_info CLI runs")
    for f in ("params_used.tsv", "versions.tsv", "sample_index.tsv",
              "run_manifest.tsv", "external_catalogs.tsv", "qa_report.tsv"):
        check(os.path.exists(os.path.join(out, f)), "wrote %s" % f)

    _, qa = R.read_tsv(os.path.join(out, "qa_report.tsv"))
    q = {r["metric"]: r["value"] for r in qa}
    check(q["n_samples"] == "2", "qa n_samples")
    check(q["n_tier1_static_genes"] == "1", "qa n_tier1_static_genes")
    check(q["n_tier1_final_genes"] == "2507", "qa n_tier1_final_genes from sensitivity baseline")
    check(q["strandedness_fidelity_min"] == "0.8100" and q["strandedness_fidelity_max"] == "0.8800",
          "qa fidelity range")
    check(q["thresholds_status"] == "thresholds_not_specific", "qa thresholds status")
    check(q["novelty_status"] == "evaluated", "qa novelty_status comes from CHECK_NOVELTY")

    _, ver = R.read_tsv(os.path.join(out, "versions.tsv"))
    check(len(ver) == 2 and ver[0]["process"] == "STRINGTIE2_MERGE",
          "versions.tsv has one row per emitting process")
    check("stringtie 3.0.0" in ver[0]["version"] and
          ver[0]["container"] == "quay.io/x:stringtie-3.0.0",
          "versions.tsv keeps version + container (colon preserved)")
    rows = R.parse_versions_yml([os.path.join(d, "versions_a.yml")])
    check(rows[1][0] == "DESEQ2_DE" and "DESeq2" in rows[1][2],
          "parse_versions_yml reads multi-block files")

    _, cats = R.read_tsv(os.path.join(out, "external_catalogs.tsv"))
    check(cats[0]["name"] == "lncipedia" and cats[0]["present"] == "FALSE",
          "external catalog absent -> present FALSE")


def test_versions_dedup():
    d = tempfile.mkdtemp()
    block = ('"QUANTIFICATION:DESEQ2_DE":\n'
             '    tool: "DESeq2/apeglm"\n'
             '    version: "R 4.3.2; DESeq2 1.42.0"\n'
             '    container: "img/r_bioc:1"\n')
    for i in (1, 2, 3):
        write(os.path.join(d, "vdir", "%d.versions.yml" % i), block)
    write(os.path.join(d, "vdir", "4.versions.yml"),
          '"SALMON_QUANTIFY":\n'
          '    tool: "salmon"\n'
          '    version: "salmon 1.10.0"\n'
          '    container: "img/discovery:1"\n')
    rows = R.parse_versions_yml([os.path.join(d, "vdir", "*.versions.yml")])
    check(len(rows) == 2, "one row per process despite one file per task (%d)" % len(rows))
    check([r[0] for r in rows] == ["QUANTIFICATION:DESEQ2_DE", "SALMON_QUANTIFY"],
          "dedup keeps file order and the first occurrence")


def test_versions_yml_cli():
    """PUBLISH_RUN_INFO passes the flag once per file, not once with N values."""
    d = tempfile.mkdtemp()
    write(os.path.join(d, "a.yml"),
          '"PROCA":\n    tool: "a"\n    version: "1"\n    container: "x"\n')
    write(os.path.join(d, "b.yml"),
          '"PROCB":\n    tool: "b"\n    version: "2"\n    container: "y"\n')
    write(os.path.join(d, "s.tsv"), "sample\tnorm_group\tcondition\tbatch\n")
    out = os.path.join(d, "out")
    old = sys.argv
    sys.argv = ["build_run_info.py", "--outdir", out, "--samples", os.path.join(d, "s.tsv"),
                "--versions-yml", os.path.join(d, "a.yml"),
                "--versions-yml", os.path.join(d, "b.yml")]
    try:
        R.main()
        ok = True
    finally:
        sys.argv = old
    _, ver = R.read_tsv(os.path.join(out, "versions.tsv"))
    check(ok and [r["process"] for r in ver] == ["PROCA", "PROCB"],
          "repeated --versions-yml flags are all read (got %s)"
          % ([] if not ok else [r["process"] for r in ver]))


def test_readme():
    contract = {
        "contract": "t", "version": "1",
        "files": [
            {"path": "catalog/catalog_genes.tsv", "kind": "table", "min_rows": 1,
             "columns": [{"name": "gene_id", "type": "string"}]},
            {"path": "catalog/annotation/lncrna_catalog.fa", "kind": "fasta"},
        ],
    }
    text = D.build_deliverables_readme(contract, [["n_samples", "2"]], "core_raw",
                                       fdr=0.05, lfc=0.0)
    check("catalog/catalog_genes.tsv" in text, "README lists contract path")
    check("labels, not filters" in text, "README carries honesty note")
    check("| n_samples | 2 |" in text, "README includes run summary")
    check("padj < 0.05" in text and "|log2FoldChange| >= 0" in text,
          "README cites the fdr/lfc defaults")
    check("SEPT2" in text and "Importar desde texto" in text,
          "README warns about Excel text import")
    rtext = D.build_run_info_readme(contract, [["n_samples", "2"]], "core_raw")
    check("versions.tsv" in rtext and "run_info/" in rtext, "run_info README points to audit set")

    d = tempfile.mkdtemp()
    old = sys.argv
    sys.argv = ["build_deliverables_readme.py", "--contract", os.path.join(d, "c.json"),
                "--outdir", d, "--qa-report", os.path.join(d, "missing.tsv")]
    write(os.path.join(d, "c.json"), __import__("json").dumps(contract))
    try:
        D.main()
        ok = True
    finally:
        sys.argv = old
    check(ok and os.path.exists(os.path.join(d, "__deliverables_README.md")),
          "README CLI writes deliverables README")


def test_validator_kinds():
    contract = {
        "contract": "t", "version": "1",
        "files": [
            {"path": "catalog/x.tsv", "kind": "table", "min_rows": 1,
             "columns": [{"name": "gene_id", "type": "string"}]},
            {"path": "catalog/annotation/x.fa", "kind": "fasta"},
        ],
    }
    d = tempfile.mkdtemp()
    write(os.path.join(d, "catalog/x.tsv"), "gene_id\nG1\n")
    write(os.path.join(d, "catalog/annotation/x.fa"), ">t1\nACGT\n")
    check(V.validate(d, contract) == [], "validator passes table + fasta present")

    d2 = tempfile.mkdtemp()
    write(os.path.join(d2, "catalog/x.tsv"), "gene_id\nG1\n")
    errs = V.validate(d2, contract)
    check(any("no file matches" in e and "x.fa" in e for e in errs),
          "validator reports missing fasta")

    d3 = tempfile.mkdtemp()
    write(os.path.join(d3, "catalog/x.tsv"), "gene_id\nG1\n")
    write(os.path.join(d3, "catalog/annotation/x.fa"), "")
    errs = V.validate(d3, contract)
    check(any("empty" in e for e in errs), "validator reports empty fasta")


def main():
    for fn in (test_run_info, test_versions_dedup, test_versions_yml_cli,
               test_readme, test_validator_kinds):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL RUN INFO TESTS PASSED")


if __name__ == "__main__":
    main()
