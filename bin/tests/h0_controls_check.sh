#!/usr/bin/env bash
# 4) H0 real-data check: run the corrected deseq2_contrasts.R on the OLD
# NB4/MSCline tximport counts, with the name-aligned coldata and --allow-no-offset,
# and report DDX41 / DHX34 across the 8 sh contrasts.
# Writes into a temp dir OUTSIDE results/.
#
#   apptainer exec $R_BIOC_IMG bash bin/tests/h0_controls_check.sh <run_dir> <outdir>
set -eu

RUN=${1:?usage: h0_controls_check.sh <lncrna_run_dir> <outdir>}
OUT=${2:?usage: h0_controls_check.sh <lncrna_run_dir> <outdir>}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
SCRIPT=$REPO/bin/deseq2_contrasts.R
FILT=$RUN/results/raw/tximport_filter
COLD=$RUN/results/raw/per_group_coldata
mkdir -p "$OUT"

for g in NB4 MSCline; do
    grep "^${g}_" "$RUN/comparisons.csv" \
        | awk -F, '{print $1"\t"$2"\t"$3"\t"}' > "$OUT/contrasts_${g}.tsv"
    echo "== ${g}: $(wc -l < "$OUT/contrasts_${g}.tsv") contrasts =="
    Rscript "$SCRIPT" \
        --group "$g" \
        --counts "$FILT/${g}_tximport_counts.rds" \
        --coldata "$COLD/coldata_${g}.csv" \
        --contrasts "$OUT/contrasts_${g}.tsv" \
        --design-formula '~ batch + condition' \
        --fdr 0.05 --allow-no-offset true \
        --outdir "$OUT/run_${g}"
    grep -m1 length_offset "$OUT/run_${g}/${g}_deseq2_results/${g}_length_offset.txt" || true
done

python3 - "$OUT" <<'PY'
import csv, glob, os, sys
out = sys.argv[1]
genes = {"ENSG00000183258": "DDX41", "ENSG00000134815": "DHX34"}
order = []
rows = []
for g in ("NB4", "MSCline"):
    d = os.path.join(out, "run_%s" % g, "%s_deseq2_results" % g)
    for f in sorted(glob.glob(os.path.join(d, "*_DE.tsv"))):
        ct = os.path.basename(f)[:-len("_DE.tsv")]
        rec = {r["gene_id"]: r for r in csv.DictReader(open(f), delimiter="\t")}
        for gid, sym in genes.items():
            r = rec.get(gid, {})
            lfc = r.get("log2FoldChange", "NA")
            padj = r.get("padj", "NA")
            rows.append((g, ct, sym, lfc, padj))

hdr = ["group", "contrast", "gene", "log2FoldChange", "padj"]
with open(os.path.join(out, "h0_controls.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(hdr)
    w.writerows(rows)

print("\n%-9s %-34s %-6s %14s %12s" % tuple(hdr))
ok = 0
total = 0
for g, ct, sym, lfc, padj in rows:
    # only the applicable control: DDX41 in DDX41sh*, DHX34 in DHX34sh*
    if sym == "DDX41" and "DDX41sh" in ct or sym == "DHX34" and "DHX34sh" in ct:
        total += 1
        try:
            good = float(lfc) < 0 and float(padj) < 0.05
        except ValueError:
            good = False
        ok += 1 if good else 0
        flag = "ok" if good else "FAIL"
        print("%-9s %-34s %-6s %14s %12s   %s" % (g, ct, sym, lfc, padj, flag))
print("\nH0 controls correct: %d/%d" % (ok, total))
PY
