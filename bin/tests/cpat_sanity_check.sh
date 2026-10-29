#!/usr/bin/env bash
# 3.2 sanity check: does the CPAT stock human model separate Ensembl 115
# protein_coding transcripts from lncRNA transcripts at the configured threshold?
# Reports sensitivity / specificity / accuracy so the 0.364 cutoff can be
# validated against THIS annotation (never a threshold from memory).
#
# Run inside the coding_potential image (needs gffread, cpat.py, python3):
#   apptainer exec $CODING_IMG bash bin/tests/cpat_sanity_check.sh \
#       <annotation.gtf> <genome.fa> [n_per_class] [threshold]
#
# Env overrides: CPAT_HEXAMER, CPAT_LOGIT.
set -eu

GTF=${1:?usage: cpat_sanity_check.sh <gtf> <genome.fa> [n] [threshold]}
FA=${2:?usage: cpat_sanity_check.sh <gtf> <genome.fa> [n] [threshold]}
N=${3:-300}
THR=${4:-0.364}

REPO=$(cd "$(dirname "$0")/../.." && pwd)
HEX=${CPAT_HEXAMER:-$REPO/assets/cpat/Human_Hexamer.tsv}
LOGIT=${CPAT_LOGIT:-$REPO/assets/cpat/Human_logitModel.RData}
for f in "$GTF" "$FA" "$HEX" "$LOGIT"; do
    [ -s "$f" ] || { echo "missing input: $f" >&2; exit 2; }
done

work=$(mktemp -d)
cd "$work"
echo "workdir: $work"
echo "threshold: $THR  n_per_class: $N"

# ---- pick transcript ids of each class ----
awk -F'\t' '$3=="transcript" && $9 ~ /transcript_biotype "protein_coding"/' "$GTF" \
    | head -n "$N" | grep -o 'transcript_id "[^"]*"' | sed 's/.*"\(.*\)"/\1/' \
    | sort -u > coding.ids
awk -F'\t' '$3=="transcript" && $9 ~ /transcript_biotype "lncRNA"/' "$GTF" \
    | head -n "$N" | grep -o 'transcript_id "[^"]*"' | sed 's/.*"\(.*\)"/\1/' \
    | sort -u > lncrna.ids
echo "ids: coding=$(wc -l < coding.ids) lncrna=$(wc -l < lncrna.ids)"

# ---- pull the full feature set for those transcripts ----
extract() {  # ids gtf out.gtf
    awk -F'\t' -v idf="$1" '
        NR==FNR { id[$1]; next }
        {
            if (match($9, /transcript_id "[^"]+"/)) {
                tid = substr($9, RSTART+15, RLENGTH-16)
                if (tid in id) print
            }
        }' "$1" "$2" > "$3"
}
extract coding.ids "$GTF" coding.gtf
extract lncrna.ids "$GTF" lncrna.gtf

gffread -w coding.fa -g "$FA" coding.gtf
gffread -w lncrna.fa -g "$FA" lncrna.gtf
n_coding=$(grep -c '^>' coding.fa || true)
n_lncrna=$(grep -c '^>' lncrna.fa || true)
echo "sequences: coding=$n_coding lncrna=$n_lncrna"

cpat.py -g coding.fa -o cpat_coding -d "$LOGIT" -x "$HEX" --antisense
cpat.py -g lncrna.fa -o cpat_lncrna -d "$LOGIT" -x "$HEX" --antisense

python3 - "$THR" "$n_coding" "$n_lncrna" <<'PY'
import sys, csv
thr = float(sys.argv[1]); n_coding = int(sys.argv[2]); n_lncrna = int(sys.argv[3])

def called_coding(best):
    ids = set()
    try:
        with open(best) as fh:
            r = csv.reader(fh, delimiter="\t")
            next(r, None)  # header
            for row in r:
                if len(row) >= 11 and row[10] not in ("", "NA"):
                    if float(row[10]) >= thr:
                        ids.add(row[0].upper())
    except FileNotFoundError:
        pass
    return ids

tp_ids = called_coding("cpat_coding.ORF_prob.best.tsv")
fp_ids = called_coding("cpat_lncrna.ORF_prob.best.tsv")
tp = len(tp_ids); fn = n_coding - tp
fp = len(fp_ids); tn = n_lncrna - fp
sens = tp / (tp + fn) if (tp + fn) else float("nan")
spec = tn / (tn + fp) if (tn + fp) else float("nan")
acc = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) else float("nan")
print("RESULT\tthreshold\tn_coding\tn_lncrna\tTP\tFN\tFP\tTN\tsensitivity\tspecificity\taccuracy")
print("RESULT\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f"
      % (thr, n_coding, n_lncrna, tp, fn, fp, tn, sens, spec, acc))
PY
