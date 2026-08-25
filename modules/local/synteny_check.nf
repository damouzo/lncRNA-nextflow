// Synteny check per frozen lncRNA transcript.
// Reporting-only: transcript_id, syntenic_locus (TRUE/FALSE/NA),
// syntenic_target_gene_id. NA = did not lift (informative, not a failure).
// Skips entirely (info log, empty table) when synteny_chain_file or
// synteny_target_gtf is unset — same contract as CATALOG_OVERLAP.
//
// Source BED is the exonic span (BED12) of each frozen transcript in genome
// coordinates. liftOver lifts to the target build; bedtools intersect marks
// whether the lifted locus overlaps a target-species lncRNA gene (biotypes
// reused from params.annotated_lncrna_biotypes).

process SYNTENY_CHECK {

    label 'process_medium'
    memory 12.GB

    input:
    path frozen_gtf         // frozen novel lncRNA GTF (transcript/exon features)
    path chain_file, optional: true     // staged liftOver chain; absent when unset
    path target_gtf, optional: true     // target-species annotated GTF; absent when unset

    output:
    path "synteny_scores.tsv", emit: scores
    path "synteny_unmapped.bed", emit: unmapped

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Coerce Groovy's "null" string to a real empty for the skip test; the
    # staged chain/target inputs override the params string when files are passed.
    # Both resources are required — skip cleanly if either is missing.
    chain="${chain_file}"
    target="${target_gtf}"
    if [ -z "\$chain" ] || [ "\$chain" = "null" ] || [ "\$chain" = "[]" ] || \
       [ -z "\$target" ] || [ "\$target" = "null" ] || [ "\$target" = "[]" ] || \
       [ -z "${params.synteny_chain_file ?: ''}" ] || [ -z "${params.synteny_target_gtf ?: ''}" ]; then
        echo "INFO: synteny_chain_file or synteny_target_gtf not configured; skipping SYNTENY_CHECK"
        printf 'transcript_id\\tsyntenic_locus\\tsyntenic_target_gene_id\\n' > synteny_scores.tsv
        : > synteny_unmapped.bed
        exit 0
    fi

    # ---- 1. Source BED12 from the frozen GTF (exonic span, genome coords) ----
    python3 - <<'PYEOF'
import re

trans = {}

def extract_attr(attrs, key):
    m = re.search(r'(?:^|; )' + key + r' "([^"]+)"', attrs)
    return m.group(1) if m else None

with open('${frozen_gtf}') as fh:
    for line in fh:
        if not line or line.startswith('#'):
            continue
        p = line.rstrip('\\n').split('\\t')
        if len(p) < 9:
            continue
        ftype, attrs = p[2], p[8]
        chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
        if ftype == 'transcript':
            tid = extract_attr(attrs, 'transcript_id')
            if tid and tid not in trans:
                trans[tid] = {'chrom': chrom, 'exons': [], 'strand': strand}
        elif ftype == 'exon':
            tid = extract_attr(attrs, 'transcript_id')
            if tid in trans:
                trans[tid]['exons'].append((start, end))

# UCSC BED12: chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
with open('synteny_source.bed', 'w') as out:
    for tid, t in trans.items():
        if not t['exons']:
            continue
        exons = sorted(set(t['exons']))
        tstart = min(e[0] for e in exons)
        tend   = max(e[1] for e in exons)
        starts = [s - tstart for s, _ in exons]
        sizes  = [e - s for s, e in exons]
        block  = ','.join(map(str, sizes)) + ','
        bstarts = ','.join(map(str, starts)) + ','
        out.write('\\t'.join([
            t['chrom'], str(tstart - 1), str(tend), tid, '0', t['strand'],
            str(tstart - 1), str(tend), '0',
            str(len(exons)), block, bstarts
        ]) + '\\n')
PYEOF

    # ---- 2. liftOver ----
    liftOver synteny_source.bed "\$chain" synteny_lifted.bed synteny_unmapped.bed

    # ---- 3. Target-species lncRNA features as BED6 (from target GTF) ----
    python3 - <<'PYEOF'
import re
BIOTYPES = {${params.annotated_lncrna_biotypes.collect { "\"$it\"" }.join(", ")}}
def extract_attr(attrs, key):
    m = re.search(r'(?:^|; )' + key + r' "([^"]+)"', attrs)
    return m.group(1) if m else None

with open('${target_gtf}') as fh, open('target_lnc_exons.bed', 'w') as out:
    for line in fh:
        if not line or line.startswith('#'):
            continue
        p = line.rstrip('\\n').split('\\t')
        if len(p) < 9:
            continue
        if p[2] != 'exon':
            continue
        attrs = p[8]
        biotype = extract_attr(attrs, 'gene_biotype')
        if biotype not in BIOTYPES:
            continue
        gid = extract_attr(attrs, 'gene_id')
        if not gid:
            continue
        out.write('\\t'.join([p[0], str(int(p[3]) - 1), p[4], gid, '0', p[6]]) + '\\n')
PYEOF

    # ---- 4. Intersect lifted loci with target lncRNA exons; build output ----
    # -wao keeps single-line-per-transcript (BED12 in A) + the matched target
    # gene (cols 13-18) even with zero overlap, so we can fill TRUE/FALSE/NA.
    bedtools intersect -wao -a synteny_lifted.bed -b target_lnc_exons.bed > synteny_intersect.bed

    python3 - <<'PYEOF'
lifted = {}
with open('synteny_lifted.bed') as fh:
    for line in fh:
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) >= 4:
            lifted[parts[3]] = True

# transcript id -> matched target gene id from the -wao file (B fields cols 13-17)
gene_of = {}
with open('synteny_intersect.bed') as fh:
    for line in fh:
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) < 13:
            continue
        tid = parts[3]
        # B (target) fields start after A's 12 fields: index 12=chrom,13=start,
        # 14=end, 15=name (target gene id), 16=score, 17=strand
        bgene = parts[15] if len(parts) > 15 else None
        if bgene and bgene != '.' and tid not in gene_of:
            gene_of[tid] = bgene

with open('synteny_scores.tsv', 'w') as out:
    out.write('transcript_id\\tsyntenic_locus\\tsyntenic_target_gene_id\\n')
    with open('synteny_source.bed') as src:
        for line in src:
            tid = line.rstrip('\\n').split('\\t')[3]
            if tid not in lifted:
                out.write(f'{tid}\\tNA\\tNA\\n')
            elif tid in gene_of:
                out.write(f'{tid}\\tTRUE\\t{gene_of[tid]}\\n')
            else:
                out.write(f'{tid}\\tFALSE\\tNA\\n')
PYEOF
    """
}