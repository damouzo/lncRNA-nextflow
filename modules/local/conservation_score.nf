// Conservation scoring per frozen lncRNA transcript.
// Reporting-only: adds transcript_id, mean_score, max_score, pct_bases_conserved.
// Skips entirely (info log, no error) when conservation_bigwig is unset —
// same contract as CATALOG_OVERLAP. Uses pyBigWig over exonic intervals
// from the frozen GTF (BED12 exon extraction, not the span locus).

process CONSERVATION_SCORE {

    label 'process_medium'
    memory 16.GB

    input:
    path frozen_gtf         // frozen novel lncRNA GTF (transcript/exon features)
    path bigwig, optional: true     // staged conservation bigWig; absent when unset

    output:
    path "conservation_scores.tsv", emit: scores

    script:
    """
    #!/usr/bin/env python3
    import sys

    bw_path = '${bigwig}'.strip()
    frozen  = '${frozen_gtf}'

    # An optional no-file input can render as empty, "null" or "[]" — treat
    # all as "not configured".
    if not bw_path or bw_path in ('null', 'none', '[]'):
        print('INFO: conservation_bigwig not configured; skipping CONSERVATION_SCORE')
        # Still emit an empty table so BUILD_GENE_CATALOG's join can rely on
        # its presence and produce all-NA columns instead of failing.
        with open('conservation_scores.tsv', 'w') as fh:
            fh.write('transcript_id\\tmean_score\\tmax_score\\tpct_bases_conserved\\n')
        sys.exit(0)

    import pyBigWig

    # Parse only transcript/exon features from the frozen GTF into BED12-style
    # intervals (chrom, start, end, strand, exons list). One record per transcript.
    trans   = {}   # tx -> {'chrom','start','end','strand','exons':[(s,e),...]}

    def extract_attr(attrs, key):
        import re
        m = re.search(r'(?:^|; )' + key + r' "([^"]+)"', attrs)
        return m.group(1) if m else None

    with open(frozen_gtf) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            p = line.rstrip('\\n').split('\\t')
            if len(p) < 9:
                continue
            ftype = p[2]
            attrs = p[8]
            chrom = p[0]
            start = int(p[3])
            end   = int(p[4])
            strand = p[6]
            if ftype == 'transcript':
                tid = extract_attr(attrs, 'transcript_id')
                if tid and tid not in trans:
                    trans[tid] = {'chrom': chrom, 'exons': [], 'strand': strand, 'start': start, 'end': end}
            elif ftype == 'exon':
                tid = extract_attr(attrs, 'transcript_id')
                if tid in trans:
                    trans[tid]['exons'].append((start, end))

    n_tx_with_exons = sum(1 for t in trans.values() if t['exons'])
    print(f'INFO: {len(trans)} transcripts, {n_tx_with_exons} with exons', file=sys.stderr)

    bw = pyBigWig.open(bw_path)
    chrom_sizes = bw.chroms()   # cached: cheap now, avoids recomputing per exon
    out = open('conservation_scores.tsv', 'w')
    out.write('transcript_id\\tmean_score\\tmax_score\\tpct_bases_conserved\\n')

    cutoff = float('${params.conservation_min_score}')
    for tid, t in sorted(trans.items()):
        if not t['exons']:
            out.write(tid + '\\tNA\\tNA\\tNA\\n')
            continue
        # Exon intervals (0-based half-open for pyBigWig)
        vals = []
        for (s, e) in t['exons']:
            if t['chrom'] not in chrom_sizes:
                vals.append([])
                continue
            vals.append(bw.values(t['chrom'], s - 1, e))
        flat = [v for sub in vals for v in sub if v is not None]
        if not flat:
            out.write(tid + '\\tNA\\tNA\\tNA\\n')
            continue
        mean = sum(flat) / len(flat)
        mx   = max(flat)
        pct  = sum(1 for v in flat if v >= cutoff) / len(flat)
        out.write(f'{tid}\\t{mean:.6f}\\t{mx:.6f}\\t{pct:.4f}\\n')

    out.close()
    bw.close()
    """
}