// REFCOMPAT: Guard that external conservation/synteny resources match
// the reference genome (params.genome) before any bigWig/chain is read.
//
// Runs very early in Phase A, only when at least one of
// conservation_bigwig / synteny_chain_file is configured, and fails fast
// (reading headers only) rather than after hours of compute. Nothing here
// is species-specific: comparisons normalize the chr prefix and compare
// chromosome lengths, which is what actually catches an assembly mismatch.
//
// Outputs a published log that also serves as methods/supplementary evidence.

process CHECK_REFERENCE_COMPATIBILITY {

    label 'process_single'
    memory 2.GB

    input:
    path genome_fa              // reference fasta; .fai generated in the task dir
    path conservation_bw, optional: true   // staged conservation bigWig; absent when unset
    path synteny_chain, optional: true     // staged liftOver chain; absent when unset

    output:
    path "reference_compatibility.log", emit: compatibility_log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    echo "reference_compatibility check: \$(date -Is)" > reference_compatibility.log
    echo "species        = ${params.species}" | tee -a reference_compatibility.log
    echo "genome_build   = ${params.genome_build ?: 'UNSET'}" | tee -a reference_compatibility.log

    # Normalized copies (Groovy renders an unset param as the literal "null",
    # and an optional no-file path input can render as empty, "null" or "[]").
    # Reduce all of them to a single empty string here, then decide the SKIP.
    bigwig="${params.conservation_bigwig ?: ''}"
    chain="${params.synteny_chain_file ?: ''}"
    # Paths staged by Nextflow for container visibility override the params
    # string whenever the file was actually passed in.
    if [ -n "${conservation_bw}" ] && [ "${conservation_bw}" != "null" ] && [ "${conservation_bw}" != "[]" ]; then bigwig="${conservation_bw}"; fi
    if [ -n "${synteny_chain}" ] && [ "${synteny_chain}" != "null" ] && [ "${synteny_chain}" != "[]" ]; then chain="${synteny_chain}"; fi

    # Chromosome size index of the reference genome (created locally, never
    # mutates the reference directory even if it is writable).
    samtools faidx "${genome_fa}"
    awk -v OFS='\\t' '\$1 !~ /^#/ {print \$1, \$2}' "${genome_fa}".fai > genome_sizes.tsv

    if [ -z "\$bigwig" ] && [ -z "\$chain" ]; then
        # Nothing external configured: nothing to validate.
        echo "SKIP: neither conservation_bigwig nor synteny_chain_file set" | tee -a reference_compatibility.log
        exit 0
    fi

    if [ -z "${params.genome_build ?: ''}" ]; then
        echo "ERROR: genome_build is required to validate external references." | tee -a reference_compatibility.log
        echo "Set params.genome_build (e.g. 'GRCh38.p14') in your params.yaml." >> reference_compatibility.log
        exit 1
    fi

    # Staged paths (from Nextflow inputs) override the params string; export so
    # the python heredoc (single-quoted, no bash expansion) can read them.
    export COMPAT_BIGWIG="$bigwig"
    export COMPAT_CHAIN="$chain"

    python3 - <<'PYEOF'
import os

def norm(name):
    # Cosmetic only: strip the UCSC 'chr' prefix for cross-checks.
    return name[3:] if name.startswith('chr') else name

# Staged paths (or empty string) exported by the surrounding bash script; single
# quotes here would block bash expansion, so we read them from the environment.
# An optional no-file input may surface as the literal "null" — coerce to empty.
bigwig = os.environ.get('COMPAT_BIGWIG', '').strip()
chain  = os.environ.get('COMPAT_CHAIN', '').strip()
bigwig = '' if bigwig in ('none', '[]', 'null') else bigwig
chain  = '' if chain in ('none', '[]', 'null') else chain

fai = {}
with open('genome_sizes.tsv') as fh:
    for line in fh:
        c, s = line.rstrip('\n').split('\t')
        fai[c] = int(s)

def warn(m):
    with open('reference_compatibility.log', 'a') as out:
        out.write(m + '\n')

if bigwig:
    import pyBigWig
    bw = pyBigWig.open(bigwig)
    chroms = bw.chroms()  # name -> size
    bw.close()
    warn('  conservation_bigwig chroms (n=%d)' % len(chroms))

    # Map by normalized (chr-stripped) name; keep original for reporting.
    by_norm = {norm(k): (k, v) for k, v in chroms.items()}
    mism = False
    for c, gsize in fai.items():
        if norm(c) not in by_norm:
            continue
        bname, bsize = by_norm[norm(c)]
        if bsize != gsize:
            mism = True
            warn('  MISMATCH chr %s: genome=%d bigwig=%d (%s)' % (c, gsize, bsize, bname))
    if mism:
        warn('FATAL: chromosome lengths differ between genome.fa and the conservation bigWig — they are not the same assembly.')
        raise SystemExit(1)
    warn('  conservation bigWig chromosome lengths match genome.fa after chr-prefix normalization.')

if chain:
    import gzip
    chain_f = chain
    opener = gzip.open if chain_f.endswith(('.gz', '.bgz')) else open
    first = None
    with opener(chain_f, 'rt') as fh2:
        for line in fh2:
            if line.startswith('chain '):
                first = line.split()
                break
    if first is None:
        warn('  WARNING: no chain blocks found in synteny_chain_file; cannot validate')
    else:
        # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd
        sides = {
            't': (first[1], int(first[2])),
            'q': (first[5], int(first[7])),
        }
        fai_norm = {norm(c): gsize for c, gsize in fai.items()}
        matched = None
        for side in ('t', 'q'):
            sname, sz = sides[side]
            if norm(sname) in fai_norm and fai_norm[norm(sname)] == sz:
                matched = side
                break
        if matched is None:
            warn('  FATAL: synteny chain lengths do not match genome.fa on either side — chain is for a different build.')
            raise SystemExit(1)
        target = 'q' if matched == 't' else 't'
        tn2, tsz2 = sides[target]
        warn('  chain matches reference genome on %s side (%s=%d)' % (matched, sides[matched][0], sides[matched][1]))
        warn('  chain declares target side: %s (size %d)' % (tn2, tsz2))

print('reference_compatibility OK')
PYEOF

    # Chain target log line — only write it when a chain was actually checked,
    # so the published log can be trusted as evidence for methods/supplementary.
    if [ -n "$chain" ]; then
        echo "  chain check complete (see python log above)" >> reference_compatibility.log
        echo "  synteny target species (params): ${params.synteny_target_species ?: 'unset'}" >> reference_compatibility.log
    fi
    """
}