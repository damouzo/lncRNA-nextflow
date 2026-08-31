#!/usr/bin/env python3
"""Shared chr-prefix normalization for conservation/synteny/catalog-overlap.

Single source of truth for chr-prefix detection and chromosome name
normalization.  All three modules (CONSERVATION_SCORE, SYNTENY_CHECK,
CATALOG_OVERLAP) use this script so the logic never diverges.

Usage:
  normalize_chrom.py detect <tsv-file>     print "chr" or "no_chr"
  normalize_chrom.py norm <tsv-file> <name>  print <name> normalised to file convention
  normalize_chrom.py chain-qname <chain>   print "chr" or "no_chr" for target side of chain
"""

import sys
import gzip


def _has_chr_prefix(name):
    return name.startswith('chr')


def _normalise(has_chr, name):
    name_has_chr = _has_chr_prefix(name)
    if has_chr and not name_has_chr:
        return 'chr' + name
    if not has_chr and name_has_chr:
        return name[3:]
    return name


def detect_has_chr(path):
    """Check whether a tab-separated file uses chr-prefixed chromosome names.

    For GTF files scans all transcript lines.  For BED and other non-GTF files
    samples the first 20 lines and uses majority vote so a single atypical
    scaffold on line 1 cannot flip the answer.
    """
    opener = gzip.open if path.endswith(('.gz', '.bgz')) else open
    chroms = []
    is_gtf = False
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 3:
                continue
            chrom = p[0]
            if not chrom:
                continue
            if p[2] == 'transcript':
                is_gtf = True
                chroms.append(chrom)
            elif not is_gtf and len(chroms) < 20:
                chroms.append(chrom)
    if not chroms:
        return False
    n_chr = sum(1 for c in chroms if _has_chr_prefix(c))
    return n_chr > len(chroms) / 2


def detect_chain_qname(chain_path):
    """Read the first chain header and report whether qName has chr prefix."""
    opener = gzip.open if chain_path.endswith(('.gz', '.bgz')) else open
    with opener(chain_path, 'rt') as fh:
        for line in fh:
            if line.startswith('chain '):
                qname = line.split()[7]
                return _has_chr_prefix(qname)
    return False


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(1)

    cmd = sys.argv[1]

    if cmd == 'detect':
        print('chr' if detect_has_chr(sys.argv[2]) else 'no_chr')

    elif cmd == 'norm':
        has_chr = detect_has_chr(sys.argv[2])
        print(_normalise(has_chr, sys.argv[3]))

    elif cmd == 'chain-qname':
        print('chr' if detect_chain_qname(sys.argv[2]) else 'no_chr')

    else:
        sys.exit(1)