#!/usr/bin/env python3
"""F0.9 guard: validate a Salmon quant dir against the library-type contract.

Pure Python (stdlib). The logic lives here so it can be unit-tested with synthetic
JSON fixtures; the Nextflow module is a thin wrapper.

Criteria
--------
Hard failures:
  H1  `lib_format_counts.expected_format == cmd_info.libType`
      (detects Salmon ignoring/overriding the `-l` it was given — the `SR` bug).
  H2  `meta_info.library_types == [cmd_info.libType]`.
  H3  end-type coherence: a paired run (two read files) needs a paired libtype, a
      single-end run a single-end libtype (catches the archived `SR` on PE data).
  H4  (optional) `cmd_info.libType == --declared` (the contract libtype).

Reported as QC only (never pass/fail):
  - `consistent_fraction` = num_frags_with_concordant_consistent_mappings /
    num_assigned_fragments (0.945–0.979 in the ISR run);
  - `compatible_fragment_ratio` (= 1.000 under `-l ISR`, so it cannot flag
    orientation);
  - the dominant observed format.

Orientation is NOT decided here: that is Guard 1 (`CHECK_STRANDEDNESS`, RSeQC).

CLI:
  check_salmon_libtype.py --quant <dir> [--declared ISR]
"""

import argparse
import json
import os
import sys

PAIRED_ORIENT = set("IMO")
SINGLE = {"U", "SF", "SR"}
LIB_KEYS = ["MSF", "OSF", "ISF", "MSR", "OSR", "ISR", "SF", "SR", "MU", "OU", "IU", "U"]


def is_paired_libtype(lt):
    """Paired libtypes are <I|M|O><SF|SR|U> (e.g. ISR, ISF, IU, MU); single-end
    are U / SF / SR."""
    return (bool(lt) and lt not in SINGLE and len(lt) >= 2
            and lt[0] in PAIRED_ORIENT and lt[1:] in {"SF", "SR", "U"})


def count_read_files(lfc):
    rf = lfc.get("read_files")
    if not rf:
        return None
    inner = str(rf).strip()
    if inner.startswith("["):
        inner = inner[1:-1]
        return len([p for p in inner.split(",") if p.strip()])
    return 1


def consistent_fraction(lfc):
    a = lfc.get("num_assigned_fragments")
    c = lfc.get("num_frags_with_concordant_consistent_mappings")
    if a and c is not None and a > 0:
        return c / a
    return None


def dominant_format(lfc):
    counts = {k: lfc[k] for k in LIB_KEYS if k in lfc}
    if not counts:
        return None, {}
    return max(counts, key=counts.get), counts


def check(qdir, declared=None):
    """-> (errors: list, qc: dict). No sys.exit; callers decide."""
    errors = []
    cmd_path = os.path.join(qdir, "cmd_info.json")
    lfc_path = os.path.join(qdir, "lib_format_counts.json")
    if not os.path.exists(cmd_path):
        return ["missing %s" % cmd_path], {}
    if not os.path.exists(lfc_path):
        return ["missing %s" % lfc_path], {}
    with open(cmd_path) as fh:
        cmd = json.load(fh)
    with open(lfc_path) as fh:
        lfc = json.load(fh)
    meta = None
    for p in ("aux_info/meta_info.json", "meta_info.json"):
        fp = os.path.join(qdir, p)
        if os.path.exists(fp):
            with open(fp) as fh:
                meta = json.load(fh)
            break

    libtype = cmd.get("libType")
    expected = lfc.get("expected_format")

    if expected != libtype:
        errors.append("H1: declared libType '%s' but expected_format '%s'" % (libtype, expected))
    if meta is not None:
        lts = meta.get("library_types") or []
        if lts != [libtype]:
            errors.append("H2: meta library_types %s != declared '%s'" % (lts, libtype))
    nrf = count_read_files(lfc)
    if nrf is not None:
        if nrf == 2 and not is_paired_libtype(libtype or ""):
            errors.append("H3: %d read files (paired) but libtype '%s' is single-end" % (nrf, libtype))
        if nrf == 1 and is_paired_libtype(libtype or ""):
            errors.append("H3: 1 read file (single) but libtype '%s' is paired-end" % libtype)
    if declared is not None and libtype != declared:
        errors.append("H4: cmd libType '%s' != contract '%s'" % (libtype, declared))

    dom, counts = dominant_format(lfc)
    qc = {
        "cmd_libtype": libtype,
        "expected_format": expected,
        "consistent_fraction": consistent_fraction(lfc),
        "compatible_fragment_ratio": lfc.get("compatible_fragment_ratio"),
        "dominant_format": dom,
        "counts": counts,
    }
    return errors, qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant", required=True, help="Salmon quant directory")
    ap.add_argument("--declared", default=None, help="contract libtype (e.g. ISR)")
    args = ap.parse_args()
    sample = os.path.basename(args.quant.rstrip("/"))
    errors, qc = check(args.quant, args.declared)
    cf = qc.get("consistent_fraction")
    print("CHECK_SALMON_LIBTYPE [%s]: cmd=%s expected=%s dominant=%s "
          "consistent_fraction=%s compatible_fragment_ratio=%s"
          % (sample, qc.get("cmd_libtype"), qc.get("expected_format"),
             qc.get("dominant_format"),
             "NA" if cf is None else "%.4f" % cf,
             qc.get("compatible_fragment_ratio")))
    for e in errors:
        print("ERROR CHECK_SALMON_LIBTYPE [%s]: %s" % (sample, e), file=sys.stderr)
    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
