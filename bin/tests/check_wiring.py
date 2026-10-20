#!/usr/bin/env python3
"""Wiring guard (static, no Nextflow): conf/modules.config vs the publication contract.

Checks, all offline:
  1. Every publishDir in conf/modules.config matches docs/publication_map.json
     (same destination set per process, and no unmapped process either way).
  2. Every destination root is one of the map `roots` (raw | run_info | deliverables).
  3. No destination escapes the tree: no f1/, discovery/, quantification/, reports/,
     pipeline_info/ roots.
  4. Every deliverables/** destination is covered by a path in
     docs/deliverables_contract.json.
  5. `audit_only_full` destinations are gated to publish_level == 'full'.
  6. Script flags: every --flag a module passes to bin/<script>.py exists in that
     script's argparse (catches typos that only show up at run time).

Run: python3 bin/tests/check_wiring.py
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MODULES_CONFIG = os.path.join(ROOT, "conf", "modules.config")
PUBLICATION_MAP = os.path.join(ROOT, "docs", "publication_map.json")
CONTRACT = os.path.join(ROOT, "docs", "deliverables_contract.json")

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def read(path):
    with open(path) as fh:
        return fh.read()


def parse_modules_config(text):
    """-> {process_name: [{'dest': rel, 'enabled': expr_text}, ...]}"""
    heads = list(re.finditer(r"withName:\s*([A-Za-z0-9_]+)\s*\{", text))
    out = {}
    for i, m in enumerate(heads):
        name = m.group(1)
        start = m.end()
        end = heads[i + 1].start() if i + 1 < len(heads) else len(text)
        block = text[start:end]
        entries = []
        # split the block into publishDir entries (each list of maps)
        for pm in re.finditer(r"path:\s*\{\s*\"\$\{params\.outdir\}/([^\"}]*)\"\s*\}", block):
            dest = pm.group(1)
            # find an `enabled:` in the following ~250 chars
            window = block[pm.end():pm.end() + 250]
            em = re.search(r"enabled:\s*\{\s*([^}]*)\}", window)
            enabled = em.group(1) if em else "true"
            entries.append({"dest": dest, "enabled": enabled})
        if entries:
            out.setdefault(name, []).extend(entries)
    return out


def dest_root(dest):
    return dest.split("/", 1)[0]


def enabled_level(expr):
    if "==" in expr and "'full'" in expr:
        return "full"
    if "core_raw" in expr:
        return "core_raw"
    return "always"


def norm_contract_paths(contract):
    return [re.sub(r"<[^>]+>", "*", f["path"]) for f in contract.get("files", [])]


def check_destinations(pub, cfg):
    roots = set(pub.get("roots", []))
    forbidden = ("f1/", "discovery/", "quantification/", "reports/", "pipeline_info/")
    for name, entries in sorted(cfg.items()):
        for e in entries:
            dest = e["dest"]
            check(dest_root(dest) in roots,
                  "%s: dest '%s' root in %s" % (name, dest, sorted(roots)))
            check(not dest.startswith(forbidden),
                  "%s: dest '%s' not in a forbidden root" % (name, dest))


def check_matches_map(pub, cfg):
    mapped = {p["module"]: p for p in pub.get("processes", [])}
    check(set(cfg) <= set(mapped),
          "every publishDir process is in the map (unknown: %s)"
          % sorted(set(cfg) - set(mapped)))
    check(set(mapped) <= set(cfg),
          "every mapped process has a publishDir (missing: %s)"
          % sorted(set(mapped) - set(cfg)))
    for name, entries in sorted(cfg.items()):
        if name not in mapped:
            continue
        got = sorted({e["dest"] for e in entries})
        want = sorted(mapped[name]["dest"])
        check(got == want, "%s: destinations match map (%s)" % (name, got))


def check_contract_coverage(cfg, contract):
    # Contract paths are relative to <outdir>/deliverables.
    paths = norm_contract_paths(contract)
    for name, entries in sorted(cfg.items()):
        for e in entries:
            dest = e["dest"]
            if dest != "deliverables" and not dest.startswith("deliverables/"):
                continue
            rel = dest[len("deliverables"):].lstrip("/")
            covered = (rel == "" or
                       any(p == rel or p.startswith(rel.rstrip("/") + "/") for p in paths))
            check(covered, "%s: deliverables dest '%s' covered by contract"
                  % (name, dest))


def check_audit_levels(pub, cfg):
    audit = set(pub.get("audit_only_full", []))
    for name, entries in sorted(cfg.items()):
        for e in entries:
            lvl = enabled_level(e["enabled"])
            if e["dest"] in audit:
                check(lvl == "full",
                      "%s: audit dest '%s' gated to full (is %s)" % (name, e["dest"], lvl))
            else:
                check(lvl in ("core_raw", "always"),
                      "%s: non-audit dest '%s' not full-only" % (name, e["dest"]))


def _py_flags(py_path):
    if not os.path.isfile(py_path):
        return None
    flags = set()
    for m in re.finditer(r"add_argument\(\s*[\"'](--[A-Za-z0-9][A-Za-z0-9_-]*)", read(py_path)):
        flags.add(m.group(1))
    return flags


def check_script_flags():
    nf_files = [os.path.join(ROOT, "main.nf")]
    for d in ("modules/local", "workflows"):
        full = os.path.join(ROOT, d)
        if os.path.isdir(full):
            nf_files += [os.path.join(full, f) for f in sorted(os.listdir(full))
                         if f.endswith(".nf")]
    for path in nf_files:
        if not os.path.isfile(path):
            continue
        text = read(path)
        current = None
        by_script = {}
        for m in re.finditer(r"([A-Za-z0-9_]+\.py)|(--[A-Za-z0-9][A-Za-z0-9_-]*)", text):
            if m.group(1):
                current = m.group(1)
            elif m.group(2) and current:
                by_script.setdefault(current, set()).add(m.group(2))
        for script, flags in sorted(by_script.items()):
            known = _py_flags(os.path.join(ROOT, "bin", script))
            if known is None:
                continue
            unknown = sorted(f for f in flags if f not in known)
            check(not unknown,
                  "%s: flags known to bin/%s (unknown: %s)" % (path.replace(ROOT + "/", ""),
                                                               script, unknown))


def main():
    cfg = parse_modules_config(read(MODULES_CONFIG))
    pub = json.loads(read(PUBLICATION_MAP))
    contract = json.loads(read(CONTRACT))
    print("wiring: %d publishDir processes" % len(cfg))
    check_destinations(pub, cfg)
    check_matches_map(pub, cfg)
    check_contract_coverage(cfg, contract)
    check_audit_levels(pub, cfg)
    check_script_flags()
    print()
    if FAILS:
        print("WIRING FAILED: %d" % len(FAILS))
        return 1
    print("ALL WIRING CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
