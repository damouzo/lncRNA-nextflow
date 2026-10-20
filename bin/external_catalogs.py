#!/usr/bin/env python3
"""F2.2: external catalog registration (and optional download).

The download is deliberately NOT part of the pipeline: this helper is run by hand.
It records, per source, the pinned version, URL, date, md5 and feature count in
`run_info/external_catalogs.tsv` so the overlap (F2.1) is reproducible.

Sources are described by a JSON spec the investigator maintains (pinned versions
and URLs). Modes:
  --check     : md5 + feature count for every source whose file is present
  --download  : fetch `url` into the dest dir (only when `url` is set)
  --register  : write external_catalogs.tsv (always printed for --check)

Stdlib only. Feature counts use bin/catalog_overlap.parse_bed (BED6/BED12).

CLI:
  external_catalogs.py --spec docs/external_catalogs.json --dest assets/catalogs \
      [--mode check|download] --out run_info/external_catalogs.tsv
"""

import argparse
import hashlib
import json
import os
import sys
import urllib.request
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalog_overlap as O  # noqa: E402


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def n_features(path):
    return len(O.parse_bed(path))


def register(entries, dest):
    rows = []
    for e in entries:
        path = e.get("file") or os.path.join(dest, "%s.bed" % e["name"])
        exists = os.path.exists(path)
        rows.append({
            "name": e["name"],
            "version": e.get("version", ""),
            "url": e.get("url", ""),
            "date": e.get("date", str(date.today())),
            "md5": md5sum(path) if exists else "",
            "n_features": n_features(path) if exists else 0,
            "file": path,
            "present": "TRUE" if exists else "FALSE",
        })
    return rows


def download(entries, dest):
    os.makedirs(dest, exist_ok=True)
    for e in entries:
        if not e.get("url"):
            print("  %s: no url in spec — skipped" % e["name"])
            continue
        target = os.path.join(dest, "%s.bed" % e["name"])
        print("  %s: downloading %s -> %s" % (e["name"], e["url"], target))
        urllib.request.urlretrieve(e["url"], target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--dest", default="assets/catalogs")
    ap.add_argument("--mode", choices=["check", "download"], default="check")
    ap.add_argument("--out", default="external_catalogs.tsv")
    args = ap.parse_args()

    with open(args.spec) as fh:
        entries = json.load(fh)["sources"]

    if args.mode == "download":
        download(entries, args.dest)

    rows = register(entries, args.dest)
    cols = ["name", "version", "url", "date", "md5", "n_features", "file", "present"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_present = sum(1 for r in rows if r["present"] == "TRUE")
    abs_hits = sum(r["n_features"] for r in rows)
    print("external_catalogs: %d/%d configured, %d features total -> %s"
          % (n_present, len(rows), abs_hits, args.out))


if __name__ == "__main__":
    main()
