#!/usr/bin/env bash
# Emit one versions.yml block per R package, in a single Rscript call.
# Never fails: a missing package records "not_available".
#
# usage: emit_r_versions.sh <process> <container> <pkg> [<pkg>...]
set -u

proc="${1:-unknown}"
cont="${2:-}"
shift 2 || true

if [ "$#" -eq 0 ] || ! command -v Rscript >/dev/null 2>&1; then
    exit 0
fi

# One line per package: "name<TAB>version"; packageVersion() does not need a load.
Rscript -e 'for (p in commandArgs(TRUE)) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not_available")
    cat(p, "\t", v, "\n", sep = "")
}' "$@" 2>/dev/null | while IFS=$'\t' read -r name ver; do
    [ -n "$name" ] || continue
    emit_versions.sh "$proc" "$cont" "$name" printf '%s' "$ver"
done
