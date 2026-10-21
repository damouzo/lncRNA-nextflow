#!/usr/bin/env bash
# Emit a versions.yml block for run_info: process, tool, version, container.
# Never fails: a missing tool or a tool without --version records "not_available",
# so this can be appended to any process without risking the task.
#
# usage: emit_versions.sh <process> <container> <tool> [<version-cmd>...]
set -u

proc="${1:-unknown}"
cont="${2:-}"
if [ "$#" -ge 2 ]; then
    shift 2
elif [ "$#" -eq 1 ]; then
    shift
fi

tool="${1:-unknown}"
if [ "$#" -ge 1 ]; then
    shift
fi

ver="not_available"
if [ "$#" -gt 0 ] && command -v "$1" >/dev/null 2>&1; then
    out="$("$@" 2>&1 | head -1 || true)"
    if [ -n "$out" ]; then
        ver="$out"
    fi
fi

# Keep the YAML single-line and quote-safe.
ver="$(printf '%s' "$ver" | tr '\t"' '  ')"
cont="$(printf '%s' "$cont" | tr '"' "'")"

{
    printf '"%s":\n' "$proc"
    printf '    tool: "%s"\n' "$tool"
    printf '    version: "%s"\n' "$ver"
    printf '    container: "%s"\n' "${cont:-none}"
} >> versions.yml

printf 'versions.yml: %s -> %s (%s)\n' "$proc" "$tool" "$ver"
