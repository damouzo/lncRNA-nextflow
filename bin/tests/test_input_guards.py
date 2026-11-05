#!/usr/bin/env python3
"""0.3 static guard: validate_input.nf must reject inconsistent input.

Nextflow is not available in the offline suite, so this locks the guard
expressions in place: a future edit that drops the lane consistency or the
numerator/denominator membership check fails here instead of silently letting a
typo reach quantification.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = open(os.path.join(ROOT, "workflows", "validate_input.nf")).read()

fails = 0


def check(cond, msg):
    global fails
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        fails += 1


check("condition differs across lanes" in SRC,
      "lanes must agree on condition")
check("batch differs across lanes" in SRC,
      "lanes must agree on batch")
check("numerator" in SRC and "is not a condition in the samplesheet" in SRC,
      "contrast numerator must exist in the samplesheet")
check("denominator" in SRC and SRC.count("is not a condition in the samplesheet") == 2,
      "contrast denominator must exist in the samplesheet")

# Freeze gate: novel_classes is its own param and Tier 1 must be a subset of it.
check("params.novel_classes" in SRC and "are not in novel_classes" in SRC,
      "tier1_classes must be a subset of novel_classes")

FREEZE = open(os.path.join(ROOT, "modules", "local", "annotation_freeze.nf")).read()
check("params.novel_classes" in FREEZE and "params.tier1_classes" not in FREEZE,
      "ANNOTATION_FREEZE gates on novel_classes, not tier1_classes")
check("class_code %in% NOVEL_CLASSES" in FREEZE,
      "ANNOTATION_FREEZE filters the frozen set by class_code")

CONFIG = open(os.path.join(ROOT, "nextflow.config")).read()


def config_list(name):
    import re
    m = re.search(r"^\s*%s\s*=\s*\[([^\]]*)\]" % name, CONFIG, re.M)
    return {x.strip().strip("'\"") for x in m.group(1).split(",") if x.strip()} if m else None


novel, tier1 = config_list("novel_classes"), config_list("tier1_classes")
check(novel is not None and tier1 is not None and tier1 <= novel,
      "config defaults: tier1_classes subset of novel_classes")
check(novel is not None and not (novel & {"c", "k", "m", "n", "o"}),
      "default novel_classes excludes reference-derived (c/k/m/n) and o")

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL INPUT GUARD TESTS PASSED")
