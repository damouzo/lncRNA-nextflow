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

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL INPUT GUARD TESTS PASSED")
