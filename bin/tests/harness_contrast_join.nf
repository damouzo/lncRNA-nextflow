#!/usr/bin/env nextflow
// Regression harness for the contrasts channel bug (no containers, no data).
//
// One base item per group combined with N contrast items for the same key must
// yield N pairs. The previous `.join()` emitted a single pair per key and dropped
// the rest. Run with:  nextflow run bin/tests/harness_contrast_join.nf
nextflow.enable.dsl = 2

workflow {
    base = Channel.of(tuple('g', 'BASE'))
    cts  = Channel.of(tuple('g', 'c1'), tuple('g', 'c2'),
                      tuple('g', 'c3'), tuple('g', 'c4'))

    pairs = base.combine(cts, by: 0)
        .map { g, b, ct -> tuple(g, ct, b) }

    pairs
        .toList()
        .map { list ->
            if (list.size() != 4) {
                error("harness_contrast_join: expected 4 pairs, got ${list.size()}")
            }
            "harness_contrast_join: OK (${list.size()} pairs)"
        }
        .view()
}
