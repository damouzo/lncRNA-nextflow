// Single strandedness contract for the whole pipeline.
//
// `library_strandedness` (reverse | forward | unstranded) is the single source of
// truth; Salmon, featureCounts and StringTie derive their flags from it. The
// legacy explicit params (`salmon_libtype`, `featurecounts_strandness`) survive
// only as coherence overrides: if set they must match what the contract derives,
// otherwise the pipeline fails instead of silently running a mismatched library.
//
// reverse  -> ISR (PE) / SR (SE), featureCounts -s 2, StringTie --rf
// forward  -> ISF (PE) / SF (SE), featureCounts -s 1, StringTie --fr
// unstranded -> IU (PE) / U (SE), featureCounts -s 0, StringTie (no flag)
class Strandedness {

    static final List MODES = ['reverse', 'forward', 'unstranded']

    private static final Map SALMON_PE = [reverse: 'ISR', forward: 'ISF', unstranded: 'IU']
    private static final Map SALMON_SE = [reverse: 'SR',  forward: 'SF',  unstranded: 'U']
    private static final Map FC        = [reverse: 2,     forward: 1,     unstranded: 0]
    private static final Map STRINGTIE = [reverse: '--rf', forward: '--fr', unstranded: '']

    static String mode(Object value) {
        def m = (value == null ? 'reverse' : value.toString().trim().toLowerCase())
        if (m in ['reverse-stranded', 'reversely-stranded', 'fr-firststrand', 'fr-first-strand', 'rf', 'reverse_stranded']) m = 'reverse'
        if (m in ['forward-stranded', 'fr-secondstrand', 'fr-second-strand', 'fr', 'stranded', 'yes', 'forward_stranded']) m = 'forward'
        if (m in ['unstranded', 'none', 'no', 'false']) m = 'unstranded'
        if (!MODES.contains(m)) {
            throw new IllegalArgumentException(
                "library_strandedness='${value}' is not valid; use one of ${MODES}")
        }
        return m
    }

    private static String derivedSalmon(String m, boolean paired) {
        return paired ? SALMON_PE[m] : SALMON_SE[m]
    }

    // Returns the Salmon -l string, honouring an explicit override only when coherent.
    static String salmonLibtype(Object modeValue, Object override, boolean paired) {
        def m        = mode(modeValue)
        def derived  = derivedSalmon(m, paired)
        def o        = (override == null ? '' : override.toString().trim())
        if (o == '' || o == 'null') return derived
        if (o != derived) {
            throw new IllegalArgumentException(
                "salmon_libtype='${o}' conflicts with library_strandedness='${m}' " +
                "(expected '${derived}' for ${paired ? 'paired-end' : 'single-end'} data). " +
                "Remove the override or make it coherent.")
        }
        return derived
    }

    static int featurecountsStrand(Object modeValue, Object override) {
        def m       = mode(modeValue)
        def derived = FC[m] as int
        def o       = (override == null ? '' : override.toString().trim())
        if (o == '' || o == 'null') return derived
        int oi
        try {
            oi = Integer.parseInt(o)
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("featurecounts_strandness='${o}' is not an integer")
        }
        if (oi != derived) {
            throw new IllegalArgumentException(
                "featurecounts_strandness=${oi} conflicts with library_strandedness='${m}' " +
                "(expected ${derived}). Remove the override or make it coherent.")
        }
        return derived
    }

    // Full StringTie flag line ('' when unstranded). Callers place it on its own
    // line so the command stays byte-identical to the pre-contract pipeline.
    static String stringtieFlag(Object modeValue) {
        return STRINGTIE[mode(modeValue)]
    }
}
