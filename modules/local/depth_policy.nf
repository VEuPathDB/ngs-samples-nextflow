import groovy.transform.Field

/*
 * Pure depth-policy functions. No params access, no file IO — everything is passed in
 * so this can be unit tested without a pipeline run.
 *
 * Units: every read count here is a FRAGMENT count, matching metrics.totalReads, which
 * MEASURE_SAMPLE derives from R1 records alone. A paired fragment contributes
 * mateCount * readLength bases, so coverage math must go through basesPerFragment
 * rather than readLength.
 */

// Cost/runtime sanity bounds: a floor so tiny genomes still get usable depth,
// a ceiling so huge genomes don't blow up downstream alignment cost.
@Field final long MIN_TARGET_FRAGMENTS = 1000000L
@Field final long MAX_TARGET_FRAGMENTS = 100000000L

// Assays whose depth is quoted as a flat fragment count rather than genome coverage,
// so these targets need no pairing adjustment. RNA-seq depth scales with transcriptome
// complexity; ChIP-seq depth accrues under peaks. Neither scales with genome size, which
// is why running them through the coverage path below would use the wrong axis.
// ChIP-seq sits higher than RNA-seq because these libraries are predominantly broad
// histone marks, which need more depth to separate domain enrichment from background
// than a point-source factor does.
@Field final Map FIXED_FRAGMENT_TARGETS = [
    RNASeq : 20000000L,
    ChipSeq: 30000000L,
]

@Field final List VALID_ASSAY_TYPES = ["DNASeq", "RNASeq", "ChipSeq"]

def targetOnTargetFragments(Map policy) {
    if (!VALID_ASSAY_TYPES.contains(policy.assayType)) {
        throw new IllegalArgumentException(
            "Unrecognized assayType '${policy.assayType}'; valid values are ${VALID_ASSAY_TYPES}"
        )
    }
    if (FIXED_FRAGMENT_TARGETS.containsKey(policy.assayType)) {
        return FIXED_FRAGMENT_TARGETS[policy.assayType]
    }
    // Truncate rather than round: a deliberately conservative (never over-) estimate.
    long raw = (long) ((policy.genomeSize * policy.targetCoverage) / policy.basesPerFragment)
    return Math.max(MIN_TARGET_FRAGMENTS, Math.min(MAX_TARGET_FRAGMENTS, raw))
}

def depthPlan(Map metrics, Map policy) {
    if (metrics.readLength == null || (metrics.readLength as double) <= 0) {
        throw new IllegalArgumentException(
            "Invalid metrics.readLength '${metrics.readLength}'; must be greater than 0"
        )
    }
    // No default: silently assuming single-end is exactly the bug this field exists to
    // prevent, and it would halve reported coverage on every paired library.
    if (metrics.mateCount == null || !([1, 2].contains(metrics.mateCount as int))) {
        throw new IllegalArgumentException(
            "Invalid metrics.mateCount '${metrics.mateCount}'; must be 1 (single-end) or 2 (paired)"
        )
    }
    if (metrics.totalReads == null || (metrics.totalReads as double) <= 0) {
        throw new IllegalArgumentException(
            "Invalid metrics.totalReads '${metrics.totalReads}'; must be greater than 0"
        )
    }
    if (metrics.onTargetFraction == null || (metrics.onTargetFraction as double) < 0.0d || (metrics.onTargetFraction as double) > 1.0d) {
        throw new IllegalArgumentException(
            "Invalid metrics.onTargetFraction '${metrics.onTargetFraction}'; must be within [0.0, 1.0]"
        )
    }
    if (policy.minOnTargetFraction == null || (policy.minOnTargetFraction as double) <= 0.0d) {
        throw new IllegalArgumentException(
            "Invalid policy.minOnTargetFraction '${policy.minOnTargetFraction}'; must be greater than 0.0"
        )
    }

    long basesPerFragment = (metrics.readLength as long) * (metrics.mateCount as long)
    def policyWithBases = policy + [basesPerFragment: basesPerFragment]
    long targetOnTarget = targetOnTargetFragments(policyWithBases)

    double observed = metrics.onTargetFraction as double
    double effective = Math.max(observed, policy.minOnTargetFraction as double)

    long rawFragments = Math.min(
        (long) Math.ceil(targetOnTarget / effective),
        metrics.totalReads as long
    )

    Double estimatedCoverage = FIXED_FRAGMENT_TARGETS.containsKey(policy.assayType)
        ? null
        : (rawFragments * observed * basesPerFragment) / (policy.genomeSize as double)

    return [
        targetOnTargetFragments: targetOnTarget,
        effectiveFraction      : effective,
        rawFragments           : rawFragments,
        basesPerFragment       : basesPerFragment,
        flagged                : observed < (policy.minPlausibleFraction as double),
        estimatedCoverage      : estimatedCoverage
    ]
}
