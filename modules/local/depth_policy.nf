import groovy.transform.Field

/*
 * Pure depth-policy functions. No params access, no file IO — everything is passed in
 * so this can be unit tested without a pipeline run.
 */

// Cost/runtime sanity bounds: a floor so tiny genomes still get usable depth,
// a ceiling so huge genomes don't blow up downstream alignment cost.
@Field final long MIN_TARGET_READS = 1000000L
@Field final long MAX_TARGET_READS = 100000000L
@Field final long RNASEQ_TARGET_READS = 20000000L

@Field final List VALID_ASSAY_TYPES = ["DNASeq", "RNASeq", "ChipSeq"]

def targetOnTargetReads(Map policy) {
    if (!VALID_ASSAY_TYPES.contains(policy.assayType)) {
        throw new IllegalArgumentException(
            "Unrecognized assayType '${policy.assayType}'; valid values are ${VALID_ASSAY_TYPES}"
        )
    }
    if (policy.assayType == "RNASeq") {
        return RNASEQ_TARGET_READS
    }
    // Truncate rather than round: a deliberately conservative (never over-) estimate.
    long raw = (long) ((policy.genomeSize * policy.targetCoverage) / policy.readLength)
    return Math.max(MIN_TARGET_READS, Math.min(MAX_TARGET_READS, raw))
}

def depthPlan(Map metrics, Map policy) {
    if (metrics.readLength == null || (metrics.readLength as double) <= 0) {
        throw new IllegalArgumentException(
            "Invalid metrics.readLength '${metrics.readLength}'; must be greater than 0"
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

    def policyWithReadLength = policy + [readLength: metrics.readLength]
    long targetOnTarget = targetOnTargetReads(policyWithReadLength)

    double observed = metrics.onTargetFraction as double
    double effective = Math.max(observed, policy.minOnTargetFraction as double)

    long rawReads = Math.min(
        (long) Math.ceil(targetOnTarget / effective),
        metrics.totalReads as long
    )

    Double estimatedCoverage = policy.assayType == "RNASeq"
        ? null
        : (rawReads * observed * metrics.readLength) / (policy.genomeSize as double)

    return [
        targetOnTarget   : targetOnTarget,
        effectiveFraction: effective,
        rawReads         : rawReads,
        flagged          : observed < (policy.minPlausibleFraction as double),
        estimatedCoverage: estimatedCoverage
    ]
}
