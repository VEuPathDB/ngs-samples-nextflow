/*
 * Pure depth-policy functions. No params access, no file IO — everything is passed in
 * so this can be unit tested without a pipeline run.
 */

def targetOnTargetReads(Map policy) {
    if (policy.assayType == "RNASeq") {
        return 20000000L
    }
    long raw = (long) ((policy.genomeSize * policy.targetCoverage) / policy.readLength)
    return Math.max(1000000L, Math.min(100000000L, raw))
}

def depthPlan(Map metrics, Map policy) {
    def resolved = policy + [readLength: metrics.readLength]
    long targetOnTarget = targetOnTargetReads(resolved)

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
