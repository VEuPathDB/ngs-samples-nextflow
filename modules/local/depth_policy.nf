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

// Each assay sets its target on a different axis, so none of these need pairing
// adjustment beyond basesPerFragment below.

// RNA-seq depth tracks transcriptome complexity - gene count and expression dynamic
// range - which moves far less than genome size across these organisms, so it is flat.
@Field final long RNASEQ_TARGET_FRAGMENTS = 20000000L

// ChIP-seq depth does scale with genome size, but sublinearly: reads-per-peak is
// genome-independent, peak count grows slowly, and only background depth is linear.
// modENCODE's worm/fly minimums and ENCODE's human ones independently imply an exponent
// near 0.5. The anchor sits above those published figures because they count uniquely
// mapped reads while this target is measured before alignment and duplicate removal.
// The floor only guards the smallest genomes: it binds below ~24Mb, where it and the curve
// agree closely, so depth is set by the curve for essentially every organism here.
@Field final long   CHIPSEQ_FLOOR_FRAGMENTS  =  3000000L
@Field final long   CHIPSEQ_ANCHOR_FRAGMENTS =  7500000L
@Field final long   CHIPSEQ_ANCHOR_GENOME    = 150000000L
@Field final double CHIPSEQ_GENOME_EXPONENT  = 0.5d

@Field final List VALID_ASSAY_TYPES = ["DNASeq", "RNASeq", "ChipSeq"]

def targetOnTargetFragments(Map policy) {
    if (!VALID_ASSAY_TYPES.contains(policy.assayType)) {
        throw new IllegalArgumentException(
            "Unrecognized assayType '${policy.assayType}'; valid values are ${VALID_ASSAY_TYPES}"
        )
    }
    if (policy.assayType == "RNASeq") {
        return RNASEQ_TARGET_FRAGMENTS
    }
    if (policy.assayType == "ChipSeq") {
        double scaled = CHIPSEQ_ANCHOR_FRAGMENTS * Math.pow(
            (policy.genomeSize as double) / CHIPSEQ_ANCHOR_GENOME, CHIPSEQ_GENOME_EXPONENT)
        return Math.min(MAX_TARGET_FRAGMENTS, Math.max(CHIPSEQ_FLOOR_FRAGMENTS, (long) scaled))
    }
    // DNASeq, the only coverage-denominated assay. Truncate rather than round: a
    // deliberately conservative (never over-) estimate.
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

    Double estimatedCoverage = policy.assayType == "DNASeq"
        ? (rawFragments * observed * basesPerFragment) / (policy.genomeSize as double)
        : null

    return [
        targetOnTargetFragments: targetOnTarget,
        effectiveFraction      : effective,
        rawFragments           : rawFragments,
        basesPerFragment       : basesPerFragment,
        flagged                : observed < (policy.minPlausibleFraction as double),
        estimatedCoverage      : estimatedCoverage
    ]
}
