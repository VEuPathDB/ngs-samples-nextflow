/*
 * Pure helpers for read layout. Downstream, mates are identified by position only:
 * R1 files and R2 files travel as separate, index-aligned lists.
 */

/*
 * Nextflow hands a `path('*.fastq.gz')` output back as a bare Path when the glob matches
 * one file and as a List when it matches several, and Path.size() is the file size in
 * bytes, so normalise before counting. fasterq-dump always names mates <acc>_1 / <acc>_2,
 * so sorting by name is the one place a filename reliably encodes mate order.
 */
def normalizeReads(reads, String sampleId = "unknown") {
    def files = reads instanceof Collection ? reads.toList() : [reads]
    if (files.size() != 1 && files.size() != 2) {
        throw new IllegalStateException("Sample ${sampleId}: fasterq-dump produced ${files.size()} files; expected 1 (single-end) or 2 (paired). A 3-file split indicates unpaired reads mixed with pairs, which this pipeline cannot currently concatenate safely.")
    }
    return files.sort { it.name }
}

/*
 * Merges the runs of one sample (the output of groupTuple) into [ meta, r1Files, r2Files ].
 * Paired-ness is derived here, once, from whether runs carry R2 files.
 */
def mergeRuns(List metas, List r1Lists, List r2Lists) {
    def sampleId = metas[0].id
    def pairedFlags = r2Lists.collect { !it.isEmpty() }.toSet()
    if (pairedFlags.size() > 1) {
        throw new IllegalStateException("Sample ${sampleId}: runs must be all paired or all single-end; found a mix.")
    }
    def meta = metas[0] + [hasPairedReads: pairedFlags.first()]
    return [ meta, r1Lists.collectMany { it }, r2Lists.collectMany { it } ]
}
