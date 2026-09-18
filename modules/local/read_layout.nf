/*
 * Pure helper for the fasterq-dump output. Nextflow hands a `path('*.fastq.gz')` output
 * back as a bare Path when the glob matches one file and as a List when it matches
 * several, and Path.size() is the file size in bytes, so the layout check must
 * normalise before counting.
 */

def normalizeReads(reads, String sampleId = "unknown") {
    def files = reads instanceof Collection ? reads.toList() : [reads]
    if (files.size() != 1 && files.size() != 2) {
        throw new IllegalStateException("Sample ${sampleId}: fasterq-dump produced ${files.size()} files; expected 1 (single-end) or 2 (paired). A 3-file split indicates unpaired reads mixed with pairs, which this pipeline cannot currently concatenate safely.")
    }
    return files
}
