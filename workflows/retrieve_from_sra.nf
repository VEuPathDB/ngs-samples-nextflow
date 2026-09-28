include { SRATOOLS_FASTERQDUMP } from '../modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH    } from '../modules/nf-core/sratools/prefetch/main'
include { PREPARE_SAMPLES      } from './prepare_samples'
include { normalizeReads; mergeRuns } from '../modules/local/read_layout'

workflow RETRIEVE_FROM_SRA {

    take:
    samples
    reference_sig
    policy

    main:
    individual_sra_samples = samples.flatMap { meta, sra_ids ->
        sra_ids.collect { sra_id -> [ meta, sra_id ] }
    }

    SRATOOLS_PREFETCH(individual_sra_samples, [], [])
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra, [], [])

    grouped_reads = SRATOOLS_FASTERQDUMP.out.reads
        .map { meta, reads ->
            def files = normalizeReads(reads, meta.id)
            return [ meta.id, meta, files.take(1), files.drop(1) ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, r1_lists, r2_lists -> mergeRuns(metas, r1_lists, r2_lists) }

    PREPARE_SAMPLES(grouped_reads, reference_sig, policy)

    emit:
    formattedInput = PREPARE_SAMPLES.out.formattedInput
}
