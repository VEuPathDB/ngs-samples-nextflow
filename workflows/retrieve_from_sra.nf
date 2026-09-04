include { SRATOOLS_FASTERQDUMP } from '../modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH    } from '../modules/nf-core/sratools/prefetch/main'
include { PREPARE_SAMPLES      } from './prepare_samples'

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
            meta.hasPairedReads = reads.size() == 2
            return [ meta.id, meta, reads ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, read_lists ->
            def firstHasPairedReads = metas[0].hasPairedReads
            if (!metas.every { it.hasPairedReads == firstHasPairedReads }) {
                throw new IllegalStateException("SRR samples must be all paired or unpaired. Mixed results found")
            }
            return [ metas[0], read_lists.flatten() ]
        }

    PREPARE_SAMPLES(grouped_reads, reference_sig, policy)

    emit:
    formattedInput = PREPARE_SAMPLES.out.formattedInput
    flags          = PREPARE_SAMPLES.out.flags
}
