#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RETRIEVE_FROM_SRA } from './workflows/retrieve_from_sra'
include { PREPARE_SAMPLES   } from './workflows/prepare_samples'
include { SKETCH_REFERENCE  } from './modules/local/sketch_reference'
include { EXPAND_SRX_IDS    } from './modules/local/expand_srx_ids'
include { mergeRuns         } from './modules/local/read_layout'

workflow {

    if (!params.referenceFasta) {
        error "--referenceFasta is required: the target organism FASTA used to estimate " +
              "on-target fraction. Subsampling cannot be coverage-aware without it."
    }
    def referenceFile = file(params.referenceFasta, checkIfExists: true)
    if (referenceFile.size() == 0) {
        error "--referenceFasta (${params.referenceFasta}) is empty."
    }

    SKETCH_REFERENCE(Channel.value(referenceFile))

    policy = SKETCH_REFERENCE.out.stats.map { statsFile ->
        def stats = new groovy.json.JsonSlurper().parse(statsFile.toFile())
        return [
            assayType           : params.assayType,
            genomeSize          : stats.genomeSize as long,
            targetCoverage      : params.targetCoverage as int,
            minOnTargetFraction : params.minOnTargetFraction as double,
            minPlausibleFraction: params.minPlausibleFraction as double
        ]
    }

    reference_sig = SKETCH_REFERENCE.out.sig

    samples = Channel.fromPath(params.input + "/" + params.samplesheetName).splitCsv(skip: 1)

    if (params.fromSra) {
        EXPAND_SRX_IDS(samples.map { row -> [ row[0], row[1], row[3] ?: "" ] })

        grouped_sra_samples = EXPAND_SRX_IDS.out.rows
            .splitCsv()
            .map { row -> [ row[0], [id: row[0], var1: row[3] ?: ""], row[1] ] }
            .groupTuple(by: 0)
            .map { sample_id, metas, sra_ids -> [ metas[0], sra_ids ] }

        RETRIEVE_FROM_SRA(grouped_sra_samples, reference_sig, policy)
    }
    else {
        grouped_local_samples = samples
            .map { row ->
                def r1 = [ file(params.input + "/" + row[1], checkIfExists: true) ]
                def r2 = row[2] ? [ file(params.input + "/" + row[2], checkIfExists: true) ] : []
                return [ row[0], [id: row[0], var1: row[3]], r1, r2 ]
            }
            .groupTuple(by: 0)
            .map { sample_id, metas, r1_lists, r2_lists -> mergeRuns(metas, r1_lists, r2_lists) }

        PREPARE_SAMPLES(grouped_local_samples, reference_sig, policy)
    }
}
