include { FORMAT_INPUT_FROM_SRA } from '../modules/local/format_input_from_sra'
include { CONCATENATE_FASTQ     } from '../modules/local/concatenate_fastq'
include { SUBSAMPLE_FASTQ       } from '../modules/local/subsample_fastq'
include { MEASURE_SAMPLE        } from '../modules/local/measure_sample'
include { depthPlan             } from '../modules/local/depth_policy'

workflow PREPARE_SAMPLES {

    take:
    grouped_reads      // [ meta, [files] ]
    reference_sig      // path (value channel)
    policy             // map (value channel)

    main:
    CONCATENATE_FASTQ(grouped_reads)

    MEASURE_SAMPLE(CONCATENATE_FASTQ.out.reads, reference_sig)

    // Join on sample id rather than on the meta map, so the join key stays stable
    // even if a meta field is mutated upstream.
    reads_by_id = CONCATENATE_FASTQ.out.reads.map { meta, reads -> [ meta.id, meta, reads ] }
    metrics_by_id = MEASURE_SAMPLE.out.metrics.map { meta, json -> [ meta.id, json ] }

    // `policy` arrives as a value channel (genome size is measured, not known up front),
    // so it must be combined into the stream rather than dereferenced directly.
    plans = reads_by_id
        .join(metrics_by_id)
        .combine(policy)
        .map { id, meta, reads, json, policyMap ->
            def metrics = new groovy.json.JsonSlurper().parse(json.toFile())
            def plan = depthPlan(metrics, policyMap)
            if (plan.flagged) {
                log.warn "Sample ${id}: on-target fraction ${metrics.onTargetFraction} is below " +
                         "minPlausibleFraction (${policyMap.minPlausibleFraction}). If every sample " +
                         "is flagged, check that --referenceFasta is the right organism."
            }
            return [ meta, reads, metrics, plan ]
        }

    SUBSAMPLE_FASTQ(
        plans.map { meta, reads, metrics, plan -> [ meta, reads, plan.rawReads, metrics.totalReads ] }
    )

    FORMAT_INPUT_FROM_SRA(SUBSAMPLE_FASTQ.out.reads)

    formatted = FORMAT_INPUT_FROM_SRA.out.samplesheet
        .collectFile(keepHeader: true, storeDir: params.outDir, name: params.samplesheetName)

    metrics_header = "sample,on_target_fraction,total_reads,raw_reads_used," +
                     "estimated_coverage,read_length,pilot_reads,flagged\n"

    metrics_csv = plans
        .map { meta, reads, metrics, plan ->
            def coverage = plan.estimatedCoverage == null
                ? ''
                : String.format('%.2f', plan.estimatedCoverage)
            "${meta.id},${metrics.onTargetFraction},${metrics.totalReads},${plan.rawReads}," +
            "${coverage},${metrics.readLength},${metrics.pilotReads},${plan.flagged}\n"
        }
        .collectFile(name: 'sample_metrics.csv', storeDir: params.outDir,
                     seed: metrics_header, sort: true)

    emit:
    formattedInput = formatted
    sampleMetrics  = metrics_csv
    flags          = plans.map { meta, reads, metrics, plan -> plan.flagged }
}
