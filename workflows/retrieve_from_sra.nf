/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SRATOOLS_FASTERQDUMP                      } from '../modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH                         } from '../modules/nf-core/sratools/prefetch/main'
include { FORMAT_INPUT_FROM_SRA                     } from '../modules/local/format_input_from_sra'
include { CONCATENATE_FASTQ                         } from '../modules/local/concatenate_fastq'
include { SUBSAMPLE_FASTQ                           } from '../modules/local/subsample_fastq'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow RETRIEVE_FROM_SRA {

    take:
    samples
    max_reads

    main:
    // Flatten samples to process each SRA ID individually
    individual_sra_samples = samples.flatMap { meta, sra_ids ->
        sra_ids.collect { sra_id ->
            [ meta, sra_id ]
        }
    }
    
    SRATOOLS_PREFETCH(individual_sra_samples,[],[]);
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra,[],[]);

    // Group back by sample ID and concatenate FASTQ files
    grouped_reads = SRATOOLS_FASTERQDUMP.out.reads
        .map { meta, reads ->
            [ meta.id, meta, reads ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, read_lists ->
            // Use the first meta and flatten all reads
            all_reads = read_lists.flatten()
            return [ metas[0], all_reads ]
        }
    
    CONCATENATE_FASTQ(grouped_reads)
    
    // Subsample reads to limit coverage
    SUBSAMPLE_FASTQ(CONCATENATE_FASTQ.out.reads, max_reads)
    
    FORMAT_INPUT_FROM_SRA(SUBSAMPLE_FASTQ.out.reads);
    formatted = FORMAT_INPUT_FROM_SRA.out.samplesheet.collectFile(keepHeader: true, storeDir: params.outDir, name: params.samplesheetName)

    emit:
    formattedInput = formatted
}
