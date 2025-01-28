/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SRATOOLS_FASTERQDUMP                      } from '../../../modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH                         } from '../../../modules/nf-core/sratools/prefetch/main'
include { FORMAT_INPUT_FROM_SRA                     } from '../../../modules/local/format_input_from_sra'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow RETRIEVE_FROM_SRA {

    take:
    input

    main:

    samples = Channel.fromPath(input)
                     .splitCsv()
		     .map { row ->
		            return [[id: row[0]], row[1] ]
	             }

    SRATOOLS_PREFETCH(samples,[],[]);
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra,[],[]);
    FORMAT_INPUT_FROM_SRA(SRATOOLS_FASTERQDUMP.out.reads);
    formatted = FORMAT_INPUT_FROM_SRA.out.collectFile(keepHeader: true, storeDir: params.outdir, name: "formattedSraInput.csv")

    emit:
    formattedInput = formatted
}
