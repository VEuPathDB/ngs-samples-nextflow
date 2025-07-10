#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RETRIEVE_FROM_SRA       } from './workflows/retrieve_from_sra'

include { FORMAT_INPUT_FROM_SRA   } from './modules/local/format_input_from_sra'

workflow {

    samples = Channel.fromPath(params.input + "/" + params.samplesheetName )
        .splitCsv( skip:1)

    if(params.fromSra) {
        RETRIEVE_FROM_SRA(samples.map { row ->
            return [ [id: row[0], var1: row[3] ], row[1]]
        })
    }
    else {
        FORMAT_INPUT_FROM_SRA(samplesFullPath = samples.map { row ->
            fasta1 = file(params.input + "/" + row[1]);

            if(row[2]) {
                fasta2 = file(params.input + "/" + row[2])
                return [ [id: row[0], var1: row[3] ], [fasta1, fasta2] ]
            }
            return [ [id: row[0], var1: row[3] ], [fasta1] ]

            })
        formatted = FORMAT_INPUT_FROM_SRA.out.samplesheet.collectFile(keepHeader: true, storeDir: params.outDir, name: params.samplesheetName)

    }

}
