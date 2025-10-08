#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def calculateMaxReads(assayType, genomeSize) {
    // If maxReads is explicitly set, use that value
    if (assayType == "RNASeq") {
        int maxReads = 20000000;
        return maxReads.toLong()
    }
    
    def genomeSizeLong = genomeSize.toLong()
    def targetCoverage
    def readLength = 150  // Assume 150bp reads
    
    // Set target coverage based on assay type
    switch(assayType) {
        case "DNASeq":
            targetCoverage = 60  // 60x coverage for DNA-seq
            break
        case "ChipSeq":
            targetCoverage = 60  // 60x coverage for Chip-seq
            break
        default:
            targetCoverage = 60  // Default to DNA-seq coverage
            break
    }
    
    // Calculate number of reads needed: (genome_size * coverage) / read_length
    def maxReads = (genomeSizeLong * targetCoverage) / readLength
    
    // Set reasonable bounds (min 1M reads, max 100M reads)
    maxReads = Math.max(1000000L, Math.min(100000000L, maxReads.toLong()))
    
    log.info "Calculated max reads for ${assayType} (${genomeSizeLong}bp genome, ${targetCoverage}x coverage): ${maxReads}"
    
    return maxReads
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RETRIEVE_FROM_SRA       } from './workflows/retrieve_from_sra'

include { FORMAT_INPUT_FROM_SRA   } from './modules/local/format_input_from_sra'
include { CONCATENATE_FASTQ       } from './modules/local/concatenate_fastq'
include { SUBSAMPLE_FASTQ         } from './modules/local/subsample_fastq'

workflow {

    samples = Channel.fromPath(params.input + "/" + params.samplesheetName )
        .splitCsv( skip:1)
    
    // Calculate maximum reads for subsampling
    max_reads = calculateMaxReads(params.assayType, params.genomeSize)

    if(params.fromSra) {
        // Group SRA samples by ID and collect all SRA accessions per sample
        grouped_sra_samples = samples.map { row ->
            return [ row[0], [id: row[0], var1: row[3] ], row[1]]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, sra_ids ->
            // Use the first meta (they should be the same except for SRA IDs)
            return [ metas[0], sra_ids ]
        }

        RETRIEVE_FROM_SRA(grouped_sra_samples, max_reads)
    }
    else {
        // Group local file samples by ID and collect all files per sample
        grouped_local_samples = samples.map { row ->
            fasta1 = file(params.input + "/" + row[1]);
            files = [fasta1]
            
            if(row[2]) {
                fasta2 = file(params.input + "/" + row[2])
                files.add(fasta2)
            }
            
            return [ row[0], [id: row[0], var1: row[3] ], files ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, file_lists ->
            // Flatten the file lists and use the first meta
            all_files = file_lists.flatten()
            return [ metas[0], all_files ]
        }
        
        // Concatenate files if multiple entries per sample
        CONCATENATE_FASTQ(grouped_local_samples)
        
        // Subsample reads to limit coverage
        SUBSAMPLE_FASTQ(CONCATENATE_FASTQ.out.reads, max_reads)
        
        FORMAT_INPUT_FROM_SRA(SUBSAMPLE_FASTQ.out.reads)
        formatted = FORMAT_INPUT_FROM_SRA.out.samplesheet.collectFile(keepHeader: true, storeDir: params.outDir, name: params.samplesheetName)
    }

}
