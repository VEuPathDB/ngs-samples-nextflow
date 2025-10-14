process CONCATENATE_FASTQ {
    tag "$meta.id"
    label 'process_low'

    shell '/bin/bash'

    container 'docker.io/veupathdb/alpine_bash:1.0.0'

    publishDir params.outDir, mode: 'copy'

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    // Check if we have only single files or if this is mixed single/paired
    def file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    //def has_paired_files = file_list.any { it.name.contains('_1.fastq') || it.name.contains('_2.fastq') || it.name.contains('_R1') || it.name.contains('_R2') }
    
    if (!meta.hasPairedReads) {
        // Single-end case 
        """        
        zcat ${file_list.join(' ')} | gzip > ${meta.id}_concat.fastq.gz 
        """
    } else {
        // Paired-end: separate R1 and R2 files and concatenate each
        """
        # Create arrays to hold R1 and R2 files
        declare -a r1_files
        declare -a r2_files
        
        # Sort files into R1 and R2 arrays
        for file in ${file_list.join(' ')}; do
            basename_file=\$(basename "\$file")
            if [[ "\$basename_file" == *"_1.fastq"* ]] || [[ "\$basename_file" == *"_R1"* ]]; then
                r1_files+=("\$file")
            elif [[ "\$basename_file" == *"_2.fastq"* ]] || [[ "\$basename_file" == *"_R2"* ]]; then
                r2_files+=("\$file")
            else
               exit 1
            fi
        done
        
        # Concatenate R1 files if any exist
        if [ \${#r1_files[@]} -gt 0 ]; then
            zcat "\${r1_files[@]}" | gzip > ${meta.id}_concat_1.fastq.gz

        fi
        
        # Concatenate R2 files if any exist  
        if [ \${#r2_files[@]} -gt 0 ]; then
            zcat "\${r2_files[@]}"  | gzip > ${meta.id}_concat_2.fastq.gz
        fi
        """
    }

    stub:
    if (!hasPairedReads) {
        """
        touch ${meta.id}.fastq.gz
        """
    } else {
        """
        touch ${meta.id}_1.fastq.gz
        touch ${meta.id}_2.fastq.gz
        """
    }
}
