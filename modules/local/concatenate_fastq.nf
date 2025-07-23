process CONCATENATE_FASTQ {
    tag "$meta.id"
    label 'process_low'

    shell '/bin/sh'

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
    def has_paired_files = file_list.any { it.name.contains('_1.fastq') || it.name.contains('_2.fastq') || it.name.contains('_R1') || it.name.contains('_R2') }
    
    if (!has_paired_files || file_list.size() == 1) {
        // Single-end case or single file
        """        
        cat ${file_list.join(' ')} > ${meta.id}.fastq.gz
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
                # Default to R1 if we can't determine
                r1_files+=("\$file")
            fi
        done
        
        # Concatenate R1 files if any exist
        if [ \${#r1_files[@]} -gt 0 ]; then
            cat "\${r1_files[@]}" > ${meta.id}_1.fastq.gz
        fi
        
        # Concatenate R2 files if any exist  
        if [ \${#r2_files[@]} -gt 0 ]; then
            cat "\${r2_files[@]}" > ${meta.id}_2.fastq.gz
        fi
        
        # If we only have R1 files, create an empty R2 to maintain pairing
        if [ \${#r1_files[@]} -gt 0 ] && [ \${#r2_files[@]} -eq 0 ]; then
            # This is actually single-end, rename accordingly
            mv ${meta.id}_1.fastq.gz ${meta.id}.fastq.gz
        fi
        """
    }

    stub:
    def file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def has_paired_files = file_list.any { it.name.contains('_1.fastq') || it.name.contains('_2.fastq') || it.name.contains('_R1') || it.name.contains('_R2') }
    
    if (!has_paired_files || file_list.size() == 1) {
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