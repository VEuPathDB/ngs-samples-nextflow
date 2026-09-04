process CONCATENATE_FASTQ {
    tag "$meta.id"
    label 'process_low'

    shell '/bin/bash'

    container 'docker.io/veupathdb/alpine_bash:1.0.0'

    input:
    tuple val(meta), path(fastq_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    // Check if we have only single files or if this is mixed single/paired
    def file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    //def has_paired_files = file_list.any { it.name.contains('_1.fastq') || it.name.contains('_2.fastq') || it.name.contains('_R1') || it.name.contains('_R2') }

    // A single already-gzipped run needs no merging; symlinking avoids a pointless
    // decompress/recompress cycle at the pipeline's peak-disk step. Guarded so a malformed
    // sample (paired flag but only one file) falls through to the real concatenation path.
    def all_gzipped = file_list.every { it.name.endsWith('.gz') }

    if (!meta.hasPairedReads && file_list.size() == 1 && all_gzipped) {
        return """
        ln -s ${file_list[0]} ${meta.id}_concat.fastq.gz
        """
    }

    if (meta.hasPairedReads && file_list.size() == 2 && all_gzipped) {
        def r1 = file_list.find { it.name.contains('_1.fastq') || it.name.contains('_R1') }
        def r2 = file_list.find { it.name.contains('_2.fastq') || it.name.contains('_R2') }
        if (r1 && r2 && r1 != r2) {
            return """
            ln -s ${r1} ${meta.id}_concat_1.fastq.gz
            ln -s ${r2} ${meta.id}_concat_2.fastq.gz
            """
        }
    }

    if (!meta.hasPairedReads) {
        // Single-end case
        """
        read_fastq() { [[ \$(xxd -l 2 "\$1" | awk '{print \$2\$3}') == "1f8b"* ]] && zcat "\$1" || cat "\$1"; }
        { for f in ${file_list.join(' ')}; do read_fastq "\$f"; done } | gzip > ${meta.id}_concat.fastq.gz
        """
    } else {
        // Paired-end: separate R1 and R2 files and concatenate each
        """
        read_fastq() { [[ \$(xxd -l 2 "\$1" | awk '{print \$2\$3}') == "1f8b"* ]] && zcat "\$1" || cat "\$1"; }

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
            { for f in "\${r1_files[@]}"; do read_fastq "\$f"; done } | gzip > ${meta.id}_concat_1.fastq.gz
        fi

        # Concatenate R2 files if any exist
        if [ \${#r2_files[@]} -gt 0 ]; then
            { for f in "\${r2_files[@]}"; do read_fastq "\$f"; done } | gzip > ${meta.id}_concat_2.fastq.gz
        fi
        """
    }

    stub:
    if (!meta.hasPairedReads) {
        """
        touch ${meta.id}_concat.fastq.gz
        """
    } else {
        """
        touch ${meta.id}_concat_1.fastq.gz
        touch ${meta.id}_concat_2.fastq.gz
        """
    }
}
