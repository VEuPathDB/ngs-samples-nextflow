process SUBSAMPLE_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    container 'staphb/seqtk:1.4'

    publishDir params.outDir, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val max_reads

    output:
    tuple val(meta), path("${meta.id}*_subsampled.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_paired = reads instanceof List && reads.size() == 2
    def seed = 42  // Fixed seed for reproducibility
    
    if (is_paired) {
        def read1 = reads[0]
        def read2 = reads[1]
        """
        # Count reads in first file to determine if subsampling is needed
        total_reads=\$(( \$(zcat ${read1} | wc -l) / 4 ))
        
        if [ \$total_reads -gt ${max_reads} ]; then
            echo "Subsampling paired-end reads from \$total_reads to ${max_reads}"
            
            # Calculate sampling fraction
            fraction=\$(echo "scale=10; ${max_reads} / \$total_reads" | bc -l)
            
            # Subsample both files with same seed to maintain pairing
            seqtk sample -s ${seed} ${read1} \$fraction | gzip > ${meta.id}_1_subsampled.fastq.gz
            seqtk sample -s ${seed} ${read2} \$fraction | gzip > ${meta.id}_2_subsampled.fastq.gz
        else
            echo "No subsampling needed. Total reads (\$total_reads) <= max reads (${max_reads})"
            # Create symlinks with subsampled naming for consistency
            ln -s ${read1} ${meta.id}_1_subsampled.fastq.gz
            ln -s ${read2} ${meta.id}_2_subsampled.fastq.gz
        fi
        """
    } else {
        def read_file = reads instanceof List ? reads[0] : reads
        """
        # Count reads to determine if subsampling is needed
        total_reads=\$(( \$(zcat ${read_file} | wc -l) / 4 ))
        
        if [ \$total_reads -gt ${max_reads} ]; then
            echo "Subsampling single-end reads from \$total_reads to ${max_reads}"
            
            # Calculate sampling fraction
            fraction=\$(echo "scale=10; ${max_reads} / \$total_reads" | bc -l)
            
            # Subsample the file
            seqtk sample -s ${seed} ${read_file} \$fraction | gzip > ${meta.id}_subsampled.fastq.gz
        else
            echo "No subsampling needed. Total reads (\$total_reads) <= max reads (${max_reads})"
            # Create symlink with subsampled naming for consistency
            ln -s ${read_file} ${meta.id}_subsampled.fastq.gz
        fi
        """
    }

    stub:
    def is_paired = reads instanceof List && reads.size() == 2
    
    if (is_paired) {
        """
        touch ${meta.id}_1_subsampled.fastq.gz
        touch ${meta.id}_2_subsampled.fastq.gz
        """
    } else {
        """
        touch ${meta.id}_subsampled.fastq.gz
        """
    }
}