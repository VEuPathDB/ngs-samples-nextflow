process SUBSAMPLE_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    container 'staphb/seqtk:1.4'

    publishDir params.outDir, mode: 'copy'

    input:
    tuple val(meta), path(reads), val(target_reads), val(total_reads)

    output:
    tuple val(meta), path("${meta.id}*_subsampled.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def seed = 42
    def read_list = reads instanceof List ? reads : [reads]

    if (meta.hasPairedReads) {
        """
        if [ ${total_reads} -gt ${target_reads} ]; then
            fraction=\$(awk -v t="${target_reads}" -v n="${total_reads}" 'BEGIN {printf "%.10f", t/n}')
            seqtk sample -s ${seed} ${read_list[0]} \$fraction | gzip > ${meta.id}_1_subsampled.fastq.gz
            seqtk sample -s ${seed} ${read_list[1]} \$fraction | gzip > ${meta.id}_2_subsampled.fastq.gz
        else
            ln -s ${read_list[0]} ${meta.id}_1_subsampled.fastq.gz
            ln -s ${read_list[1]} ${meta.id}_2_subsampled.fastq.gz
        fi
        """
    } else {
        """
        if [ ${total_reads} -gt ${target_reads} ]; then
            fraction=\$(awk -v t="${target_reads}" -v n="${total_reads}" 'BEGIN {printf "%.10f", t/n}')
            seqtk sample -s ${seed} ${read_list[0]} \$fraction | gzip > ${meta.id}_subsampled.fastq.gz
        else
            ln -s ${read_list[0]} ${meta.id}_subsampled.fastq.gz
        fi
        """
    }

    stub:
    if (meta.hasPairedReads) {
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
