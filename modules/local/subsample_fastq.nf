process SUBSAMPLE_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    container 'staphb/seqtk:1.4'

    publishDir params.outDir, mode: 'copy'

    // Both counts are fragments (R1 records), so their ratio is the per-file sampling
    // fraction for single- and paired-end alike. One seed keeps mates in step.
    // format ("fastq" or "fasta", as detected by MEASURE_SAMPLE) sets the output extension.
    input:
    tuple val(meta), path(reads), val(target_fragments), val(total_fragments), val(format)

    output:
    tuple val(meta), path("${meta.id}*_subsampled.${format}.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def seed = 42
    def read_list = reads instanceof List ? reads : [reads]

    if (meta.hasPairedReads) {
        """
        if [ ${total_fragments} -gt ${target_fragments} ]; then
            fraction=\$(awk -v t="${target_fragments}" -v n="${total_fragments}" 'BEGIN {printf "%.10f", t/n}')
            seqtk sample -s ${seed} ${read_list[0]} \$fraction | gzip > ${meta.id}_1_subsampled.${format}.gz
            seqtk sample -s ${seed} ${read_list[1]} \$fraction | gzip > ${meta.id}_2_subsampled.${format}.gz
        else
            ln -s ${read_list[0]} ${meta.id}_1_subsampled.${format}.gz
            ln -s ${read_list[1]} ${meta.id}_2_subsampled.${format}.gz
        fi
        """
    } else {
        """
        if [ ${total_fragments} -gt ${target_fragments} ]; then
            fraction=\$(awk -v t="${target_fragments}" -v n="${total_fragments}" 'BEGIN {printf "%.10f", t/n}')
            seqtk sample -s ${seed} ${read_list[0]} \$fraction | gzip > ${meta.id}_subsampled.${format}.gz
        else
            ln -s ${read_list[0]} ${meta.id}_subsampled.${format}.gz
        fi
        """
    }

    stub:
    if (meta.hasPairedReads) {
        """
        touch ${meta.id}_1_subsampled.${format}.gz
        touch ${meta.id}_2_subsampled.${format}.gz
        """
    } else {
        """
        touch ${meta.id}_subsampled.${format}.gz
        """
    }
}
