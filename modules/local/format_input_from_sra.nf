process FORMAT_INPUT_FROM_SRA {
    tag "$meta.id"
    label 'process_single'

    shell '/bin/sh'

    container 'docker.io/veupathdb/alpine_bash:v1.0.0'

    publishDir params.outDir, mode: 'copy', pattern: "*fastq*"

    input:
    tuple val(meta), path(sra)

    output:
    path('formattedInput.csv'), emit: samplesheet
    path(sra), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "sample,fastq_1,fastq_2,strand" > formattedInput.csv
    if [ -f "${meta.id}_2.fastq.gz" ]; then
    echo "${meta.id},${params.outDir}/${meta.id}_1.fastq.gz,${params.outDir}/${meta.id}_2.fastq.gz,${meta.strand}" >> formattedInput.csv
    else
        echo "${meta.id},${params.outDir}/$sra,,${meta.strand}" >> formattedInput.csv
    fi
    """
    stub:
    """
    echo "sample,fastq_1,fastq_2,strand" > formattedInput.csv
    """
}
