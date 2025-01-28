process FORMAT_INPUT_FROM_SRA {
    tag "$meta.id"
    label 'process_single'

    shell '/bin/sh'

    container 'docker.io/veupathdb/alpine_bash:latest'

    input:
    tuple val(meta), path(sra)

    output:
    path('formattedInput.csv')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "sample,fastq_1,fastq_2" > formattedInput.csv
    if [ -f "${meta.id}_2.fastq.gz" ]; then
        echo "${meta.id},\$(pwd)/${meta.id}_1.fastq.gz,\$(pwd)/${meta.id}_2.fastq.gz" >> formattedInput.csv
    else
        echo "${meta.id},\$(pwd)/$sra" >> formattedInput.csv
    fi
    """
    stub:
    """
    echo "sample,fastq_1,fastq_2" > formattedInput.csv
    """
}
