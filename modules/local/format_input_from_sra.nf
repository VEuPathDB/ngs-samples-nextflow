process FORMAT_INPUT_FROM_SRA {
    tag "$meta.id"
    label 'process_single'

    shell '/bin/sh'

    container 'docker.io/veupathdb/alpine_bash:1.0.0'

    publishDir params.outDir, mode: 'copy',  saveAs: { file ->
        file == "formattedInput.csv" ? null : file 
    }

    input:
    tuple val(meta), path(sra)

    output:
    path('formattedInput.csv'), emit: samplesheet
    path(sra), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def var1 = ""
    if(meta.var1) {
        var1 = meta.var1
    }
    if(params.fromSra) {
        """
        echo "sample,fastq_1,fastq_2,var1" > formattedInput.csv
        file_count=\$(ls -l ${sra} | grep -v '^d' | wc -l)
        if [ \$file_count == 1 ]; then
            echo "${meta.id},${params.outDir}/${meta.id}_1.fastq.gz,${params.outDir}/${meta.id}_2.fastq.gz,${var1}" >> formattedInput.csv
        else
            echo "${meta.id},${params.outDir}/${meta.id}_1.fastq.gz,,${var1}" >> formattedInput.csv
        fi
        """
    }
    else {
        """
        echo "sample,fastq_1,fastq_2,var1" > formattedInput.csv
        file_count=\$(ls -l ${sra} | grep -v '^d' | wc -l)
        if [ \$file_count == 1 ]; then
            echo "${meta.id},${params.outDir}/${sra[0]},,${var1}" >> formattedInput.csv
        else
            echo "${meta.id},${params.outDir}/${sra[0]},${params.outDir}/${sra[1]},${var1}" >> formattedInput.csv
        fi
        """
    }
    stub:
    """
    echo "sample,fastq_1,fastq_2,var1" > formattedInput.csv
    """
}
