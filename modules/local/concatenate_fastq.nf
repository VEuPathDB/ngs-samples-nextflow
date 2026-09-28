process CONCATENATE_FASTQ {
    tag "$meta.id"
    label 'process_low'

    shell '/bin/bash'

    container 'docker.io/veupathdb/alpine_bash:1.0.0'

    input:
    // Mates are identified by position: r1_files[i] pairs with r2_files[i]. Filenames are
    // never inspected, since naming conventions vary too widely to infer mate from.
    tuple val(meta), path(r1_files, stageAs: "r1/?/*"), path(r2_files, stageAs: "r2/?/*")

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def asList = { it instanceof Collection ? it.toList() : (it ? [it] : []) }
    def r1 = asList(r1_files)
    def r2 = asList(r2_files)

    def layoutError = null
    if (r1.isEmpty()) {
        layoutError = "no R1 files"
    } else if (meta.hasPairedReads && r2.size() != r1.size()) {
        layoutError = "${r1.size()} R1 file(s) but ${r2.size()} R2 file(s)"
    } else if (!meta.hasPairedReads && r2) {
        layoutError = "single-end sample was given ${r2.size()} R2 file(s)"
    }
    if (layoutError) {
        return """
        echo "CONCATENATE_FASTQ ${meta.id}: ${layoutError}" >&2
        exit 1
        """
    }

    def mates = meta.hasPairedReads
        ? [ [r1, "${meta.id}_concat_1.fastq.gz"], [r2, "${meta.id}_concat_2.fastq.gz"] ]
        : [ [r1, "${meta.id}_concat.fastq.gz"] ]

    // A single already-gzipped run needs no merging; symlinking avoids a pointless
    // decompress/recompress cycle at the pipeline's peak-disk step.
    def commands = mates.collect { files, target ->
        files.size() == 1 && files[0].name.endsWith('.gz')
            ? "ln -s ${files[0]} ${target}"
            : "{ for f in ${files.join(' ')}; do read_fastq \"\$f\"; done } | gzip > ${target}"
    }

    """
    read_fastq() { [[ \$(xxd -l 2 "\$1" | awk '{print \$2\$3}') == "1f8b"* ]] && zcat "\$1" || cat "\$1"; }
    ${commands.join('\n    ')}
    """

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
