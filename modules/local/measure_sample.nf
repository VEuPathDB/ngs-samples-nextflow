process MEASURE_SAMPLE {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0'

    input:
    tuple val(meta), path(reads)
    path reference_sig

    output:
    tuple val(meta), path("${meta.id}.metrics.json"), emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def read1 = reads instanceof List ? reads[0] : reads
    def pilot_size = params.pilotSize
    """
    set -euo pipefail

    total_reads=\$(( \$(zcat -f ${read1} | wc -l) / 4 ))
    read_length=\$(zcat -f ${read1} | sed -n '2p' | tr -d '\\n' | wc -c)

    step=\$(( total_reads / ${pilot_size} ))
    if [ "\$step" -lt 1 ]; then step=1; fi

    zcat -f ${read1} \\
      | awk -v step="\$step" 'NR%4==1{keep=(int((NR-1)/4)%step==0)} keep' > pilot.fastq
    pilot_reads=\$(( \$(wc -l < pilot.fastq) / 4 ))

    sourmash sketch dna -p k=31,scaled=1000,abund --name '${meta.id}' pilot.fastq -o pilot.sig

    # --threshold-bp 0 is required: the default 50kbp threshold makes gather write no result
    # row for low-overlap samples, which is exactly the case this measurement exists for.
    sourmash gather pilot.sig ${reference_sig} --threshold-bp 0 -o gather.csv || true

    TOTAL_READS=\$total_reads READ_LENGTH=\$read_length PILOT_READS=\$pilot_reads \\
    OUT='${meta.id}.metrics.json' python3 -c '
import csv, json, os
rows = []
if os.path.exists("gather.csv"):
    with open("gather.csv") as fh:
        rows = list(csv.DictReader(fh))
fraction = float(rows[0]["f_unique_weighted"]) if rows else 0.0
json.dump({"onTargetFraction": fraction,
           "totalReads": int(os.environ["TOTAL_READS"]),
           "readLength": int(os.environ["READ_LENGTH"]),
           "pilotReads": int(os.environ["PILOT_READS"])},
          open(os.environ["OUT"], "w"), indent=2)
'
    """

    stub:
    """
    echo '{"onTargetFraction": 0.5, "totalReads": 1000000, "readLength": 150, "pilotReads": 100000}' > ${meta.id}.metrics.json
    """
}
