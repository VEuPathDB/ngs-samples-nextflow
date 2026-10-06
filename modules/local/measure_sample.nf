process MEASURE_SAMPLE {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0'

    input:
    tuple val(meta), path(reads)
    path reference_sig
    val min_fragments

    output:
    tuple val(meta), path("${meta.id}.metrics.json"), emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def read1 = reads instanceof List ? reads[0] : reads
    def pilot_size = params.pilotSize
    """
    set -euo pipefail

    # Single decompression pass: reservoir-sample up to pilot_size reads (Algorithm R)
    # while counting total records and accumulating sequence length for the mean, so we
    # never need to know total_reads before deciding what to keep. Below pilot_size every
    # read is kept (no randomness, no bias); above it, a uniform random subset spanning the
    # whole file is kept rather than just a fixed stride or the head of the file.
    # FASTA input is only counted: it is never subsampled, so it needs no pilot.
    stats=\$(zcat -f ${read1} | awk -v pilot_size="${pilot_size}" '
      BEGIN { srand(42); format = "fastq" }
      NR == 1 {
        c = substr(\$0, 1, 1)
        if (c == ">") format = "fasta"
        else if (c != "@") format = "unknown"
      }
      format == "unknown" { next }
      format == "fasta" {
        if (substr(\$0, 1, 1) == ">") total++
        else sumlen += length(\$0)
        next
      }
      {
        rec[(NR - 1) % 4] = \$0
        if ((NR - 1) % 4 == 3) {
          total++
          seqlen = length(rec[1])
          full = rec[0] "\\n" rec[1] "\\n" rec[2] "\\n" rec[3]
          if (total <= pilot_size) {
            reservoir[total] = full
            lens[total] = seqlen
            sumlen += seqlen
            count = total
          } else {
            j = int(rand() * total) + 1
            if (j <= pilot_size) {
              sumlen += seqlen - lens[j]
              reservoir[j] = full
              lens[j] = seqlen
            }
          }
        }
      }
      END {
        if (format == "unknown") { print format; exit }
        if (format == "fasta") {
          print format, total + 0, 0, (total > 0) ? int(sumlen / total + 0.5) : 0
          exit
        }
        printf "" > "pilot.fastq"
        for (i = 1; i <= count; i++) print reservoir[i] > "pilot.fastq"
        close("pilot.fastq")
        mean_len = (count > 0) ? int(sumlen / count + 0.5) : 0
        print format, total + 0, count + 0, mean_len
      }
    ')
    read format total_reads pilot_reads read_length <<< "\$stats"

    if [ "\$format" = "unknown" ]; then
        echo "ERROR: sample ${meta.id} is neither FASTQ nor FASTA (first character is not '@' or '>')." >&2
        exit 1
    fi
    if [ "\$total_reads" -eq 0 ]; then
        echo "ERROR: sample ${meta.id} contains zero reads. Check the input FASTQ is not empty or truncated." >&2
        exit 1
    fi

    skip_reason=""
    if [ "\$format" = "fasta" ]; then
        skip_reason="fasta_input"
    elif [ "\$total_reads" -le ${min_fragments} ]; then
        skip_reason="below_min_fragments"
    fi

    if [ -z "\$skip_reason" ]; then
        # Assumes meta.id contains no single quotes/shell metacharacters; it comes straight from
        # the input samplesheet column and is not sanitized upstream (see main.nf CSV parsing).
        sourmash sketch dna -p k=31,scaled=1000,abund --name '${meta.id}' pilot.fastq -o pilot.sig

        # --threshold-bp 0 is required: the default 50kbp threshold makes gather write no result
        # row for low-overlap samples, which is exactly the case this measurement exists for.
        # No `|| true` here: a real sourmash gather (this container's version) exits 0 both when
        # it finds a match and when it finds none (writing no CSV in the latter case), so a
        # nonzero exit is a genuine failure (OOM, corrupt sig, version change) that must not be
        # silently reported as 0.0 on-target.
        sourmash gather pilot.sig ${reference_sig} --threshold-bp 0 -o gather.csv
    fi

    TOTAL_READS=\$total_reads READ_LENGTH=\$read_length PILOT_READS=\$pilot_reads \\
    SKIP_REASON=\$skip_reason FORMAT=\$format OUT='${meta.id}.metrics.json' python3 -c '
import csv, json, os
rows = []
if os.path.exists("gather.csv"):
    with open("gather.csv") as fh:
        rows = list(csv.DictReader(fh))
skip_reason = os.environ["SKIP_REASON"] or None
if skip_reason:
    fraction = None
else:
    fraction = float(rows[0]["f_unique_weighted"]) if rows else 0.0
json.dump({"onTargetFraction": fraction,
           "skipReason": skip_reason,
           "format": os.environ["FORMAT"],
           "totalReads": int(os.environ["TOTAL_READS"]),
           "readLength": int(os.environ["READ_LENGTH"]),
           "pilotReads": int(os.environ["PILOT_READS"])},
          open(os.environ["OUT"], "w"), indent=2)
'
    """

    stub:
    // Wiring-test placeholders only; these are not meaningful measurements.
    """
    echo '{"onTargetFraction": 0.5, "format": "fastq", "totalReads": 1000000, "readLength": 150, "pilotReads": 100000}' > ${meta.id}.metrics.json
    """
}
