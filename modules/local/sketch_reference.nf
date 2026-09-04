process SKETCH_REFERENCE {
    label 'process_low'

    container 'quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0'

    input:
    path fasta

    output:
    path("reference.sig"),  emit: sig
    path("reference.json"), emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    n_seqs=\$(grep -c '^>' ${fasta} || true)
    if [ "\$n_seqs" -lt 1 ]; then
        echo "ERROR: ${fasta} contains no FASTA headers. Expected a nucleotide FASTA." >&2
        exit 1
    fi

    genome_size=\$(grep -v '^>' ${fasta} | tr -cd 'ACGTNacgtn' | wc -c)
    if [ "\$genome_size" -lt 1000 ]; then
        echo "ERROR: ${fasta} has only \$genome_size nucleotide bases." >&2
        exit 1
    fi

    sourmash sketch dna -p k=31,scaled=1000 --name-from-first ${fasta} -o reference.sig

    GENOME_SIZE=\$genome_size N_SEQS=\$n_seqs python3 -c '
import json, os
json.dump({"genomeSize": int(os.environ["GENOME_SIZE"]),
           "sequences": int(os.environ["N_SEQS"])},
          open("reference.json", "w"), indent=2)
'
    """

    stub:
    """
    touch reference.sig
    echo '{"genomeSize": 200000, "sequences": 1}' > reference.json
    """
}
