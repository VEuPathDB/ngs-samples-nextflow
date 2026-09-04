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

    # POSIX case, not bash [[ ]], so the container's /bin/sh can run this unmodified.
    case "${fasta}" in
        *.gz) reader="zcat" ;;
        *)    reader="cat" ;;
    esac

    n_seqs=\$(\$reader ${fasta} | grep -c '^>' || true)
    if [ "\$n_seqs" -lt 1 ]; then
        echo "ERROR: ${fasta} contains no FASTA headers. Expected a nucleotide FASTA." >&2
        exit 1
    fi

    # Sequence lines go to a temp file, not a shell variable, so a real multi-GB genome
    # isn't held in memory as one bash string. grep -v matching nothing (e.g. a
    # headers-only file) would otherwise exit 1 and, under set -e, abort with an
    # unexplained pipeline failure instead of the diagnostics below.
    seq_file=\$(mktemp)
    \$reader ${fasta} | grep -v '^>' > "\$seq_file" || true

    # Composition check (Fix 1): strict ACGTN only, deliberately excluding IUPAC ambiguity
    # codes. Those codes (R,Y,S,W,K,M,B,D,H,V) are ALSO valid amino acid letters, so adding
    # them here would let a protein FASTA score ~70% and defeat this check. Genome SIZE
    # below uses a wider set - that's a different, non-composition, count.
    total_chars=\$(tr -cd 'A-Za-z' < "\$seq_file" | wc -c)
    if [ "\$total_chars" -lt 1 ]; then
        echo "ERROR: ${fasta} has no sequence characters after the FASTA headers." >&2
        exit 1
    fi

    strict_nt_chars=\$(tr -cd 'ACGTNacgtn' < "\$seq_file" | wc -c)
    pct=\$(( strict_nt_chars * 100 / total_chars ))
    if [ "\$pct" -lt 90 ]; then
        echo "ERROR: ${fasta} is only \${pct}% strict nucleotide characters (ACGTN)." >&2
        echo "Is this a protein FASTA? SKETCH_REFERENCE expects a nucleotide genome FASTA." >&2
        exit 1
    fi

    # Genome size includes IUPAC ambiguity codes (real, if rare, in genome assemblies) -
    # a deliberately wider set than the strict composition check above.
    genome_size=\$(tr -cd 'ACGTNRYSWKMBDHVacgtnryswkmbdhv' < "\$seq_file" | wc -c)
    if [ "\$genome_size" -lt 1000 ]; then
        echo "ERROR: ${fasta} has only \$genome_size nucleotide bases." >&2
        exit 1
    fi
    rm -f "\$seq_file"

    # k=31 is sourmash's standard DNA k-mer size, tolerant of strain-level divergence;
    # scaled=1000 subsamples k-mers to keep the sketch small. sourmash reads gzip natively.
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
    # Arbitrary placeholder for wiring tests only - not derived from data. Do not couple
    # any test's expected genomeSize to this number.
    echo '{"genomeSize": 200000, "sequences": 1}' > reference.json
    """
}
