process EXPAND_SRX_IDS {
    tag "$sample_id:$accession"
    label 'process_single'

    container 'quay.io/biocontainers/entrez-direct:21.6--he881be0_0'

    input:
    tuple val(sample_id), val(accession), val(var1)

    output:
    path("expanded_rows.txt"), emit: rows

    script:
    def safe_var1 = var1 ?: ""
    """
    if echo "${accession}" | grep -qi '^SRX'; then
        esearch -db sra -query "${accession}[Accession]" \\
            | efetch -format runinfo \\
            | grep -v '^Run,' \\
            | grep -v '^\$' \\
            | cut -d',' -f1 \\
            | while read -r srr_id; do
                printf '%s,%s,,%s\\n' "${sample_id}" "\${srr_id}" "${safe_var1}"
              done > expanded_rows.txt
        if [ ! -s expanded_rows.txt ]; then
            echo "ERROR: No SRR IDs found for accession ${accession}" >&2
            exit 1
        fi
    else
        printf '%s,%s,,%s\\n' "${sample_id}" "${accession}" "${safe_var1}" > expanded_rows.txt
    fi
    """
}
