# ngs-samples-nextflow

A Nextflow pipeline that prepares NGS FASTQ samples — from SRA or local files — into a standardized, subsampled samplesheet for downstream VEuPathDB analysis pipelines.

## Overview

Downstream VEuPathDB genomics pipelines (e.g. `dnaseq-nextflow`, `bulk-rnaseq-nextflow`) expect a consistent samplesheet of per-sample FASTQ files at a manageable read depth. This pipeline is the preparation step that produces that input: it either downloads raw reads from NCBI's Sequence Read Archive or takes existing local FASTQ files, concatenates multiple files belonging to the same sample (e.g. technical replicates or multiple SRR runs under one SRX), subsamples each sample down to a target read count, and emits a formatted samplesheet with absolute paths to the processed files.

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/) (default), [Singularity](https://sylabs.io/singularity/)/Apptainer, or LSF — see `conf/docker.config`, `conf/singularity.config`, `conf/lsf.config`

## Usage

The pipeline has a single entry point (the default `workflow`), which branches into SRA-download or local-file mode based on `--fromSra`.

```bash
# Download and prepare samples from SRA
nextflow run VEuPathDB/ngs-samples-nextflow -r main \
  --fromSra true \
  --input /path/to/samplesheet_dir \
  --outDir /path/to/output \
  -resume -C <config>

# Prepare samples from local FASTQ files
nextflow run VEuPathDB/ngs-samples-nextflow -r main \
  --fromSra false \
  --input /path/to/samplesheet_dir \
  --outDir /path/to/output \
  -resume -C <config>
```

### Input samplesheet

A CSV with a header row and columns `sample, fastq_1, fastq_2, var1`:

```csv
sample,fastq_1,fastq_2,var1
sample1,SRR123456,,control
sample1,SRX789012,,control
sample2,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,treatment
```

- In SRA mode, `fastq_1` holds an SRA accession — either a run (SRR) or an experiment (SRX) accession. SRX accessions are expanded to their constituent SRR runs via NCBI Entrez Direct before download.
- In local mode, `fastq_1`/`fastq_2` are paths (relative to `--input`) to existing FASTQ files; `fastq_2` is left empty for single-end data.
- Multiple rows sharing the same `sample` ID are grouped and their reads concatenated into one file (or one R1/R2 pair) per sample.
- `var1` is passed through as arbitrary per-sample metadata into the output samplesheet.

## Key parameters

| Parameter | Description |
|---|---|
| `--input` | Directory containing the input samplesheet |
| `--samplesheetName` | Samplesheet filename within `--input` (default `samplesheet.csv`) |
| `--fromSra` | `true` to download reads from SRA, `false` to use local FASTQ files (default `true`) |
| `--outDir` | Output directory for the processed FASTQs and final samplesheet |
| `--assayType` | `DNASeq`, `RNASeq`, or `ChipSeq` — determines the read-subsampling target (default `DNASeq`) |
| `--genomeSize` | Genome size in base pairs, used to calculate the DNASeq/ChipSeq read cap (default `3000000000`, human) |
| `--maxDownloadSize` | Maximum SRA run size `prefetch` will download (default `50G`); raise if prefetch skips a run for exceeding sra-tools' default 20G limit |

Read subsampling is calculated automatically: RNASeq samples are capped at a fixed 20,000,000 reads, while DNASeq/ChipSeq samples are capped at `genomeSize * 60x coverage / 150bp read length`, bounded between 1,000,000 and 100,000,000 reads. Subsampling uses `seqtk sample` with a fixed seed for reproducibility, preserving read pairing for paired-end data.

## Output

The pipeline publishes to `--outDir`:

- Per-sample subsampled FASTQ file(s) — `{sample}_subsampled.fastq.gz` (single-end) or `{sample}_1_subsampled.fastq.gz` / `{sample}_2_subsampled.fastq.gz` (paired-end)
- A formatted samplesheet (named per `--samplesheetName`) with absolute paths to the processed FASTQ files, ready for downstream pipelines
