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
  --referenceFasta /path/to/target_organism.fasta \
  -resume -C <config>

# Prepare samples from local FASTQ files
nextflow run VEuPathDB/ngs-samples-nextflow -r main \
  --fromSra false \
  --input /path/to/samplesheet_dir \
  --outDir /path/to/output \
  --referenceFasta /path/to/target_organism.fasta \
  -resume -C <config>
```

`--referenceFasta` is required on every run and is validated before any process is
submitted.

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
| `--referenceFasta` | **Required.** Target organism FASTA. Genome size is measured from it and it is used to estimate each sample's on-target fraction. Gzipped FASTA is accepted |
| `--assayType` | `DNASeq`, `RNASeq`, or `ChipSeq` — determines the read-subsampling target (default `DNASeq`). Unrecognized values fail loudly |
| `--targetCoverage` | Coverage target for DNASeq/ChipSeq (default `60`) |
| `--minOnTargetFraction` | Fraction floor, which doubles as the inflation cap — never retain more than `1/minOnTargetFraction` times a clean sample's requirement (default `0.05`) |
| `--minPlausibleFraction` | Below this a sample is flagged in `sample_metrics.csv` (default `0.01`); every sample flagged usually means the wrong `--referenceFasta` |
| `--pilotSize` | Reads drawn per sample to estimate on-target fraction (default `100000`) |
| `--maxDownloadSize` | Maximum SRA run size `prefetch` will download (default `50G`); raise if prefetch skips a run for exceeding sra-tools' default 20G limit |

`--genomeSize` and `--maxReads` have been removed. Genome size is now measured from
`--referenceFasta`, and the raw-read target is always derived from the measured on-target
fraction rather than set manually. See "Contamination-aware subsampling" below.

## Outputs

The pipeline publishes to `--outDir`:

- `samplesheet.csv` — `sample,fastq_1,fastq_2,var1`. This contract is stable; downstream
  workflows can rely on the column set.
- `sample_metrics.csv` — per-sample measurements:
  `sample,on_target_fraction,total_reads,raw_reads_used,estimated_coverage,read_length,pilot_reads,flagged`.
  `estimated_coverage` is genome-relative and is left empty for RNASeq, where the depth
  target is a fixed read count. It can also legitimately read far below `--targetCoverage`
  when a sample simply doesn't contain enough reads to reach the target — `raw_reads_used`
  is capped at the reads actually available, so a low figure there can mean "small sample"
  rather than "sequencing failure".
- Subsampled FASTQ files — `{sample}_subsampled.fastq.gz` (single-end) or
  `{sample}_1_subsampled.fastq.gz` / `{sample}_2_subsampled.fastq.gz` (paired-end)

## Contamination-aware subsampling

Subsampling targets *on-target* reads, not raw reads. Each sample's on-target fraction is
estimated from a pilot draw of reads using k-mer containment against `--referenceFasta`, and
the number of raw reads retained is inflated by that fraction. A sample that is 15% target
therefore retains ~6.7x more raw reads than a clean one, and both reach the requested
coverage after alignment.

If `sample_metrics.csv` shows every sample flagged, `--referenceFasta` is almost certainly not
the organism the reads came from.
