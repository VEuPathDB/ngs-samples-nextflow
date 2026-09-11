# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline for NGS (Next Generation Sequencing) sample processing. It handles two main use cases:
1. **SRA Download**: Downloads FASTQ files from NCBI SRA and creates formatted samplesheets
2. **Local Files**: Processes existing local FASTQ files and creates absolute path samplesheets

Both use cases converge on a contamination-aware subsampling step: each sample's on-target
fraction is estimated against `--referenceFasta`, and subsampling targets on-target reads
rather than raw reads.

## Common Commands

### Running the Pipeline

```bash
# For SRA download (fromSra=true)
nextflow run main.nf --fromSra true --input /path/to/samplesheet --outDir /path/to/output \
  --referenceFasta /path/to/target_organism.fasta

# For local files (fromSra=false, default)
nextflow run main.nf --input /path/to/samplesheet --outDir /path/to/output \
  --referenceFasta /path/to/target_organism.fasta

# With custom samplesheet name
nextflow run main.nf --samplesheetName custom.csv --input /path/to/data \
  --referenceFasta /path/to/target_organism.fasta

# Specify assay type for subsampling (targetCoverage applies to DNASeq only; RNASeq and
# ChipSeq use fragment-count targets instead)
nextflow run main.nf --assayType RNASeq --input /path/to/data \
  --referenceFasta /path/to/target_organism.fasta
```

`--referenceFasta` is required on every run — the pipeline fails fast, before submitting any
process, if it is missing.

### Testing

```bash
# Requires nf-test on PATH (install: curl -fsSL https://code.askimed.com/install/nf-test | bash)
export PATH="$HOME/bin:$PATH"

# Run a specific test file
nf-test test modules/local/depth_policy_tests/main.nf.test

# Run by tag
nf-test test --tag depth_policy

# NOTE: the vendored modules/nf-core/sratools/** tests require
# params.modules_testdata_base_path and do not currently run in this repo.
```

### Development with Different Executors

```bash
# Using Docker (default configuration)
nextflow run main.nf -c conf/docker.config

# Using Singularity 
nextflow run main.nf -c conf/singularity.config

# Using LSF scheduler
nextflow run main.nf -c conf/lsf.config
```

## Architecture

### Pipeline Structure

- **main.nf**: Entry point. Validates `--referenceFasta`, runs `SKETCH_REFERENCE` once, then
  orchestrates SRA vs local mode based on the `fromSra` parameter
- **workflows/retrieve_from_sra.nf**: Workflow for downloading from SRA using prefetch → fasterqdump, then handing off to `PREPARE_SAMPLES`
- **workflows/prepare_samples.nf**: Shared subworkflow (both modes converge here) — concatenation, on-target measurement, depth-policy calculation, subsampling, and samplesheet formatting
- **modules/local/format_input_from_sra.nf**: Custom process to create properly formatted samplesheets
- **modules/local/sketch_reference.nf**: Validates and sketches `--referenceFasta`, measuring genome size
- **modules/local/measure_sample.nf**: Estimates each sample's on-target fraction via a pilot draw and k-mer containment
- **modules/local/depth_policy.nf**: Pure functions computing the raw-read target from assay type, genome size, and measured on-target fraction
- **modules/nf-core/**: Standard nf-core modules for SRA tools (prefetch, fasterqdump)

### Key Parameters

- `input`: Directory containing input samplesheet (default: `$launchDir/data/`)
- `samplesheetName`: Name of samplesheet file (default: `samplesheet.csv`)
- `fromSra`: Boolean to determine SRA download vs local files (default: `true`)
- `outDir`: Output directory (default: `$launchDir/ngs-samples-output`)
- `workDir`: Nextflow work directory (default: `$launchDir/ngs-samples-work`)
- `maxDownloadSize`: Maximum SRA run size `prefetch` will download (default: `"50G"`)

#### Subsampling Parameters
- `referenceFasta`: **Required.** Target organism FASTA. Used to estimate each sample's
  on-target fraction so subsampling targets on-target reads rather than raw reads. Genome
  size is measured from this file. Gzipped FASTA is accepted.
- `assayType`: "DNASeq", "RNASeq", or "ChipSeq" (default: `"DNASeq"`). Unrecognized values
  fail loudly rather than silently defaulting.
- `targetCoverage`: Coverage target for DNASeq only (default: `60`). RNASeq and ChipSeq set
  their targets as fragment counts instead
- `minOnTargetFraction`: Fraction floor, which doubles as the inflation cap (default: `0.05`,
  i.e. never retain more than 20x a clean sample's requirement)
- `minPlausibleFraction`: Below this a sample is flagged (default: `0.01`). All samples
  flagged usually means the wrong `referenceFasta`.
- `pilotSize`: Reads drawn per sample to estimate contamination (default: `100000`)

`genomeSize` has been removed — genome size is now measured from `referenceFasta`.

### Input Samplesheet Format

Expected CSV format with header:
- **Column 0**: Sample ID (can have multiple rows with same ID for concatenation)
- **Column 1**: FASTQ1 path (or SRA ID when fromSra=true)  
- **Column 2**: FASTQ2 path (optional, for paired-end)
- **Column 3**: Additional variable (var1)

**Multi-file concatenation**: If multiple rows have the same sample ID, their FASTQ files will be automatically concatenated:
- For single-end data: all files are concatenated into `{sample_id}_concat.fastq.gz`
- For paired-end data: R1 files are concatenated into `{sample_id}_concat_1.fastq.gz` and R2 files into `{sample_id}_concat_2.fastq.gz`

### Process Flow

`SKETCH_REFERENCE` runs once per pipeline, then both modes converge on `PREPARE_SAMPLES`:

1. **SRA Mode**: `samples` → `EXPAND_SRX_IDS` → group → `SRATOOLS_PREFETCH` →
   `SRATOOLS_FASTERQDUMP` → `PREPARE_SAMPLES`
2. **Local Mode**: `samples` → group → `PREPARE_SAMPLES`

`PREPARE_SAMPLES` = `CONCATENATE_FASTQ` → `MEASURE_SAMPLE` → depth policy →
`SUBSAMPLE_FASTQ` → `FORMAT_INPUT_FROM_SRA`.

### Container Management

- Uses Docker by default (`docker.enabled = true`)
- All images pulled from `quay.io` registry
- Custom Alpine bash container for formatting: `docker.io/veupathdb/alpine_bash:1.0.0`
- SRA tools use biocontainers images
- `SKETCH_REFERENCE`/`MEASURE_SAMPLE` use `quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0`
- `SUBSAMPLE_FASTQ` uses `staphb/seqtk:1.4`

### Error Handling

- `SRATOOLS_FASTERQDUMP` has retry logic: falls back from `fasterq-dump` to `fastq-dump` on failure
- `SRATOOLS_PREFETCH` uses retry template with exponential backoff
- Maximum 2 concurrent processes (`maxForks = 2`)
- `SKETCH_REFERENCE` fails loudly (rather than producing a bogus genome size) if
  `--referenceFasta` has no FASTA headers, no sequence characters, or is not predominantly
  nucleotide (e.g. a protein FASTA)
- A sample below `minPlausibleFraction` is warned about individually and marked in
  `sample_metrics.csv`. This is not treated as a run-level error: low on-target fractions are
  expected wherever the target organism is sequenced out of host tissue

## File Locations

- Configuration files: `conf/` directory
- Local modules: `modules/local/` (includes `concatenate_fastq.nf`, `measure_sample.nf`, `depth_policy.nf`, `sketch_reference.nf`, `subsample_fastq.nf`, and `format_input_from_sra.nf`)
- nf-core modules: `modules/nf-core/`
- Workflows: `workflows/`
- Test files: Located in each module's `tests/` (or `<module>_tests/`) subdirectory

## New Features

### FASTQ File Concatenation
- **Purpose**: Handles samplesheets with multiple rows per sample ID
- **Implementation**: Uses `CONCATENATE_FASTQ` process to merge files before formatting
- **Supported patterns**: Automatically detects R1/R2 files using `_1.fastq`, `_2.fastq`, `_R1`, `_R2` patterns
- **Output**: Single concatenated file per sample (or paired files for paired-end data)

### Contamination-Aware Subsampling
- **Purpose**: Targets on-target reads rather than raw reads, so a contaminated sample still
  reaches the requested coverage after alignment instead of being under-sampled
- **Implementation**: `MEASURE_SAMPLE` draws a pilot of reads and estimates on-target
  fraction via sourmash k-mer containment against `SKETCH_REFERENCE`'s sketch of
  `--referenceFasta`; `modules/local/depth_policy.nf` turns that fraction into a raw-read
  target; `SUBSAMPLE_FASTQ` retains that many raw reads with seqtk
- **Depth targets**: three rules, one per assay. DNASeq targets `targetCoverage` (default 60x),
  linear in genome size. RNASeq uses a flat 20,000,000-fragment target, since transcriptome
  complexity barely tracks genome size. ChipSeq scales as the square root of genome size above
  a 3,000,000-fragment floor (which binds only below ~24Mb), matching the scaling implied by the
  modENCODE worm/fly and ENCODE human guidelines
- **Fragment limits**: Bounded between 1M and 100M fragments per sample
- **Paired-end handling**: Maintains read pairing using consistent random seed
- **Container**: Uses `staphb/seqtk:1.4` Docker image
