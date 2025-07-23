# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline for NGS (Next Generation Sequencing) sample processing. It handles two main use cases:
1. **SRA Download**: Downloads FASTQ files from NCBI SRA and creates formatted samplesheets
2. **Local Files**: Processes existing local FASTQ files and creates absolute path samplesheets

## Common Commands

### Running the Pipeline

```bash
# For SRA download (fromSra=true)
nextflow run main.nf --fromSra true --input /path/to/samplesheet --outDir /path/to/output

# For local files (fromSra=false, default)
nextflow run main.nf --input /path/to/samplesheet --outDir /path/to/output

# With custom samplesheet name
nextflow run main.nf --samplesheetName custom.csv --input /path/to/data

# Specify assay type and genome size for subsampling
nextflow run main.nf --assayType RNASeq --genomeSize 120000000 --input /path/to/data

# Override automatic read calculation with manual limit
nextflow run main.nf --maxReads 10000000 --input /path/to/data
```

### Testing

```bash
# Run nf-core module tests
nextflow test modules/nf-core/sratools/fasterqdump/tests/main.nf.test
nextflow test modules/nf-core/sratools/prefetch/tests/main.nf.test
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

- **main.nf**: Entry point that orchestrates the workflow based on `fromSra` parameter
- **workflows/retrieve_from_sra.nf**: Workflow for downloading from SRA using prefetch → fasterqdump → format
- **modules/local/format_input_from_sra.nf**: Custom process to create properly formatted samplesheets
- **modules/nf-core/**: Standard nf-core modules for SRA tools (prefetch, fasterqdump)

### Key Parameters

- `input`: Directory containing input samplesheet (default: `$launchDir/data/samplesheet`)
- `samplesheetName`: Name of samplesheet file (default: `samplesheet.csv`)
- `fromSra`: Boolean to determine SRA download vs local files (default: `false`)
- `outDir`: Output directory (default: `$launchDir/ngs-samples-output`)
- `workDir`: Nextflow work directory (default: `$launchDir/ngs-samples-work`)

#### Subsampling Parameters
- `assayType`: Type of sequencing assay - "DNASeq" or "RNASeq" (default: `"DNASeq"`)
- `genomeSize`: Genome size in base pairs for coverage calculation (default: `"3000000000"` for human)
- `maxReads`: Manual override for maximum reads per sample (default: `null` - uses automatic calculation)

### Input Samplesheet Format

Expected CSV format with header:
- **Column 0**: Sample ID (can have multiple rows with same ID for concatenation)
- **Column 1**: FASTQ1 path (or SRA ID when fromSra=true)  
- **Column 2**: FASTQ2 path (optional, for paired-end)
- **Column 3**: Additional variable (var1)

**Multi-file concatenation**: If multiple rows have the same sample ID, their FASTQ files will be automatically concatenated:
- For single-end data: all files are concatenated into `{sample_id}.fastq.gz`
- For paired-end data: R1 files are concatenated into `{sample_id}_1.fastq.gz` and R2 files into `{sample_id}_2.fastq.gz`

### Process Flow

1. **SRA Mode**: `samples` → group by ID → `SRATOOLS_PREFETCH` → `SRATOOLS_FASTERQDUMP` → `CONCATENATE_FASTQ` → `SUBSAMPLE_FASTQ` → `FORMAT_INPUT_FROM_SRA` → output samplesheet
2. **Local Mode**: `samples` → group by ID → `CONCATENATE_FASTQ` → `SUBSAMPLE_FASTQ` → `FORMAT_INPUT_FROM_SRA` → output samplesheet

Both modes include automatic file concatenation and intelligent subsampling based on assay type and genome size.

### Container Management

- Uses Docker by default (`docker.enabled = true`)
- All images pulled from `quay.io` registry
- Custom Alpine bash container for formatting: `docker.io/veupathdb/alpine_bash:1.0.0`
- SRA tools use biocontainers images

### Error Handling

- `SRATOOLS_FASTERQDUMP` has retry logic: falls back from `fasterq-dump` to `fastq-dump` on failure
- `SRATOOLS_PREFETCH` uses retry template with exponential backoff
- Maximum 2 concurrent processes (`maxForks = 2`)

## File Locations

- Configuration files: `conf/` directory
- Local modules: `modules/local/` (includes `concatenate_fastq.nf` and `format_input_from_sra.nf`)
- nf-core modules: `modules/nf-core/`
- Main workflow: `workflows/`
- Test files: Located in each module's `tests/` subdirectory

## New Features

### FASTQ File Concatenation
- **Purpose**: Handles samplesheets with multiple rows per sample ID
- **Implementation**: Uses `CONCATENATE_FASTQ` process to merge files before formatting
- **Supported patterns**: Automatically detects R1/R2 files using `_1.fastq`, `_2.fastq`, `_R1`, `_R2` patterns
- **Output**: Single concatenated file per sample (or paired files for paired-end data)

### Read Subsampling
- **Purpose**: Limits the number of reads per sample to optimize downstream processing
- **Implementation**: Uses `SUBSAMPLE_FASTQ` process with seqtk for random subsampling
- **Coverage calculation**: 
  - DNASeq: 30x coverage target
  - RNASeq: 50x coverage target (higher due to expression variation)
- **Read limits**: Bounded between 1M and 100M reads per sample
- **Paired-end handling**: Maintains read pairing using consistent random seed
- **Container**: Uses `staphb/seqtk:1.4` Docker image