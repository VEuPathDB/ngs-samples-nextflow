# Contamination-Aware Subsampling

**Date:** 2026-09-04
**Status:** Approved design, not yet implemented

## Problem

The pipeline subsamples reads before any alignment. `calculateMaxReads()` in `main.nf` is a
pure function of `assayType` and `genomeSize`, so the invariant it enforces is "N raw reads."
What consumers actually need is "N on-target reads."

Those are the same number only for a pure sample. In a mixed sample with substantial host
contamination, cutting to N raw reads delivers far fewer than N on-target reads, and every
downstream workflow gets less coverage than it asked for. The pipeline is asserting a
coverage guarantee it lacks the information to make, because it has never seen an alignment.

## Why the fix belongs here

Subsampling exists to save alignment compute in child workflows and to save storage on
retained FASTQ. Both require cutting before alignment and near acquisition. Post-alignment
downsampling cannot serve either goal, and moving the policy into each child workflow would
duplicate it across repositories while preserving the same wrong invariant.

The fix is not to relocate the cut. It is to measure the on-target fraction here, so the cut
is made against a number that means something.

## Approach

Insert a measurement step between concatenation and subsampling. The user passes a reference
FASTA; the pipeline sketches it once, estimates each sample's on-target fraction from a small
pilot draw of reads, and inflates the raw read target accordingly:

```
rawReadsNeeded = targetOnTargetReads / onTargetFraction
```

A sample that is 15% target gets ~6.7x more raw reads retained than a clean one, and both
land at the requested coverage.

### Why k-mer containment

The estimate needs only the on-target fraction, so contamination can be defined as
"everything that is not the target." That requires one small artifact per target organism —
no host genome, no host database, no aligner index.

`sourmash` FracMinHash sketches at `k=31,scaled=1000` are a few hundred KB for a 23Mb
genome. `sourmash gather` reports `f_unique_weighted`: the abundance-weighted fraction of the
pilot's k-mers explained by the reference. That value is the on-target fraction directly.

**`gather` must be run with `--threshold-bp 0`.** The default `--threshold-bp 50000` causes
gather to exit without writing any result row when overlap is small — precisely the
low-fraction regime this design exists to handle. Verified against synthetic mixtures with
sourmash 4.8.14 and a 2.7Mb *L. major* contig, 100k-read pilots:

| True target fraction | default threshold | `--threshold-bp 0` |
|---|---|---|
| 10% | 0.1076 | 0.1076 |
| 1% | 0.0106 | 0.0106 |
| 0.2% | *no CSV written* | 0.0018 |

Without the flag, a 0.2% sample is silently indistinguishable from a wrong reference FASTA,
which defeats the `minPlausibleFraction` diagnostic entirely.

Kraken2 with a two-genome database would give per-read classification and identify *what* the
contaminant is, at the cost of multi-GB resident memory per task. Given the target/host size
asymmetry and `maxForks = 2`, sourmash is the better trade. The measurement is isolated in
its own process so this choice is reversible.

### Known estimator biases

- Conserved regions and rRNA shared between target and host inflate the estimate. `k=31`
  keeps this modest but nonzero.
- A strain divergent from the reference assembly deflates the estimate, causing over-pull.
  This is the safer direction — it costs storage, not coverage.
- `scaled=1000` is lossy at very low abundance, so sub-1% estimates are least trustworthy.
  This is the same regime the fraction floor governs.

### Reference scope

Per-run, via `params.referenceFasta`. Studies never mix organisms. The reference is a value
channel, sketched once, broadcast to all samples.

Passing the FASTA rather than a prebuilt sketch keeps the estimator swappable and requires no
hosting infrastructure. It also allows genome size to be measured rather than asserted,
retiring `params.genomeSize`.

## Architecture

`main.nf` and `retrieve_from_sra.nf` currently run identical
`CONCATENATE → SUBSAMPLE → FORMAT → collectFile` chains. Adding processes to that chain means
maintaining two copies of a longer chain, so the tail is extracted first.

### New: `workflows/prepare_samples.nf`

```
PREPARE_SAMPLES(grouped_reads, reference_metrics)
  CONCATENATE_FASTQ
    → MEASURE_SAMPLE
    → DepthPolicy (pure Groovy, applied in a map operator)
    → SUBSAMPLE_FASTQ
    → FORMAT_INPUT_FROM_SRA
    → collectFile
```

Called by both the local branch of `main.nf` and by `RETRIEVE_FROM_SRA`. The `fromSra` fork
reduces to what it should be — how reads are obtained — and everything downstream of reads
existing is one code path.

### New: `modules/local/sketch_reference.nf`

`SKETCH_REFERENCE` runs once per pipeline. Takes `params.referenceFasta`; emits
`reference.sig` and `reference.json` carrying the measured genome size. Nextflow's cache makes
reruns against the same genome free.

Sanity-checks its own output: an implausibly small sketch or a zero genome size means a GFF,
an index, or an empty file was passed.

Container: `quay.io/biocontainers/sourmash`.

### New: `modules/local/measure_sample.nf`

`MEASURE_SAMPLE` takes `tuple(meta, reads)` plus the reference sketch; emits
`tuple(meta, path("${meta.id}.metrics.json"))`.

Performs the pilot draw (`seqtk sample`), sketches it with abundance tracking, runs
`sourmash gather`, and records the total read count, read length, and pilot size.

The total read count moves here from `SUBSAMPLE_FASTQ`, which currently pays for a full
`zcat | wc -l`. Net cost is unchanged, but every measured fact about a sample now originates
in one process, and `SUBSAMPLE_FASTQ` becomes pure execution with nothing to decide.

Metrics JSON fields: `onTargetFraction`, `totalReads`, `readLength`, `pilotReads`.

### New: `lib/DepthPolicy.groovy`

Plain Groovy class. No process, no container. Consumes the metrics JSON plus `assayType` and
the measured genome size; returns the target raw read count. Applied in a `.map{}` operator.

Kept separate from measurement because it is the piece most likely to change, it is testable
without Nextflow or FASTQ data, and the separation keeps a future two-stage
characterize-then-acquire split reachable as an added entry point rather than a rewrite.

Two steps, distinct because they fail differently:

```groovy
// 1. How many on-target reads do we want?
targetOnTarget = assayType == "RNASeq"
    ? 20_000_000
    : clamp(genomeSize * targetCoverage / readLength, 1_000_000, 100_000_000)

// 2. How many raw reads must we retain to get them?
fraction = max(metrics.onTargetFraction, policy.minOnTargetFraction)
rawReads = min(ceil(targetOnTarget / fraction), metrics.totalReads)
```

`readLength` and `genomeSize` are now measured values, replacing the hardcoded `150` and the
hand-entered `params.genomeSize`.

### Modified processes

**`SUBSAMPLE_FASTQ`** — input becomes `tuple(meta, reads, target_reads)` instead of
`val max_reads`. Drops read counting and receives the total from upstream. Paired detection
switches from `reads.size() == 2` to `meta.hasPairedReads`, which the meta already carries.

**`CONCATENATE_FASTQ`** — short-circuits to a symlink when `file_list.size() == 1`. The
single-run case currently decompresses and recompresses to produce a byte-equivalent file.
Independent of the rest of this design; included because concat is the pipeline's actual peak
disk high-water mark (~2x pooled size, since inputs stay staged alongside the merged output).

**`main.nf`** — `calculateMaxReads()` is deleted; its logic moves into `DepthPolicy`.

## Multi-run samples

Samples with multiple SRA accessions are concatenated into a single pooled library, then cut
once. The deliverable is a sample, not a run.

A single random draw from the pool means each run contributes proportional to its size. This
is correct here: these are pooled libraries, not balanced replicates. No equal-contribution
policy is needed.

Cutting each run to `target/N` is explicitly rejected — it assumes equal-sized runs and turns
one coverage target into N approximations. A coarse proportional pre-trim before concat would
reduce peak disk for samples with very many runs; omitted as unnecessary complexity until
that case is shown to exist.

## Guardrails

**One knob, not two.** `minOnTargetFraction` *is* the inflation cap — a floor of 0.05 means
never retaining more than 20x what a clean sample would need. A separate absolute
`maxRawReads` cap would govern the same behavior, and `metrics.totalReads` is already a hard
natural ceiling. Default: **0.05**.

**The degenerate case matters more than the floor.** A near-zero fraction has two causes: a
genuinely filthy sample, or the wrong reference FASTA. Those want opposite responses, and the
floor silently treats the second as the first — quietly burning 20x storage across the whole
run for no coverage benefit.

`minPlausibleFraction` (default 0.01): samples below it are flagged.

## Error handling

**Fail-fast on the reference.** `params.referenceFasta` is validated at workflow entry —
exists, readable, non-empty — before any `prefetch` runs. `SKETCH_REFERENCE` adds its own
output sanity check. Together these catch the entire "wrong file" class at minute one.

**Wrong organism, deliberately without a barrier.** Failing the run when every sample trips
`minPlausibleFraction` requires collecting all metrics before deciding. Placing that barrier
before `SUBSAMPLE_FASTQ` would force every download to finish before any subsampling starts,
so the wrong-genome error would surface only after the full download cost was already paid.
The barrier buys nothing at the moment it is most wanted.

Instead: warn per flagged sample as it is measured, carry a `flagged` column in the metrics
sidecar, and check the all-samples-flagged condition in `workflow.onComplete`, where it errors
with a message pointing at the reference FASTA. The first warning appears in the log as soon
as the first sample is measured, which is the practical early signal.

**Degenerate cases, handled rather than crashed:**

- `sourmash gather` with no matches emits an empty result, not an error. The parser returns
  `0.0`. Because gather runs with `--threshold-bp 0`, an empty result means genuine
  zero overlap rather than sub-threshold overlap, so `0.0` is trustworthy here.
- Samples smaller than the pilot size: `seqtk sample` returns everything. `pilotReads` in the
  metrics makes a noisy estimate visible rather than implicit.
- `totalReads` below the computed ask: no subsampling, symlink passthrough as today.

## Observability

The measurement is worthless downstream if it evaporates.

The samplesheet contract stays frozen at `sample,fastq_1,fastq_2,var1`. Metrics go to a
separate `sample_metrics.csv` in `params.outDir`, with columns: `sample`,
`on_target_fraction`, `total_reads`, `raw_reads_used`, `estimated_coverage`, `read_length`,
`pilot_reads`, `flagged`.

Assembled in a `.map{}` and written with `collectFile(keepHeader: true, storeDir:
params.outDir)`. No process and no container — every value is already in Groovy by that
point.

`estimated_coverage` is computed from reads actually delivered:
`rawReadsDelivered * onTargetFraction * readLength / genomeSize`. It is genome-relative and therefore not meaningful for RNASeq, where the depth target is a fixed read count rather than a coverage figure; the column is emitted as empty for that assay type.

## Parameters

**New:**

| Param | Default | Purpose |
|---|---|---|
| `referenceFasta` | *(required)* | Target organism FASTA |
| `minOnTargetFraction` | `0.05` | Fraction floor; caps inflation at 20x |
| `minPlausibleFraction` | `0.01` | Below this, flag the sample |
| `targetCoverage` | `60` | Coverage target for DNASeq/ChipSeq |
| `pilotSize` | `100000` | Pilot draw size |

**Retired:** `genomeSize` — now measured from the reference FASTA.

## Testing

**`lib/DepthPolicy.groovy` — unit tests, table-driven.** Clean sample; 15%-target sample;
fraction below floor; fraction below plausible threshold; `totalReads` under the ask; RNASeq
fixed path; both boundary clamps. This is arithmetic that silently produces plausible-looking
wrong numbers and it is the piece most likely to be edited later. As a pure Groovy class it
tests in milliseconds with no data.

**`MEASURE_SAMPLE` — nf-test against a known mixture.** A small synthetic reference plus a
FASTQ built from 1,000 target reads and 9,000 off-target reads; assert the estimate lands near
0.1. This is a correctness test of the estimator and the regression net if sourmash is ever
swapped for Kraken2 or real alignment.

**`PREPARE_SAMPLES` — nf-test on wiring**, using stubs.

**Stub blocks on all new processes** so `-stub-run` continues to exercise the full DAG.

## Rejected alternatives

**Move subsampling into child workflows.** Duplicates the depth policy across repositories
and preserves the wrong invariant in each of them. Also cannot serve the storage goal.

**Two-pass acquisition (pilot before download).** The only option that reduces *peak* disk
rather than retained disk. Rejected on two grounds: it abandons `SRATOOLS_PREFETCH` and its
`retry_with_backoff.sh` machinery, since prefetch is all-or-nothing on the `.sra`, trading a
solved reliability problem for a storage optimization; and dumping the first N spots makes the
delivered reads flowcell-ordered rather than randomly sampled, which is a quiet statistical
asterisk on data handed to every downstream workflow. `maxForks = 2` already bounds peak
transient disk to two full runs.

**Split into separate characterize and acquire entry points.** The right shape if depth policy
later becomes contested across teams. YAGNI today — two invocations and more operator surface
for a problem that does not yet exist. The measurement/policy separation keeps it reachable.

**Expected-fraction column in the samplesheet.** Zero new dependencies, but accuracy is only
as good as a human guess, and it would not have caught the problem that motivated this work.

**Full aligner index for the estimate.** Most accurate, but the acquisition layer would need
per-organism indexes, a new config surface, and cache management.
