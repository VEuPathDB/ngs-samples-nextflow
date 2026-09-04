# Contamination-Aware Subsampling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make subsampling target "N on-target reads" instead of "N raw reads", by measuring each sample's on-target fraction with k-mer containment against a user-supplied reference FASTA and inflating the retained read count accordingly.

**Architecture:** A reference FASTA is sketched once per run. Each sample's concatenated FASTQ is systematically sub-sampled to a 100k-read pilot, sketched, and compared to the reference with `sourmash gather --threshold-bp 0`; `f_unique_weighted` is the on-target fraction. A pure Nextflow function turns that fraction plus assay type into a raw read count, which drives the existing seqtk subsampling. The duplicated `CONCATENATE → SUBSAMPLE → FORMAT` tail in `main.nf` and `retrieve_from_sra.nf` is extracted into a shared `PREPARE_SAMPLES` subworkflow first, so the new steps are added in one place.

**Tech Stack:** Nextflow DSL2 (25.10.x), nf-test, Docker, `quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0`, `staphb/seqtk:1.4`, `docker.io/veupathdb/alpine_bash:1.0.0`.

**Spec:** `docs/superpowers/specs/2026-09-04-contamination-aware-subsampling-design.md`

**Branch:** `contamination-aware-subsampling`

---

## Deviations from the spec (read before starting)

1. **`lib/DepthPolicy.groovy` becomes `modules/local/depth_policy.nf`.** The spec called for a
   Groovy class in `lib/`. Testing a `lib/` class requires a Gradle/JUnit setup this repo does
   not have. Nextflow functions defined in a `.nf` file are directly testable with nf-test's
   `nextflow_function` block, which is the toolchain the repo already vendors. The function
   stays pure — no `params` access, all config passed in — so every property the spec wanted
   is preserved.

2. **The pilot draw uses systematic sampling (`awk`, every Nth read) rather than
   `seqtk sample`.** This keeps `MEASURE_SAMPLE` in a single container. Systematic sampling
   spans the whole file, so it is unbiased for composition, which is the only property the
   estimate needs. `seqtk` is still used by `SUBSAMPLE_FASTQ` for the actual cut.

3. **Task 9 (paired-end mate accounting) is flagged optional.** It fixes a pre-existing
   inconsistency not covered by the spec. Read its preamble and decide before implementing.

4. **`PREPARE_SAMPLES` is verified by a real `-stub-run` of the pipeline rather than an
   nf-test workflow test** (Task 8 Step 5, Task 11 Step 2). The spec asked for an nf-test on
   the wiring; a stub run exercises the same channel topology plus the `main.nf` entry logic
   and param validation, which the workflow test would not cover.

5. **The symlink short-circuit covers single-run paired-end as well as single-end.** The spec
   mentioned only `file_list.size() == 1`. Single-run paired-end is the common SRA case and
   is where most of the saving is.

---

## File Structure

**Create:**

| Path | Responsibility |
|---|---|
| `nf-test.config` | nf-test root config |
| `tests/nextflow.config` | Config applied to all nf-tests (docker profile) |
| `modules/local/depth_policy.nf` | Pure functions: metrics + config → depth plan |
| `modules/local/depth_policy_tests/main.nf.test` | Unit tests for the policy |
| `modules/local/sketch_reference.nf` | `SKETCH_REFERENCE` process |
| `modules/local/sketch_reference_tests/main.nf.test` | Test for reference sketching |
| `modules/local/measure_sample.nf` | `MEASURE_SAMPLE` process |
| `modules/local/measure_sample_tests/main.nf.test` | Known-mixture estimator test |
| `workflows/prepare_samples.nf` | Shared post-acquisition subworkflow |
| `tests/fixtures/make_fixtures.py` | Generates synthetic reference + mixture FASTQs |
| `tests/fixtures/ref.fasta` | Synthetic 200kb reference (generated, committed) |
| `tests/fixtures/mix10.fastq.gz` | 10% target mixture (generated, committed) |
| `tests/fixtures/pair_1.fastq.gz`, `pair_2.fastq.gz` | Paired-end passthrough fixtures |
| `tests/fixtures/not_a_fasta.gff` | Negative fixture for reference validation |

**Modify:**

| Path | Change |
|---|---|
| `modules/local/concatenate_fastq.nf` | Symlink short-circuit for single-file samples |
| `modules/local/subsample_fastq.nf` | Take per-sample `target_reads` and `total_reads`; drop counting |
| `workflows/retrieve_from_sra.nf` | Delegate tail to `PREPARE_SAMPLES` |
| `main.nf` | Delete `calculateMaxReads()`; wire `SKETCH_REFERENCE` + `PREPARE_SAMPLES`; validate params; `onComplete` check |
| `nextflow.config` | New params; retire `genomeSize` |
| `CLAUDE.md` | Correct the bogus `nextflow test` commands; document new params |
| `README.md` | Document reference FASTA requirement and metrics sidecar |

---

## Task 1: Bootstrap nf-test

The repo vendors nf-core `.nf.test` files but has no nf-test binary and no `nf-test.config`,
so they have never run. `nextflow test` (as written in CLAUDE.md) is not a real Nextflow
subcommand. Everything downstream depends on a working test runner.

Note: the vendored `modules/nf-core/sratools/**` tests reference
`params.modules_testdata_base_path`, which is not set in this repo. They will fail. Do not fix
them — out of scope. Always run new tests by explicit path.

**Files:**
- Create: `nf-test.config`
- Create: `tests/nextflow.config`
- Modify: `.gitignore`

- [ ] **Step 1: Install the nf-test binary**

```bash
mkdir -p ~/bin
cd ~/bin && curl -fsSL https://code.askimed.com/install/nf-test | bash
export PATH="$HOME/bin:$PATH"
nf-test version
```

Expected: prints an nf-test version (2.x) and a Nextflow version.

- [ ] **Step 2: Create the nf-test root config**

Create `nf-test.config`:

```groovy
config {
    testsDir "."
    workDir ".nf-test"
    configFile "tests/nextflow.config"
}
```

**Do not add `profile "docker"`.** Verified during implementation: this repo has no
`profiles {}` block, so Nextflow fails with `Unknown configuration profile: 'docker'`.
Docker is already enabled unconditionally via `conf/docker.config`, and again in
`tests/nextflow.config`.

- [ ] **Step 3: Create the shared test config**

Create `tests/nextflow.config`:

```groovy
params {
    outDir = "$launchDir/test-output"
}

docker.enabled = true
docker.registry = 'quay.io'
```

- [ ] **Step 4: Ignore nf-test working dirs**

Append to `.gitignore` (create the file if absent):

```
.nf-test/
.nf-test.log
test-output/
```

- [ ] **Step 5: Verify the runner starts**

Run: `nf-test test --dryRun modules/nf-core/sratools/prefetch/tests/main.nf.test`
Expected: exit 0 with `SUCCESS: Executed 3 tests`. This proves the config parses and the
runner starts. (Do not use `--help`: nf-test 0.9.5 exits 2 on it.)

- [ ] **Step 6: Commit**

```bash
git add nf-test.config tests/nextflow.config .gitignore
git commit -m "test: bootstrap nf-test runner and config"
```

---

## Task 2: Depth policy (pure function)

This is the arithmetic that silently produces plausible-looking wrong numbers. It is tested
first and most heavily. It reads no `params` and touches no files.

**Contract:**

`depthPlan(Map metrics, Map policy) -> Map`

`metrics`: `onTargetFraction` (double), `totalReads` (long), `readLength` (int), `pilotReads` (int)

`policy`: `assayType` (String), `genomeSize` (long), `targetCoverage` (int),
`minOnTargetFraction` (double), `minPlausibleFraction` (double)

Returns: `targetOnTarget` (long), `effectiveFraction` (double), `rawReads` (long),
`flagged` (boolean), `estimatedCoverage` (Double or null)

**Files:**
- Create: `modules/local/depth_policy.nf`
- Test: `modules/local/depth_policy_tests/main.nf.test`

- [ ] **Step 1: Write the failing tests**

Create `modules/local/depth_policy_tests/main.nf.test`:

```groovy
nextflow_function {

    name "Test depth policy"
    script "../depth_policy.nf"
    tag "depth_policy"

    test("clean sample retains close to the on-target target") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.95d, totalReads: 200000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.targetOnTarget == 13200000L
            assert function.result.rawReads == 13894737L
            assert function.result.flagged == false
        }
    }

    test("15 percent target inflates the raw ask by the fraction") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.15d, totalReads: 200000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.rawReads == 88000000L
            assert function.result.flagged == false
        }
    }

    test("fraction below the floor is clamped to the floor") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.02d, totalReads: 500000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.effectiveFraction == 0.05d
            assert function.result.rawReads == 264000000L
        }
    }

    test("raw ask never exceeds the reads actually available") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.02d, totalReads: 100000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.rawReads == 100000000L
        }
    }

    test("fraction below the plausible threshold is flagged") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.005d, totalReads: 100000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.flagged == true
        }
    }

    test("zero fraction is flagged and does not divide by zero") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.0d, totalReads: 100000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.flagged == true
            assert function.result.rawReads == 100000000L
            assert function.result.estimatedCoverage == 0.0d
        }
    }

    test("RNASeq uses a fixed read target and reports no coverage") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 0.5d, totalReads: 200000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'RNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.targetOnTarget == 20000000L
            assert function.result.rawReads == 40000000L
            assert function.result.estimatedCoverage == null
        }
    }

    test("tiny genome clamps up to the one million read floor") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 1.0d, totalReads: 5000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 1000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.targetOnTarget == 1000000L
            assert function.result.rawReads == 1000000L
        }
    }

    test("huge genome clamps down to the hundred million read ceiling") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 1.0d, totalReads: 500000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 1000000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d]
                """
            }
        }
        then {
            assert function.success
            assert function.result.targetOnTarget == 100000000L
            assert function.result.rawReads == 100000000L
        }
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
export PATH="$HOME/bin:$PATH"   # nf-test lives at ~/bin/nf-test; not on non-interactive PATH
nf-test test modules/local/depth_policy_tests/main.nf.test
```

Expected: FAIL — the script `../depth_policy.nf` does not exist.

- [ ] **Step 3: Write the implementation**

Create `modules/local/depth_policy.nf`:

```groovy
/*
 * Pure depth-policy functions. No params access, no file IO — everything is passed in
 * so this can be unit tested without a pipeline run.
 */

def targetOnTargetReads(Map policy) {
    if (policy.assayType == "RNASeq") {
        return 20000000L
    }
    long raw = (long) ((policy.genomeSize * policy.targetCoverage) / policy.readLength)
    return Math.max(1000000L, Math.min(100000000L, raw))
}

def depthPlan(Map metrics, Map policy) {
    def resolved = policy + [readLength: metrics.readLength]
    long targetOnTarget = targetOnTargetReads(resolved)

    double observed = metrics.onTargetFraction as double
    double effective = Math.max(observed, policy.minOnTargetFraction as double)

    long rawReads = Math.min(
        (long) Math.ceil(targetOnTarget / effective),
        metrics.totalReads as long
    )

    Double estimatedCoverage = policy.assayType == "RNASeq"
        ? null
        : (rawReads * observed * metrics.readLength) / (policy.genomeSize as double)

    return [
        targetOnTarget   : targetOnTarget,
        effectiveFraction: effective,
        rawReads         : rawReads,
        flagged          : observed < (policy.minPlausibleFraction as double),
        estimatedCoverage: estimatedCoverage
    ]
}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
nf-test test modules/local/depth_policy_tests/main.nf.test
```

Expected: PASS, 9 tests.

- [ ] **Step 5: Commit**

```bash
git add modules/local/depth_policy.nf modules/local/depth_policy_tests/main.nf.test
git commit -m "feat: add pure depth policy function with unit tests"
```

---

## Task 3: Test fixtures

`MEASURE_SAMPLE` needs a reference and a mixture with a *known* target fraction, otherwise its
test proves nothing. The generator is committed alongside the fixtures so they can be
regenerated and so the expected fraction is auditable.

**Files:**
- Create: `tests/fixtures/make_fixtures.py`
- Create: `tests/fixtures/ref.fasta` (generated)
- Create: `tests/fixtures/mix10.fastq.gz` (generated)

- [ ] **Step 1: Write the generator**

Create `tests/fixtures/make_fixtures.py`:

```python
#!/usr/bin/env python3
"""Generate a synthetic reference and a mixture FASTQ with a known target fraction.

The 'host' reads are random sequence, which shares no k-mers with the reference. That makes
the expected on-target fraction exactly the mixing ratio, so the estimator test has a real
answer to check against.
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
REF_LEN = 200000
READ_LEN = 150
N_TARGET = 1000
N_HOST = 9000

random.seed(1729)


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


ref = "".join(random.choice("ACGT") for _ in range(REF_LEN))

with open(os.path.join(HERE, "ref.fasta"), "w") as out:
    out.write(">synthetic_target length=%d\n" % REF_LEN)
    for i in range(0, REF_LEN, 60):
        out.write(ref[i:i + 60] + "\n")

reads = []
for i in range(N_TARGET):
    pos = random.randint(0, REF_LEN - READ_LEN)
    seq = ref[pos:pos + READ_LEN]
    if random.random() < 0.5:
        seq = revcomp(seq)
    reads.append(("target_%d" % i, seq))

for i in range(N_HOST):
    reads.append(("host_%d" % i, "".join(random.choice("ACGT") for _ in range(READ_LEN))))

random.shuffle(reads)

# gzip.open embeds the current mtime and filename in its header, which would make the
# output non-deterministic across runs. Pin both so regenerating produces identical bytes.
with open(os.path.join(HERE, "mix10.fastq.gz"), "wb") as raw:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
        for name, seq in reads:
            out.write(("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq))).encode())

print("ref.fasta: %d bp" % REF_LEN)
print("mix10.fastq.gz: %d target / %d host = %.2f expected fraction"
      % (N_TARGET, N_HOST, N_TARGET / float(N_TARGET + N_HOST)))
```

- [ ] **Step 2: Generate the fixtures**

```bash
python3 tests/fixtures/make_fixtures.py
ls -la tests/fixtures/
```

Expected: `ref.fasta` ~203KB, `mix10.fastq.gz` ~500KB, and the printed line
`1000 target / 9000 host = 0.10 expected fraction`.

- [ ] **Step 3: Sanity check the fixture against sourmash by hand**

```bash
docker run --rm -v "$PWD/tests/fixtures":/d -w /d \
  quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0 sh -c '
    sourmash sketch dna -p k=31,scaled=1000 --name-from-first ref.fasta -o ref.sig &&
    sourmash sketch dna -p k=31,scaled=1000,abund --name mix mix10.fastq.gz -o mix.sig &&
    sourmash gather mix.sig ref.sig --threshold-bp 0 -o gather.csv &&
    cat gather.csv' | head -3
rm -f tests/fixtures/*.sig tests/fixtures/gather.csv
```

Expected: a CSV row whose `f_unique_weighted` is between 0.08 and 0.13. If it is not, the
fixture is wrong and the estimator test built on it will be meaningless — stop and fix it.
Do not widen the band to make it pass.

Measured during implementation: **0.11503**. (Higher than the 0.10 mixing ratio because a
150bp read contributes 120 distinct 31-mers while the reference is sketched at
`scaled=1000`; the weighting is over retained hashes, not reads.)

- [ ] **Step 4: Commit**

```bash
git add tests/fixtures/make_fixtures.py tests/fixtures/ref.fasta tests/fixtures/mix10.fastq.gz
git commit -m "test: add synthetic reference and known-fraction mixture fixtures"
```

---

## Task 4: SKETCH_REFERENCE

Runs once per pipeline. Emits the sketch plus a measured genome size, and fails loudly on a
file that is not a nucleotide FASTA.

**Files:**
- Create: `modules/local/sketch_reference.nf`
- Test: `modules/local/sketch_reference_tests/main.nf.test`

- [ ] **Step 1: Write the failing test**

Create `modules/local/sketch_reference_tests/main.nf.test`:

```groovy
nextflow_process {

    name "Test Process SKETCH_REFERENCE"
    script "../sketch_reference.nf"
    process "SKETCH_REFERENCE"
    tag "sketch_reference"

    test("sketches a reference and measures genome size") {
        when {
            process {
                """
                input[0] = file("\${projectDir}/tests/fixtures/ref.fasta", checkIfExists: true)
                """
            }
        }
        then {
            assert process.success
            def stats = new groovy.json.JsonSlurper().parseText(
                file(process.out.stats[0]).text)
            assert stats.genomeSize == 200000
            assert file(process.out.sig[0]).size() > 0
        }
    }

    test("fails on a file with no FASTA headers") {
        when {
            process {
                """
                input[0] = file("\${projectDir}/tests/fixtures/not_a_fasta.gff", checkIfExists: true)
                """
            }
        }
        then {
            assert process.failed
        }
    }
}
```

- [ ] **Step 2: Create the negative-case fixture**

```bash
cat > tests/fixtures/not_a_fasta.gff <<'GFF'
##gff-version 3
LmjF.01	VEuPathDB	gene	1000	2000	.	+	.	ID=LmjF.01.0010
LmjF.01	VEuPathDB	gene	3000	4000	.	-	.	ID=LmjF.01.0020
GFF
```

- [ ] **Step 3: Run test to verify it fails**

```bash
nf-test test modules/local/sketch_reference_tests/main.nf.test
```

Expected: FAIL — `../sketch_reference.nf` does not exist.

- [ ] **Step 4: Write the implementation**

Create `modules/local/sketch_reference.nf`:

```groovy
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
```

- [ ] **Step 5: Run test to verify it passes**

```bash
nf-test test modules/local/sketch_reference_tests/main.nf.test
```

Expected: PASS, 2 tests.

- [ ] **Step 6: Commit**

```bash
git add modules/local/sketch_reference.nf modules/local/sketch_reference_tests/main.nf.test tests/fixtures/not_a_fasta.gff
git commit -m "feat: add SKETCH_REFERENCE with FASTA validation"
```

---

## Task 5: MEASURE_SAMPLE

Produces every measured fact about a sample: on-target fraction, total reads, read length,
pilot size. The `--threshold-bp 0` flag is load-bearing — without it `gather` writes no result
row for low-fraction samples and the estimate silently collapses to zero.

For paired-end data only read 1 is measured. `totalReads` therefore counts *pairs*, matching
the existing convention in `subsample_fastq.nf`.

**Files:**
- Create: `modules/local/measure_sample.nf`
- Test: `modules/local/measure_sample_tests/main.nf.test`

- [ ] **Step 1: Write the failing test**

Create `modules/local/measure_sample_tests/main.nf.test`:

```groovy
nextflow_process {

    name "Test Process MEASURE_SAMPLE"
    script "../measure_sample.nf"
    process "MEASURE_SAMPLE"
    tag "measure_sample"

    test("estimates a known 10 percent mixture") {
        setup {
            run("SKETCH_REFERENCE") {
                script "../sketch_reference.nf"
                process {
                    """
                    input[0] = file("\${projectDir}/tests/fixtures/ref.fasta", checkIfExists: true)
                    """
                }
            }
        }
        when {
            params { pilotSize = 100000 }
            process {
                """
                input[0] = Channel.of([ [id: 'mix10'],
                    file("\${projectDir}/tests/fixtures/mix10.fastq.gz", checkIfExists: true) ])
                input[1] = SKETCH_REFERENCE.out.sig
                """
            }
        }
        then {
            assert process.success
            def m = new groovy.json.JsonSlurper().parseText(
                file(process.out.metrics[0][1]).text)
            assert m.totalReads == 10000
            assert m.readLength == 150
            assert m.pilotReads == 10000
            assert m.onTargetFraction > 0.08
            assert m.onTargetFraction < 0.13
        }
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
nf-test test modules/local/measure_sample_tests/main.nf.test
```

Expected: FAIL — `../measure_sample.nf` does not exist.

- [ ] **Step 3: Write the implementation**

Create `modules/local/measure_sample.nf`:

```groovy
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
nf-test test modules/local/measure_sample_tests/main.nf.test
```

Expected: PASS. The reported `onTargetFraction` should be near 0.10.

- [ ] **Step 5: Commit**

```bash
git add modules/local/measure_sample.nf modules/local/measure_sample_tests/main.nf.test
git commit -m "feat: add MEASURE_SAMPLE on-target fraction estimator"
```

---

## Task 6: CONCATENATE_FASTQ symlink short-circuit

Single-file samples currently get fully decompressed and recompressed to produce a
byte-equivalent file. Concat is the pipeline's peak-disk high-water mark, so this is the
cheapest available win.

**Files:**
- Modify: `modules/local/concatenate_fastq.nf`
- Test: `modules/local/concatenate_fastq_tests/main.nf.test`

- [ ] **Step 1: Write the failing test**

Create `modules/local/concatenate_fastq_tests/main.nf.test`:

```groovy
nextflow_process {

    name "Test Process CONCATENATE_FASTQ"
    script "../concatenate_fastq.nf"
    process "CONCATENATE_FASTQ"
    tag "concatenate_fastq"

    test("single-end single-file input is passed through without recompression") {
        when {
            process {
                """
                input[0] = Channel.of([ [id: 'solo', hasPairedReads: false],
                    [ file("\${projectDir}/tests/fixtures/mix10.fastq.gz", checkIfExists: true) ] ])
                """
            }
        }
        then {
            assert process.success
            def out = file(process.out.reads[0][1] instanceof List
                ? process.out.reads[0][1][0] : process.out.reads[0][1])
            assert out.name == 'solo_concat.fastq.gz'
            // passthrough must be byte-identical to the input, not a re-gzip
            assert out.bytes == file("\${projectDir}/tests/fixtures/mix10.fastq.gz").bytes
        }
    }

    test("paired-end single-run input is passed through without recompression") {
        when {
            process {
                """
                input[0] = Channel.of([ [id: 'pair', hasPairedReads: true],
                    [ file("\${projectDir}/tests/fixtures/pair_1.fastq.gz", checkIfExists: true),
                      file("\${projectDir}/tests/fixtures/pair_2.fastq.gz", checkIfExists: true) ] ])
                """
            }
        }
        then {
            assert process.success
            def outs = process.out.reads[0][1].collect { file(it) }.sort { it.name }
            assert outs*.name == ['pair_concat_1.fastq.gz', 'pair_concat_2.fastq.gz']
            assert outs[0].bytes == file("\${projectDir}/tests/fixtures/pair_1.fastq.gz").bytes
            assert outs[1].bytes == file("\${projectDir}/tests/fixtures/pair_2.fastq.gz").bytes
        }
    }

    test("single-end multi-file input is concatenated") {
        when {
            process {
                """
                input[0] = Channel.of([ [id: 'duo', hasPairedReads: false],
                    [ file("\${projectDir}/tests/fixtures/mix10.fastq.gz", checkIfExists: true),
                      file("\${projectDir}/tests/fixtures/mix10.fastq.gz", checkIfExists: true) ] ])
                """
            }
        }
        then {
            assert process.success
            def out = path(process.out.reads[0][1] instanceof List
                ? process.out.reads[0][1][0] : process.out.reads[0][1])
            // two copies of a 10000-read file
            assert out.linesGzip.size() / 4 == 20000
        }
    }
}
```

The paired fixtures do not exist yet. Create them first:

```bash
python3 - <<'PY'
import gzip, os, shutil
src = 'tests/fixtures/mix10.fastq.gz'
for mate in ('1', '2'):
    dst = 'tests/fixtures/pair_%s.fastq.gz' % mate
    with gzip.open(src, 'rt') as i, gzip.open(dst, 'wt') as o:
        for n, line in enumerate(i):
            o.write(line if n %% 4 else line.rstrip('\n') + '/' + mate + '\n')
    print(dst, os.path.getsize(dst))
PY
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
nf-test test modules/local/concatenate_fastq_tests/main.nf.test
```

Expected: the two passthrough tests FAIL on the byte-identity assertion — the current
implementation decompresses and re-gzips, which changes the gzip stream. The multi-file test
passes (that path is unchanged).

- [ ] **Step 3: Add the short-circuit**

In `modules/local/concatenate_fastq.nf`, insert this block immediately after the
`def file_list = ...` line and before the `if (!meta.hasPairedReads)` branch:

```groovy
    // A single already-gzipped run needs no merging; symlinking avoids a pointless
    // decompress/recompress cycle at the pipeline's peak-disk step. Guarded so that a
    // malformed sample (paired flag but only one file) still falls through to the real
    // concatenation path, which reports the error.
    def all_gzipped = file_list.every { it.name.endsWith('.gz') }

    if (!meta.hasPairedReads && file_list.size() == 1 && all_gzipped) {
        return """
        ln -s ${file_list[0]} ${meta.id}_concat.fastq.gz
        """
    }

    if (meta.hasPairedReads && file_list.size() == 2 && all_gzipped) {
        def r1 = file_list.find { it.name.contains('_1.fastq') || it.name.contains('_R1') }
        def r2 = file_list.find { it.name.contains('_2.fastq') || it.name.contains('_R2') }
        if (r1 && r2) {
            return """
            ln -s ${r1} ${meta.id}_concat_1.fastq.gz
            ln -s ${r2} ${meta.id}_concat_2.fastq.gz
            """
        }
    }
```

Single-run paired-end is the common SRA case, so this second branch is where most of the
saving actually lands.

- [ ] **Step 4: Run tests to verify they pass**

```bash
nf-test test modules/local/concatenate_fastq_tests/main.nf.test
```

Expected: PASS, 3 tests.

- [ ] **Step 5: Commit**

```bash
git add modules/local/concatenate_fastq.nf modules/local/concatenate_fastq_tests/main.nf.test
git commit -m "perf: symlink single-run samples instead of recompressing"
```

---

## Task 7: SUBSAMPLE_FASTQ takes a per-sample target

Read counting moves upstream to `MEASURE_SAMPLE`, so this process becomes pure execution: it
is told how many reads to keep out of how many, and it keeps them.

**Files:**
- Modify: `modules/local/subsample_fastq.nf`
- Test: `modules/local/subsample_fastq_tests/main.nf.test`

- [ ] **Step 1: Write the failing test**

Create `modules/local/subsample_fastq_tests/main.nf.test`:

```groovy
nextflow_process {

    name "Test Process SUBSAMPLE_FASTQ"
    script "../subsample_fastq.nf"
    process "SUBSAMPLE_FASTQ"
    tag "subsample_fastq"

    test("cuts to the requested number of reads") {
        when {
            process {
                """
                input[0] = Channel.of([ [id: 'cut', hasPairedReads: false],
                    file("\${projectDir}/tests/fixtures/mix10.fastq.gz", checkIfExists: true),
                    2000L, 10000L ])
                """
            }
        }
        then {
            assert process.success
            def out = path(process.out.reads[0][1] instanceof List
                ? process.out.reads[0][1][0] : process.out.reads[0][1])
            def n = out.linesGzip.size() / 4
            assert n > 1800 && n < 2200
        }
    }

    test("passes through untouched when target exceeds available reads") {
        when {
            process {
                """
                input[0] = Channel.of([ [id: 'keep', hasPairedReads: false],
                    file("\${projectDir}/tests/fixtures/mix10.fastq.gz", checkIfExists: true),
                    50000L, 10000L ])
                """
            }
        }
        then {
            assert process.success
            def out = path(process.out.reads[0][1] instanceof List
                ? process.out.reads[0][1][0] : process.out.reads[0][1])
            assert out.linesGzip.size() / 4 == 10000
        }
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
nf-test test modules/local/subsample_fastq_tests/main.nf.test
```

Expected: FAIL — the process currently declares two separate inputs
(`tuple val(meta), path(reads)` and `val max_reads`), not a 4-element tuple.

- [ ] **Step 3: Rewrite the process**

Replace the entire contents of `modules/local/subsample_fastq.nf` with:

```groovy
process SUBSAMPLE_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    container 'staphb/seqtk:1.4'

    publishDir params.outDir, mode: 'copy'

    input:
    tuple val(meta), path(reads), val(target_reads), val(total_reads)

    output:
    tuple val(meta), path("${meta.id}*_subsampled.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def seed = 42
    def read_list = reads instanceof List ? reads : [reads]

    if (meta.hasPairedReads) {
        """
        if [ ${total_reads} -gt ${target_reads} ]; then
            fraction=\$(awk -v t="${target_reads}" -v n="${total_reads}" 'BEGIN {printf "%.10f", t/n}')
            seqtk sample -s ${seed} ${read_list[0]} \$fraction | gzip > ${meta.id}_1_subsampled.fastq.gz
            seqtk sample -s ${seed} ${read_list[1]} \$fraction | gzip > ${meta.id}_2_subsampled.fastq.gz
        else
            ln -s ${read_list[0]} ${meta.id}_1_subsampled.fastq.gz
            ln -s ${read_list[1]} ${meta.id}_2_subsampled.fastq.gz
        fi
        """
    } else {
        """
        if [ ${total_reads} -gt ${target_reads} ]; then
            fraction=\$(awk -v t="${target_reads}" -v n="${total_reads}" 'BEGIN {printf "%.10f", t/n}')
            seqtk sample -s ${seed} ${read_list[0]} \$fraction | gzip > ${meta.id}_subsampled.fastq.gz
        else
            ln -s ${read_list[0]} ${meta.id}_subsampled.fastq.gz
        fi
        """
    }

    stub:
    if (meta.hasPairedReads) {
        """
        touch ${meta.id}_1_subsampled.fastq.gz
        touch ${meta.id}_2_subsampled.fastq.gz
        """
    } else {
        """
        touch ${meta.id}_subsampled.fastq.gz
        """
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
nf-test test modules/local/subsample_fastq_tests/main.nf.test
```

Expected: PASS, 2 tests.

- [ ] **Step 5: Commit**

```bash
git add modules/local/subsample_fastq.nf modules/local/subsample_fastq_tests/main.nf.test
git commit -m "refactor: SUBSAMPLE_FASTQ takes per-sample target instead of counting reads"
```

---

## Task 8: PREPARE_SAMPLES subworkflow and wiring

Collapses the duplicated tail in `main.nf` and `retrieve_from_sra.nf` into one subworkflow,
and inserts measurement and policy into it.

**Files:**
- Create: `workflows/prepare_samples.nf`
- Modify: `workflows/retrieve_from_sra.nf`
- Modify: `main.nf`
- Modify: `nextflow.config`

- [ ] **Step 1: Write the subworkflow**

Create `workflows/prepare_samples.nf`:

```groovy
include { FORMAT_INPUT_FROM_SRA } from '../modules/local/format_input_from_sra'
include { CONCATENATE_FASTQ     } from '../modules/local/concatenate_fastq'
include { SUBSAMPLE_FASTQ       } from '../modules/local/subsample_fastq'
include { MEASURE_SAMPLE        } from '../modules/local/measure_sample'
include { depthPlan             } from '../modules/local/depth_policy'

workflow PREPARE_SAMPLES {

    take:
    grouped_reads      // [ meta, [files] ]
    reference_sig      // path
    policy             // map: assayType, genomeSize, targetCoverage, minOnTargetFraction, minPlausibleFraction

    main:
    CONCATENATE_FASTQ(grouped_reads)

    MEASURE_SAMPLE(CONCATENATE_FASTQ.out.reads, reference_sig)

    // Join on sample id rather than on the meta map, so the join key stays stable
    // even if a meta field is mutated upstream.
    reads_by_id = CONCATENATE_FASTQ.out.reads.map { meta, reads -> [ meta.id, meta, reads ] }
    metrics_by_id = MEASURE_SAMPLE.out.metrics.map { meta, json -> [ meta.id, json ] }

    // `policy` arrives as a value channel (genome size is measured, not known up front),
    // so it must be combined into the stream rather than dereferenced directly.
    plans = reads_by_id
        .join(metrics_by_id)
        .combine(policy)
        .map { id, meta, reads, json, policyMap ->
            def metrics = new groovy.json.JsonSlurper().parse(json.toFile())
            def plan = depthPlan(metrics, policyMap)
            if (plan.flagged) {
                log.warn "Sample ${id}: on-target fraction ${metrics.onTargetFraction} is below " +
                         "minPlausibleFraction (${policyMap.minPlausibleFraction}). If every sample " +
                         "is flagged, check that --referenceFasta is the right organism."
            }
            return [ meta, reads, metrics, plan ]
        }

    SUBSAMPLE_FASTQ(
        plans.map { meta, reads, metrics, plan -> [ meta, reads, plan.rawReads, metrics.totalReads ] }
    )

    FORMAT_INPUT_FROM_SRA(SUBSAMPLE_FASTQ.out.reads)

    formatted = FORMAT_INPUT_FROM_SRA.out.samplesheet
        .collectFile(keepHeader: true, storeDir: params.outDir, name: params.samplesheetName)

    metrics_header = "sample,on_target_fraction,total_reads,raw_reads_used," +
                     "estimated_coverage,read_length,pilot_reads,flagged\n"

    metrics_csv = plans
        .map { meta, reads, metrics, plan ->
            def coverage = plan.estimatedCoverage == null
                ? ''
                : String.format('%.2f', plan.estimatedCoverage)
            "${meta.id},${metrics.onTargetFraction},${metrics.totalReads},${plan.rawReads}," +
            "${coverage},${metrics.readLength},${metrics.pilotReads},${plan.flagged}\n"
        }
        .collectFile(name: 'sample_metrics.csv', storeDir: params.outDir,
                     seed: metrics_header, sort: true)

    emit:
    formattedInput = formatted
    sampleMetrics  = metrics_csv
    flags          = plans.map { meta, reads, metrics, plan -> plan.flagged }
}
```

- [ ] **Step 2: Rewrite retrieve_from_sra to delegate**

Replace the contents of `workflows/retrieve_from_sra.nf` with:

```groovy
include { SRATOOLS_FASTERQDUMP } from '../modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH    } from '../modules/nf-core/sratools/prefetch/main'
include { PREPARE_SAMPLES      } from './prepare_samples'

workflow RETRIEVE_FROM_SRA {

    take:
    samples
    reference_sig
    policy

    main:
    individual_sra_samples = samples.flatMap { meta, sra_ids ->
        sra_ids.collect { sra_id -> [ meta, sra_id ] }
    }

    SRATOOLS_PREFETCH(individual_sra_samples, [], [])
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra, [], [])

    grouped_reads = SRATOOLS_FASTERQDUMP.out.reads
        .map { meta, reads ->
            meta.hasPairedReads = reads.size() == 2
            return [ meta.id, meta, reads ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, read_lists ->
            def firstHasPairedReads = metas[0].hasPairedReads
            if (!metas.every { it.hasPairedReads == firstHasPairedReads }) {
                throw new IllegalStateException("SRR samples must be all paired or unpaired. Mixed results found")
            }
            return [ metas[0], read_lists.flatten() ]
        }

    PREPARE_SAMPLES(grouped_reads, reference_sig, policy)

    emit:
    formattedInput = PREPARE_SAMPLES.out.formattedInput
    flags          = PREPARE_SAMPLES.out.flags
}
```

- [ ] **Step 3: Rewrite main.nf**

Replace the contents of `main.nf` with:

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RETRIEVE_FROM_SRA } from './workflows/retrieve_from_sra'
include { PREPARE_SAMPLES   } from './workflows/prepare_samples'
include { SKETCH_REFERENCE  } from './modules/local/sketch_reference'
include { EXPAND_SRX_IDS    } from './modules/local/expand_srx_ids'

def sampleFlags = Collections.synchronizedList([])

workflow {

    if (!params.referenceFasta) {
        error "--referenceFasta is required: the target organism FASTA used to estimate " +
              "on-target fraction. Subsampling cannot be coverage-aware without it."
    }
    def referenceFile = file(params.referenceFasta, checkIfExists: true)
    if (referenceFile.size() == 0) {
        error "--referenceFasta (${params.referenceFasta}) is empty."
    }

    SKETCH_REFERENCE(Channel.value(referenceFile))

    policy = SKETCH_REFERENCE.out.stats.map { statsFile ->
        def stats = new groovy.json.JsonSlurper().parse(statsFile.toFile())
        return [
            assayType           : params.assayType,
            genomeSize          : stats.genomeSize as long,
            targetCoverage      : params.targetCoverage as int,
            minOnTargetFraction : params.minOnTargetFraction as double,
            minPlausibleFraction: params.minPlausibleFraction as double
        ]
    }.first()

    reference_sig = SKETCH_REFERENCE.out.sig.first()

    samples = Channel.fromPath(params.input + "/" + params.samplesheetName).splitCsv(skip: 1)

    if (params.fromSra) {
        EXPAND_SRX_IDS(samples.map { row -> [ row[0], row[1], row[3] ?: "" ] })

        grouped_sra_samples = EXPAND_SRX_IDS.out.rows
            .splitCsv()
            .map { row -> [ row[0], [id: row[0], var1: row[3] ?: ""], row[1] ] }
            .groupTuple(by: 0)
            .map { sample_id, metas, sra_ids -> [ metas[0], sra_ids ] }

        RETRIEVE_FROM_SRA(grouped_sra_samples, reference_sig, policy)
        RETRIEVE_FROM_SRA.out.flags.subscribe { sampleFlags << it }
    }
    else {
        grouped_local_samples = samples.map { row ->
                def files = [ file(params.input + "/" + row[1], checkIfExists: true) ]
                boolean hasPairedReads = false
                if (row[2]) {
                    files.add(file(params.input + "/" + row[2], checkIfExists: true))
                    hasPairedReads = true
                }
                return [ row[0], [id: row[0], var1: row[3], hasPairedReads: hasPairedReads], files ]
            }
            .groupTuple(by: 0)
            .map { sample_id, metas, file_lists ->
                def firstHasPairedReads = metas[0].hasPairedReads
                if (!metas.every { it.hasPairedReads == firstHasPairedReads }) {
                    throw new IllegalStateException("Samples must be all paired or unpaired. Mixed results found")
                }
                return [ metas[0], file_lists.flatten() ]
            }

        PREPARE_SAMPLES(grouped_local_samples, reference_sig, policy)
        PREPARE_SAMPLES.out.flags.subscribe { sampleFlags << it }
    }
}

workflow.onComplete {
    if (sampleFlags && sampleFlags.every { it }) {
        log.error "All ${sampleFlags.size()} samples fell below minPlausibleFraction " +
                  "(${params.minPlausibleFraction}). This usually means --referenceFasta " +
                  "(${params.referenceFasta}) is not the organism these reads came from. " +
                  "Check ${params.outDir}/sample_metrics.csv."
    }
}
```

- [ ] **Step 4: Update nextflow.config**

In `nextflow.config`, replace the `genomeSize` line and add the new params so the `params`
block reads:

```groovy
params {
  input = "$launchDir/data/"
  samplesheetName = "samplesheet.csv"
  fromSra = true
  outDir = "$launchDir/ngs-samples-output"

  // Target organism FASTA. Required: on-target fraction cannot be estimated without it.
  referenceFasta = null

  // Subsampling parameters
  assayType = "DNASeq" // RNASeq uses a fixed 20M read target; others use targetCoverage
  targetCoverage = 60

  // Fraction floor doubles as the inflation cap: 0.05 means never retain more than 20x
  // what a clean sample would need.
  minOnTargetFraction = 0.05
  // Below this a sample is flagged; all samples flagged usually means the wrong reference.
  minPlausibleFraction = 0.01
  pilotSize = 100000

  // Maximum SRA file size prefetch will download. Raise this if prefetch skips a
  // large run with "is larger than maximum allowed" (default sra-tools limit is 20G).
  maxDownloadSize = "50G"
}
```

- [ ] **Step 5: Verify the DAG builds with stubs**

```bash
mkdir -p /tmp/nf-smoke/data
printf 'sample,fastq_1,fastq_2,var1\nsolo,mix10.fastq.gz,,x\n' > /tmp/nf-smoke/data/samplesheet.csv
cp tests/fixtures/mix10.fastq.gz /tmp/nf-smoke/data/
nextflow run main.nf -stub-run \
  --fromSra false \
  --input /tmp/nf-smoke/data \
  --referenceFasta "$PWD/tests/fixtures/ref.fasta" \
  --outDir /tmp/nf-smoke/out
```

Expected: completes successfully; `SKETCH_REFERENCE`, `CONCATENATE_FASTQ`, `MEASURE_SAMPLE`,
`SUBSAMPLE_FASTQ`, and `FORMAT_INPUT_FROM_SRA` all appear in the trace.

- [ ] **Step 6: Run the pipeline for real on the fixture**

```bash
rm -rf /tmp/nf-smoke/out /tmp/nf-smoke/work
nextflow run main.nf \
  --fromSra false \
  --input /tmp/nf-smoke/data \
  --referenceFasta "$PWD/tests/fixtures/ref.fasta" \
  --outDir /tmp/nf-smoke/out \
  -w /tmp/nf-smoke/work
cat /tmp/nf-smoke/out/sample_metrics.csv
cat /tmp/nf-smoke/out/samplesheet.csv
```

Expected: `sample_metrics.csv` has a header plus one row for `solo` with
`on_target_fraction` near 0.10 and `flagged` false. `samplesheet.csv` still has exactly the
four columns `sample,fastq_1,fastq_2,var1`.

- [ ] **Step 7: Verify the samplesheet contract is unchanged**

```bash
head -1 /tmp/nf-smoke/out/samplesheet.csv
```

Expected exactly: `sample,fastq_1,fastq_2,var1`

- [ ] **Step 8: Commit**

```bash
git add workflows/prepare_samples.nf workflows/retrieve_from_sra.nf main.nf nextflow.config
git commit -m "feat: wire contamination-aware subsampling through shared PREPARE_SAMPLES"
```

---

## Task 9 (OPTIONAL — decide before implementing): paired-end mate accounting

**Read this before doing it.** The existing pipeline counts read 1 only, so for paired-end
data `totalReads` means *pairs*, and the coverage calculation treats each pair as contributing
`readLength` bases when it actually contributes `2 * readLength`. The result is that
paired-end samples retain roughly twice the data needed to hit `targetCoverage`.

This is a pre-existing inconsistency, not something this feature introduced, and it is not in
the spec. Fixing it will halve retained data for every paired-end sample, which changes
outputs for existing studies. **Skip this task unless John has explicitly approved it.**

**Files:**
- Modify: `modules/local/depth_policy.nf`
- Modify: `modules/local/depth_policy_tests/main.nf.test`
- Modify: `workflows/prepare_samples.nf`

- [ ] **Step 1: Add the failing test**

Append to `modules/local/depth_policy_tests/main.nf.test`, inside the `nextflow_function` block:

```groovy
    test("paired reads contribute two mates worth of bases") {
        function "depthPlan"
        when {
            function {
                """
                input[0] = [onTargetFraction: 1.0d, totalReads: 200000000L, readLength: 150, pilotReads: 100000]
                input[1] = [assayType: 'DNASeq', genomeSize: 33000000L, targetCoverage: 60,
                            minOnTargetFraction: 0.05d, minPlausibleFraction: 0.01d,
                            matesPerRead: 2]
                """
            }
        }
        then {
            assert function.success
            // 33Mb * 60 / (150 * 2) = 6,600,000 pairs
            assert function.result.targetOnTarget == 6600000L
        }
    }
```

- [ ] **Step 2: Run to verify it fails**

```bash
nf-test test modules/local/depth_policy_tests/main.nf.test
```

Expected: FAIL — `targetOnTarget` is 13200000, not 6600000.

- [ ] **Step 3: Honor matesPerRead in the policy**

In `modules/local/depth_policy.nf`, change `targetOnTargetReads` to:

```groovy
def targetOnTargetReads(Map policy) {
    if (policy.assayType == "RNASeq") {
        return 20000000L
    }
    int mates = (policy.matesPerRead ?: 1) as int
    long raw = (long) ((policy.genomeSize * policy.targetCoverage) / (policy.readLength * mates))
    return Math.max(1000000L, Math.min(100000000L, raw))
}
```

and change the coverage line in `depthPlan` to:

```groovy
    int mates = (policy.matesPerRead ?: 1) as int
    Double estimatedCoverage = policy.assayType == "RNASeq"
        ? null
        : (rawReads * observed * metrics.readLength * mates) / (policy.genomeSize as double)
```

- [ ] **Step 4: Pass matesPerRead per sample**

In `workflows/prepare_samples.nf`, change the `depthPlan` call inside the `plans` map to:

```groovy
            def samplePolicy = policy + [matesPerRead: meta.hasPairedReads ? 2 : 1]
            def plan = depthPlan(metrics, samplePolicy)
```

- [ ] **Step 5: Run tests to verify all pass**

```bash
nf-test test modules/local/depth_policy_tests/main.nf.test
```

Expected: PASS, 10 tests. The nine existing tests pass unchanged because `matesPerRead`
defaults to 1 when absent.

- [ ] **Step 6: Commit**

```bash
git add modules/local/depth_policy.nf modules/local/depth_policy_tests/main.nf.test workflows/prepare_samples.nf
git commit -m "fix: account for both mates when computing paired-end coverage target"
```

---

## Task 10: Documentation

`CLAUDE.md` currently documents `nextflow test <path>`, which is not a real command, and
documents `genomeSize`, which no longer exists.

**Files:**
- Modify: `CLAUDE.md`
- Modify: `README.md`

- [ ] **Step 1: Fix the testing section in CLAUDE.md**

Replace the `### Testing` code block in `CLAUDE.md` with:

```bash
# Requires nf-test on PATH (install: curl -fsSL https://code.askimed.com/install/nf-test | bash)

# Run a specific test file
nf-test test modules/local/depth_policy_tests/main.nf.test

# Run by tag
nf-test test --tag depth_policy

# NOTE: the vendored modules/nf-core/sratools/** tests require
# params.modules_testdata_base_path and do not currently run in this repo.
```

- [ ] **Step 2: Update the parameters section in CLAUDE.md**

Replace the `#### Subsampling Parameters` block with:

```markdown
#### Subsampling Parameters
- `referenceFasta`: **Required.** Target organism FASTA. Used to estimate each sample's
  on-target fraction so subsampling targets on-target reads rather than raw reads.
- `assayType`: "DNASeq", "RNASeq", or "ChipSeq" (default: `"DNASeq"`)
- `targetCoverage`: Coverage target for non-RNASeq assays (default: `60`)
- `minOnTargetFraction`: Fraction floor, which doubles as the inflation cap (default: `0.05`,
  i.e. never retain more than 20x a clean sample's requirement)
- `minPlausibleFraction`: Below this a sample is flagged (default: `0.01`). All samples
  flagged usually means the wrong `referenceFasta`.
- `pilotSize`: Reads drawn to estimate contamination (default: `100000`)

`genomeSize` has been removed — genome size is now measured from `referenceFasta`.
```

- [ ] **Step 3: Update the process flow section in CLAUDE.md**

Replace the `### Process Flow` block with:

```markdown
### Process Flow

`SKETCH_REFERENCE` runs once per pipeline, then both modes converge on `PREPARE_SAMPLES`:

1. **SRA Mode**: `samples` → `EXPAND_SRX_IDS` → group → `SRATOOLS_PREFETCH` →
   `SRATOOLS_FASTERQDUMP` → `PREPARE_SAMPLES`
2. **Local Mode**: `samples` → group → `PREPARE_SAMPLES`

`PREPARE_SAMPLES` = `CONCATENATE_FASTQ` → `MEASURE_SAMPLE` → depth policy →
`SUBSAMPLE_FASTQ` → `FORMAT_INPUT_FROM_SRA`.
```

- [ ] **Step 4: Document outputs in README.md**

Add this section to `README.md`:

```markdown
## Outputs

- `samplesheet.csv` — `sample,fastq_1,fastq_2,var1`. This contract is stable; downstream
  workflows can rely on the column set.
- `sample_metrics.csv` — per-sample measurements:
  `sample,on_target_fraction,total_reads,raw_reads_used,estimated_coverage,read_length,pilot_reads,flagged`.
  `estimated_coverage` is genome-relative and is left empty for RNASeq, where the depth
  target is a fixed read count.
- Subsampled FASTQ files.

## Contamination-aware subsampling

Subsampling targets *on-target* reads, not raw reads. Each sample's on-target fraction is
estimated from a 100k-read pilot using k-mer containment against `--referenceFasta`, and the
number of raw reads retained is inflated by that fraction. A sample that is 15% target
therefore retains ~6.7x more raw reads than a clean one, and both reach the requested
coverage after alignment.

If `sample_metrics.csv` shows every sample flagged, `--referenceFasta` is almost certainly not
the organism the reads came from.
```

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md README.md
git commit -m "docs: document reference FASTA, metrics sidecar, and real test commands"
```

---

## Task 11: Full-suite verification

- [ ] **Step 1: Run every new test**

```bash
export PATH="$HOME/bin:$PATH"
nf-test test \
  modules/local/depth_policy_tests/main.nf.test \
  modules/local/sketch_reference_tests/main.nf.test \
  modules/local/measure_sample_tests/main.nf.test \
  modules/local/concatenate_fastq_tests/main.nf.test \
  modules/local/subsample_fastq_tests/main.nf.test
```

Expected: all PASS. Record the actual counts; do not claim success without reading the output.

- [ ] **Step 2: Confirm the stub run still builds both branches**

```bash
nextflow run main.nf -stub-run --fromSra false \
  --input /tmp/nf-smoke/data \
  --referenceFasta "$PWD/tests/fixtures/ref.fasta" \
  --outDir /tmp/nf-smoke/out-stub
```

Expected: exit 0.

- [ ] **Step 3: Confirm the missing-reference error is clear**

```bash
nextflow run main.nf --fromSra false --input /tmp/nf-smoke/data --outDir /tmp/nf-smoke/out-err
```

Expected: fails immediately with the `--referenceFasta is required` message, before any
process is submitted.

- [ ] **Step 4: Clean up scratch**

```bash
rm -rf /tmp/nf-smoke
```

- [ ] **Step 5: Push the branch**

```bash
git push -u origin contamination-aware-subsampling
```
