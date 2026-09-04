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

# Paired-end fixtures derived from the mix10 read set: same reads, with /1 and /2
# appended to the read NAME (header line only, not sequence/+/quality). Read back
# mix10's own lines so the mate files are guaranteed to match it read-for-read.
# Same mtime/filename pinning as mix10.fastq.gz above for determinism.
with gzip.open(os.path.join(HERE, "mix10.fastq.gz"), "rb") as f:
    mix10_lines = f.read().splitlines()


def write_mates(path, mate_suffix):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
            for n, line in enumerate(mix10_lines):
                if n % 4 == 0:
                    line = line + mate_suffix.encode()
                out.write(line + b"\n")


write_mates(os.path.join(HERE, "pair_1.fastq.gz"), "/1")
write_mates(os.path.join(HERE, "pair_2.fastq.gz"), "/2")

print("pair_1.fastq.gz / pair_2.fastq.gz: mix10 reads with /1 and /2 mate suffixes")

# gzip of the existing reference, for testing that SKETCH_REFERENCE accepts
# gzipped FASTA input. Same mtime/filename pinning as mix10.fastq.gz above.
with open(os.path.join(HERE, "ref.fasta.gz"), "wb") as raw:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
        with open(os.path.join(HERE, "ref.fasta"), "rb") as ref_in:
            out.write(ref_in.read())

print("ref.fasta.gz: gzip of ref.fasta")

# A protein FASTA, long enough to clear the old (pre-fix) 1000-character floor.
# A/C/G/T/N are all valid amino acid codes, so random protein sequence still
# has ~20-25% "nucleotide-looking" characters by chance - this is the fixture
# that proves the composition-ratio check (not just a length floor) is needed.
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
PROTEIN_LEN = 5000

protein_seq = "".join(random.choice(AA_ALPHABET) for _ in range(PROTEIN_LEN))

with open(os.path.join(HERE, "protein.fasta"), "w") as out:
    out.write(">synthetic_protein length=%d\n" % PROTEIN_LEN)
    for i in range(0, PROTEIN_LEN, 60):
        out.write(protein_seq[i:i + 60] + "\n")

print("protein.fasta: %d aa" % PROTEIN_LEN)

# Headers only, no sequence lines at all - exercises the guard around the
# genome-size grep, which otherwise dies on grep -v finding zero matches.
with open(os.path.join(HERE, "headers_only.fasta"), "w") as out:
    out.write(">seq1 description one\n")
    out.write(">seq2 description two\n")
    out.write(">seq3 description three\n")

print("headers_only.fasta: 3 headers, no sequence")
