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
