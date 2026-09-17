# File formats

Lightweight readers for FASTA and FASTQ, and a couple of format-validation helpers. These are plain, dependency-free classes — useful when you want to stream a file record-by-record without pulling in a heavier library.

## FASTA

```python
import pybio.data
f = pybio.data.Fasta("sequences.fasta")  # also accepts .gz and .bz2
while f.read():
    print(f.id, len(f.sequence))
```

`read()` returns `False` at end of file. `f.id` is the header line (everything after `>`, untouched — split on whitespace yourself if you only want the accession), `f.sequence` is the concatenated sequence with newlines stripped.

## FASTQ

```python
f = pybio.data.Fastq("reads.fastq.gz")
while f.read():
    print(f.id, f.sequence, f.quality)
```

Same read-loop pattern as `Fasta`. `f.quality` is the raw quality string, unconverted — use `ord(c)` yourself to get Phred scores, or see `fastq_qminmax()` below for a ready-made min/max helper.

## Format checks

Two small helpers for validating files before you commit to processing them:

```python
pybio.data.fasta_check("sequences.fasta")
# (True, "File in FASTA format")

pybio.data.fastq_qminmax("reads.fastq.gz")
# (qmin, qmax) — the min and max Phred quality ord() values found in the file
```

`fasta_check(filename, allowed_chars=["A","C","T","G","N"])` reads the whole file and confirms every sequence uses only the allowed characters (and isn't empty), returning `(False, reason)` on the first violation. `fastq_qminmax()` scans a FASTQ file and returns the observed `(qmin, qmax)` quality ordinal range — useful for guessing a FASTQ's quality encoding (Phred+33 vs. Phred+64) before aligning.
