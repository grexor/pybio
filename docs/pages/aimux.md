# Demultiplexing (aimux)

`aimux` demultiplexes short-read sequencing data using a samples + barcodes annotation table. The input is a TAB-delimited file of samples and barcodes, for example (`annotation.tab`):

```text
sample_id   barcode1   barcode2
sample1     ATTCGT     ACC
sample2     AGGTCC     ATT
...
```

Having R1, R2, I1 and I2 fastq files, aimux can be run with:

```bash
pybio aimux -r1 r1.fastq.gz -r2 r2.fastq.gz -i1 i1.fastq.gz -i2 i2.fastq.gz \
    -annotation annotation.tab \
    -barcodes barcode1:i1:RRRRRR_0_m1,barcode2:i2:RRR_0 \
    -stats aimux.stats \
    -output samples
```

This matches the I1 sequence against a reversed `barcode1` of length 6 (`RRRRRR`), starting at position 0 in the I1 read (`_0`) and allowing 1 mismatch (`_m1`). At the same time, it checks the I2 sequence against `barcode2`, requiring a perfect match (no `_m` given) and again starting at position 0. Demultiplexed output is written under `samples` (`-output`), and matching statistics to `aimux.stats`.
