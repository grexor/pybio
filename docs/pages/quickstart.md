# Quick Start

pybio is tightly integrated with Ensembl and lets you search any genomic position for the annotated features that overlap it — genes, their transcripts, and each transcript's exons and 5'/3' UTRs.

First, download and prepare a genome (here, human):

```bash
pybio genome homo_sapiens
```

Then query a genomic position from Python with `pybio.core.genomes.annotate(species, chr, strand, pos)`:

```python
import pybio
result = pybio.core.genomes.annotate("homo_sapiens", "1", "+", 11012344)
genes, transcripts, exons, utr3, utr5 = result
```

This returns five lists of feature objects — genes, transcripts, exons, 3'-UTRs and 5'-UTRs — all overlapping the given position. (See [pybio/core/genomes.py](https://github.com/grexor/pybio/blob/master/pybio/core/genomes.py) for the full `Gene`/`Transcript`/`Exon`/`Utr5`/`Utr3` class definitions.)

!!! note
    Positions are 0-based, left+right inclusive. See [Genomic coordinates](coordinates.md).

To list the genes that span the position:

```python
for gene in genes:
    print(gene.gene_id, gene.gene_name, gene.start, gene.stop)
```

...and every transcript of each gene:

```python
for gene in genes:
    print(gene.gene_id, gene.gene_name, gene.start, gene.stop)
    for transcript in gene.transcripts:
        print(transcript.transcript_id)
```

You can also start from `transcripts` directly, and walk back up to the gene each one belongs to:

```python
for transcript in transcripts:
    print(transcript.gene.gene_id, transcript.transcript_id)
```

## Feature relationships

Feature objects link to each other in both directions, so you can navigate the annotation graph from whichever level is most convenient:

```text
gene <-> transcript_1 <-> exon_1
                          <-> exon_2
                          ...
                          <-> utr5
                          <-> utr3
         <-> transcript_2 <-> exon_1
                          <-> exon_2
                          ...
                          <-> utr5
                          <-> utr3
```

| Object | Attribute | Points to |
|---|---|---|
| `gene` | `.transcripts` | list of all transcripts of the gene |
| `transcript` | `.gene` | the gene the transcript belongs to |
| `transcript` | `.exons` | list of all exons of the transcript |
| `transcript` | `.utr5` / `.utr3` | the transcript's 5'-UTR / 3'-UTR object |
| `exon` | `.transcript` | the transcript the exon belongs to |
| `utr5` / `utr3` | `.transcript` | the transcript the UTR belongs to |

## Next steps

- [Genomes](genomes.md) — downloading more species, custom genomes, and where the data is stored.
- [Sequences & annotation](sequences.md) — pulling raw sequence for a region, and more on `annotate()`.
- [Read mapping](mapping.md) — aligning FASTQ reads to the genomes you've downloaded.
