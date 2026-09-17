# pybio: basic genomics toolset

![pybio](assets/pybio_logo.png)

**pybio** is a Python framework for common genomics operations, built around a direct interface to Ensembl genome assemblies and annotations.

```bash
# install from PyPI
pip install pybio

# download the human genome (assembly + annotation)
pybio genome homo_sapiens
```

With a genome downloaded, you can query it straight from Python:

```python
import pybio
genes, transcripts, exons, utr3, utr5 = pybio.core.genomes.annotate("homo_sapiens", "1", "+", 11012344)
```

## What's included

- **Genome download & management** — fetch assemblies and annotations from Ensembl (or Ensembl Fungi/Plants/Protists/Metazoa) with one command, and manage custom genomes from your own FASTA/GTF files. See [Genomes](genomes.md).
- **Position annotation & sequence retrieval** — resolve a genomic position to its overlapping genes, transcripts, exons and UTRs, or pull out raw sequence for any region. See [Sequences & annotation](sequences.md).
- **Read mapping** — build STAR and salmon indices and align FASTQ reads to a genome with a single command. See [Read mapping](mapping.md).
- **Sequence & motif tools** — IUPAC-aware motif search, reverse complementation, signal smoothing and nucleotide composition plots. See [Sequence & motif tools](motifs.md).
- **CLIP / positional signal data** — load, cluster and query per-position bedGraph-style data (e.g. CLIP crosslink sites). See [CLIP & interval data](bedgraph.md).
- **FASTA/FASTQ/bedGraph file support** — lightweight readers and format checks. See [File formats](fileformats.md).
- **Barcode demultiplexing** — split multiplexed FASTQ files by sample using a barcode annotation table. See [Demultiplexing (aimux)](aimux.md).

## Where to start

New to pybio? Read [Installation](installation.md) and then [Quick Start](quickstart.md) — together they get you from a fresh install to your first annotated genomic position in a couple of minutes. Everything else in these docs is reference material to dig into as you need it.
