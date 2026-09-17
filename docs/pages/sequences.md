# Sequences & annotation

Two core operations once you have a genome downloaded: pulling out raw sequence for a region, and resolving a genomic position to the annotated features that overlap it.

## Retrieving genomic sequence

`seq(species, chr, strand, pos, start, stop)` returns the sequence around `pos`, spanning `[pos+start .. pos+stop]`:

```python
import pybio
seq = pybio.core.genomes.seq("homo_sapiens", "1", "+", 450000, -20, 20)
```

This fetches chromosome 1 from `450000-20` to `450000+20` — 41 bases: `TACCCTGATTCTGAAACGAAAAAGCTTTACAAAATCCAAGA`.

On the minus strand, `seq()` returns the reverse complement, and `start`/`stop` are interpreted relative to the minus-strand reading direction (`gstart = pos-stop`, `gstop = pos-start`) so the result still reads 5'→3' on that strand.

If you already have absolute start/stop coordinates rather than a `pos` + offsets, call `seq_direct(species, chr, strand, start, stop)` directly — it's what `seq()` calls internally, and it's the faster option in a loop since it skips the offset arithmetic. Both are 0-based, left+right inclusive (see [Genomic coordinates](coordinates.md)); out-of-range positions are padded with `N` (configurable via the `flank` argument) rather than raising an error.

!!! note "How sequence lookups work"
    Under the hood, pybio doesn't re-parse the FASTA file for every lookup. When a genome is downloaded, each chromosome's sequence is written out as a flat, newline-free `.string` file, so `seq_direct()` can `seek()` straight to a position and read exactly the bytes it needs — no indexing library, no loading the chromosome into memory.

## Annotating genomic positions

`annotate(species, chr, strand, pos)` resolves a single position to every gene, transcript, exon and UTR that overlaps it:

```python
result = pybio.core.genomes.annotate("homo_sapiens", "1", "+", 11012344)
genes, transcripts, exons, utr3, utr5 = result
for gene in genes:
    print(gene.gene_id, gene.gene_name, gene.start, gene.stop)
# [pybio] loading genome annotation for homo_sapiens with genome version ensembl109
# ENSG00000120948, TARDBP, 11012343, 11030527
```

If several genes, transcripts or exons overlap the position — common with overlapping genes or alternative transcripts — all of them are returned. You can walk from a gene to its transcripts (`gene.transcripts`) and from a transcript to its exons (`transcript.exons`) as shown in [Quick Start](quickstart.md#feature-relationships).

The first call for a given `(species, genome_version)` pair loads and caches that genome's annotation in memory; subsequent calls for the same genome are fast. See the note below if you're calling `annotate()` for genomes other than the configured default.

!!! tip
    `annotate()` (and `seq()`/`seq_direct()`) default to `ensembl_version_latest` from your [config](genomes.md#configuration) when `genome_version` isn't given. If you work with more than one genome version, always pass `genome_version` explicitly to avoid ambiguity.

## Finding genes by name

`find_genes(search)` looks up genes by name against the currently loaded genome, returning `(gene_id, gene_name)` pairs — exact matches first, then substring matches:

```python
pybio.core.genomes.load("homo_sapiens")
pybio.core.genomes.find_genes("TARDBP")
# [('ENSG00000120948', 'TARDBP')]
```

Call `pybio.core.genomes.load(species, genome_version)` first if you haven't already called `annotate()` for that genome in this process — `find_genes()` searches whichever genome is currently loaded and does not load one itself.
