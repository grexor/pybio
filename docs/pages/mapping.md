# Read mapping

pybio wraps **STAR** for aligning reads to a genome, and can build the indices STAR and salmon need automatically.

!!! note
    STAR and/or salmon must be installed and on your `PATH`. See [Dependencies](dependencies.md).

## Building indices

Indices are built for you as part of [downloading a genome](genomes.md#downloading-ensembl-genomes):

```bash
pybio genome homo_sapiens          # also builds a STAR index and a salmon index
pybio genome homo_sapiens -nostar  # skip the STAR index
pybio genome homo_sapiens -nosalmon # skip the salmon index
```

The STAR index is written to `<species>.assembly.<genome_version>.star`, alongside the assembly, and the salmon index to `<species>.transcripts.<genome_version>.salmon` — see [Genomes](genomes.md#where-the-data-lives) for the full folder layout.

## Aligning reads with STAR

Once a genome has a STAR index, align FASTQ reads with `pybio star`:

```bash
# single-end
pybio star homo_sapiens r1.fastq.gz output.bam

# paired-end
pybio star homo_sapiens r1.fastq.gz r2.fastq.gz output.bam
```

This runs STAR against the genome's index, converts the resulting SAM to a sorted, indexed BAM (via `sam2bam`, below), and cleans up STAR's intermediate files, leaving:

```text
output.bam            # sorted, indexed alignment
output.stats.txt       # STAR's Log.final.out
output.log.txt         # STAR's Log.out
output.progress.txt    # STAR's Log.progress.out
output.splice.tab      # STAR's splice junctions (SJ.out.tab)
```

Useful flags:

- `-threads n` — number of threads (default 1).
- `-genome_version` — pick a specific downloaded genome version instead of the latest.
- `-alignIntronMax n` — cap STAR's maximum intron size.
- `-genomeSAindexNbases`, `-genomeChrBinNbits` — passed through to `STAR --runMode genomeGenerate` when building an index for a small genome.
- any other flag not recognized by `pybio star` is passed straight through to the underlying `STAR` command.

## Converting SAM to BAM

`pybio sam2bam` wraps the samtools view/sort/index pipeline (mapped reads only, `-F 4`):

```bash
pybio sam2bam input.sam output.bam -threads 4
```

## STARsolo for Visium spatial data

`pybio.core.genomes.starsolo()` runs STARsolo configured specifically for 10x Visium spatial RNA-seq — `CB_UMI_Simple` mode with the Visium barcode geometry (16 bp cell barcode + 12 bp UMI):

```python
import pybio
pybio.core.genomes.starsolo(
    genome_dir="/genomes/homo_sapiens.assembly.ensembl109.star",
    r2_files=["cdna_R2.fastq.gz"],       # cDNA reads
    r1_files=["barcode_umi_R1.fastq.gz"], # barcode + UMI reads
    whitelist="visium-v1.txt",            # barcode whitelist, one 16bp barcode per line
    out_prefix="/data/sample_",
    threads=8,
)
```

STAR writes a coordinate-sorted BAM to `{out_prefix}Aligned.sortedByCoord.out.bam`; the function returns STAR's exit code. This is a thin, opinionated wrapper for the Visium case — for other single-cell barcode geometries, call STAR's `--soloType` options directly.
