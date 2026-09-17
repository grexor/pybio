# Command-line reference

Every `pybio` invocation prints its config file path and active genomes folder first, then runs the requested command. Run `pybio -help` for the built-in summary; this page is the fuller reference.

## Genomes

| Command | Description |
|---|---|
| `pybio <species>` | Shorthand for `pybio genome <species>` — download/prepare a genome directly. |
| `pybio genome <species> [version]` | Download and prepare an Ensembl genome (or import a custom one with `-fasta`/`-gtf`). See [Genomes](genomes.md). |
| `pybio species [text]` | List available Ensembl species, optionally filtered by a search term. Alias: `pybio search`. |
| `pybio path <species> [version]` | Print the FASTA/GTF/GFF3 file paths for an already-downloaded genome. |
| `pybio config [folder]` | Show (no argument) or change (`folder`) the genomes storage folder in `~/.pybio`. |

## Mapping

| Command | Description |
|---|---|
| `pybio star <species> r1.fastq.gz [r2.fastq.gz] output.bam` | Align reads to a genome's STAR index and produce a sorted, indexed BAM. See [Read mapping](mapping.md). |
| `pybio sam2bam input.sam output.bam` | Convert, sort and index a SAM file into BAM. |

## Other tools

| Command | Description |
|---|---|
| `pybio aimux -r1 ... -r2 ... -annotation ... -barcodes ... -stats ... -output ...` | Demultiplex FASTQ reads by barcode. See [Demultiplexing (aimux)](aimux.md). |
| `pybio gff4jbrowse input.gff output.gff` | Rewrite a GFF3 file for JBrowse2: drops gene records and moves `Parent=gene:` onto transcripts' `Name` property. |

## Global options

These apply to `pybio genome`/`pybio <species>` and, where relevant, `pybio star`:

| Option | Description |
|---|---|
| `-genome_version <v>` | Use a specific genome version instead of the latest Ensembl release. |
| `-fasta <file>`, `-gtf <file>` | Assembly/annotation files for importing a custom genome. |
| `-nostar` | Skip building the STAR index. |
| `-nosalmon` | Skip building the salmon index. |
| `-threads n` (or `-t n`) | Number of threads to use (default 1). |
| `-alignIntronMax n` | STAR's maximum intron size, passed through to `pybio star`. |
| `-genomeSAindexNbases`, `-genomeChrBinNbits` | Passed through to `STAR --runMode genomeGenerate` when building an index. |
| `-version` | Print the installed pybio version and exit. |
| `-help` | Print the built-in usage summary. |

Any option `pybio star` doesn't recognize is forwarded as-is to the underlying `STAR` command, so STAR-specific flags not listed above still work.
