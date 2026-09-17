# Genomes

pybio downloads, prepares and locates genome data for you. This page covers Ensembl genomes, custom genomes, and where everything is stored.

## Downloading Ensembl genomes

```bash
# downloads the homo_sapiens assembly and annotation (latest Ensembl version)
pybio genome homo_sapiens

# downloads a specific Ensembl version
pybio genome homo_sapiens 109
```

This downloads the FASTA sequence and GTF annotation. If **STAR** and **salmon** are installed on your system, pybio also builds a STAR genome index and a salmon transcriptome index (skip either with `-nostar` / `-nosalmon`).

pybio also supports Ensembl Genomes (Fungi, Plants, Protists, Metazoa) through the same command — just provide the species name:

```bash
# search for genomes with "dicty" in the species name or description;
# if exactly one match is found, download it directly
pybio genome dicty

# the exact species id always works too
pybio genome dictyostelium_discoideum

# another Ensembl Genomes example
pybio genome arabidopsis_thaliana
```

### Finding a species

If you're not sure of the exact species id, search first:

```bash
pybio species arabidopsis
```
```text
arabidopsis_halleri   Ahal2.2    ensemblgenomes  plants  ensemblgenomes56
arabidopsis_lyrata    v.1.0      ensemblgenomes  plants  ensemblgenomes56
arabidopsis_thaliana  TAIR10     ensemblgenomes  plants  ensemblgenomes56
```

`pybio search <text>` is an alias for the same lookup. Run `pybio species` with no search term to see where the full list of available Ensembl species is cached (`ensembl.json`, inside your genomes folder).

If you pass an ambiguous or partial species name to `pybio genome`, pybio runs the same fuzzy match and either downloads the single hit or prints the list of candidates for you to choose from.

## Adding custom genomes

To register a genome that isn't on Ensembl, provide your own FASTA (assembly) and GTF (annotation) files along with a genome version label of your choosing:

```bash
pybio genome species_name -fasta sample.fasta -gtf sample.gtf -genome_version v1
```

`species_name` and `v1` (the genome version) become the identifiers you pass to every other pybio function, exactly like `"homo_sapiens"` and `"ensembl109"` for an Ensembl genome.

## Where the data lives

Genome data is stored under the folder configured as `genomes_folder` (see [Configuration](#configuration) below), organized per species and version:

```text
homo_sapiens.assembly.ensembl109           # FASTA files of the genome
homo_sapiens.annotation.ensembl109         # annotation, GTF and pybio's own TAB/pickle format
homo_sapiens.assembly.ensembl109.star      # STAR index, GTF-annotation aware
homo_sapiens.transcripts.ensembl109        # transcriptome, Ensembl cDNA FASTA
homo_sapiens.transcripts.ensembl109.salmon # salmon index of the transcriptome
```

`pybio path <species>` prints the resolved FASTA/GTF/GFF3 paths for a species you've already downloaded — handy for feeding other tools:

```bash
pybio path homo_sapiens
```
```text
/home/user/genomes/homo_sapiens.assembly.ensembl109/homo_sapiens.fasta
/home/user/genomes/homo_sapiens.annotation.ensembl109/homo_sapiens.gtf.gz
/home/user/genomes/homo_sapiens.annotation.ensembl109/homo_sapiens.gff3.gz
```

## Configuration

pybio reads its configuration from `~/.pybio`, created automatically on first run from a bundled template. The two settings that matter are:

```text
genomes_folder="~/genomes"
ensembl_version_latest=ensembl113
```

- `genomes_folder` — where all downloaded and custom genomes are stored.
- `ensembl_version_latest` — the Ensembl release used whenever a function or command is called without an explicit `genome_version`.

To change the genomes folder without hand-editing the file, run:

```bash
pybio config /path/to/genomes
```

This updates `~/.pybio` and prints the new location. Run `pybio config` with no argument to see the current config file path (printed at the top of every `pybio` command's output, along with the genomes folder in use).
