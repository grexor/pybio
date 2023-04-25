<picture><img src="media/pybio.png" width="70"/></picture>
## pybio: basic genomics toolset

pybio is a Python framework for handling genomics operations and a direct interface to Ensembl genomes.

Downloading an Ensembl genome+annotation takes one line: `pybio genome homo_sapiens`.

Features include genome+annotation download from Ensembl and processing with STAR and salmon, support of Fasta, Fastq, bedGraph and other file formats, motif sequence searches, and specific APA feautes like alternative polyadenylation site-pair classification (same-exon, skipped-exon, composite-exon) and more.

## Contents

* [Installation](#installation)
* [Quick Start](#Quick-Start)
* [Documentation](#documentation)
  * [Downloading Ensembl genomes](#downloading-Ensembl-genomes)
  * [Retrieving genomic sequences](#retrieving-genomic-sequences)
  * [Annotate genomic position](#annotate-genomic-position)
  * [Importing genome annotation](#importing-genome-annotation)
  * [File formats](#file-formats)
  * [Genomic Coordinates](#Genomic-Coordinates)
* [Authors](#authors)
* [Reporting problems](#reporting-problems)

### Installation

The easiest way to install pybio is running:

`pip install pybio`

If you would like instead to install the latest developmental version from this repository:

```
# clone pybio repository locally
git clone https://github.com/grexor/pybio.git

# build and install
./build.sh
```

### Quick Start

pybio is strongly integrated with Ensembl and provides genomic loci search for diverse annotated features (genes -> transcripts -> exons + 5UTR + 3UTR).

Let's say we are interested in the human genome. First download and prepare the genome with a single command:

```
pybio genome homo_sapiens
```

Searching a genomic position for features is easy in python, for example:

```
import pybio
genes, transcripts, exons, UTR5, UTR3 = annotate("homo_sapiens", "1", "+", 11012344)
```

This will return a list of feature objects (genes, transctipts, exons, 3'-UTR and 5'-UTR) (check [core/genomes.py](core/genomes.py) classes to see details of these objects).

If you would like to know all genes that span the provided position, you could then write:

```
for gene in genes:
   print(gene.gene_id, gene.gene_name, gene.start, gene.stop)
```

And to list all transcripts of each gene, you could extend the code like this:

```
for gene in genes:
   print(gene.gene_id, gene.gene_name, gene.start, gene.stop)
   for transcript in gene.transcripts:
      print(transcript.transcript_id)
```

However you could also start directly with transcripts, and print to which genes are transcripts assigned to:

```
for transcript in transcripts:
  print(transcript.gene.gene_id, transcript.transcript_id)
```

And an intuitive graph representation of relationships between feature objects:

```
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

Representation of relationships between feature objects:

```
                gene = Gene instance object
    gene.transcripts = list of all transcript objects of the gene
          transcript = Transcript instance object
     transcript.gene = points to the gene of the transcript
    transcript.exons = list of all exon objects of the transcript
transcript.utr5/utr3 = points to the UTR5 / UTR3 of the transcript
                exon = Exon instance object
     exon.transcript = points to the transcript of the exon
           utr5/utr3 = Utr5 / Utr3 instance object
utr5/utr3.transcript = points to the transcript of the UTR5/UTR3
```

### Documentation

Here we provide basic `pybio` usage examples.

#### Downloading Ensembl genomes

To download Ensembl genomes simply run a few commands on the command line. For example:

```
$ pybio genome homo_sapiens      # downloads the latest version of Ensembl homo_sapiens assembly and annotation
$ pybio genome homo_sapiens 109  # downloads a specific version (in this case, v109) of Ensembl homo_sapiens assembly and annotation
$ pybio genome elephant          # will list all available elephant genomes and make you choose which one to download
```

The above will download the FASTA sequence and GTF annotation. If you have STAR and salmon installed on your system, pybio will also build an index of the genome for both.

Data will be stored in the folder specified in the file `pybio.config`. The genomes folder structure is as follows:

```
homo_sapiens.assembly.ensembl109             # FASTA files of the genome, each chromosome in a separate file
homo_sapiens.annotation.ensembl109           # Annotation in GTF and TAB format
homo_sapiens.assembly.ensembl109.star        # STAR index, GTF annotation aware
homo_sapiens.transcripts.ensembl109          # transcriptome, this is the Ensembl "cDNA" file in FASTA format
homo_sapiens.transcripts.ensembl109.salmon   # Salmon index of the transcriptome
```

pybio also supports the download of Ensembl Genomes (Ensembl Fungi, Ensembl Plants, Ensembl Protists, Ensembl Metazoa). You simply provide the name of the species on the command line to automagically download the genome, the assembly and prepare STAR and salmon indices.

For example, to download the latest version of the Dictyostelium discoideum genome, you would write:

```
$ pybio genome dicty                      # would search for genomes with "dicty" in the name of the species or description
$ pybio genome dictyostelium_discoideum   # also directly providing the exact genome species works
```

Another example is download the latest Arabidopsis thaliana genome:

```
$ pybio genome arabidopsis_thaliana
```

To see all available species, simply run `pybio species`. Moreover, to see all available `arabidopsis` genomes, you could run:

```
$ pybio species arabidopsis
arabidopsis_halleri	Ahal2.2	ensemblgenomes	plants	ensemblgenomes56
arabidopsis_lyrata	v.1.0	ensemblgenomes	plants	ensemblgenomes56
arabidopsis_thaliana	TAIR10	ensemblgenomes	plants	ensemblgenomes56
```

Voila.

#### Retrieving genomic sequences

To retrieve stretches of genomic sequence, we use the seq(genome, chr, strand, position, upstream, downstream) method:

```
import pybio
seq = pybio.core.genomes.seq("homo_sapiens", "1", "+", 450000, -20, 20)
```

The above command fetches the chr 1 sequence from 450000-20..450000+20, the resulting sequence is of length 41, `TACCCTGATTCTGAAACGAAAAAGCTTTACAAAATCCAAGA`.

#### Annotating genomic positions

Given a genomic position, we can quickly retrieve the gene, transcript, exon and utr5/3 information at the given position. If there are several features (genes, transcripts, exons, UTR regions) at the specified position, they are all returned.

```
# annotate position
genes, transcripts, exons, utr5, utr3 = pybio.genomes.annotate("hg38", "1", "+", 11012344)

# print all genes that cover the position
for gene in genes:
   print(gene.gene_id, gene.gene_name, gene.start, gene.stop)
```

The above command would return:

```
[pybio] loading genome annotation for homo_sapiens with Ensembl version 109
ENSG00000120948, TARDBP, 11012343, 11030527
```

You can also easily access all transcripts of each gene with `gene.transcripts` and all exons of each transcript with `transcript.exons`.

#### Dependencies

Basic dependencies include [pysam](https://pysam.readthedocs.io/en/latest/api.html), [numpy](https://numpy.org/) and [samtools](http://www.htslib.org) and should be installed automatically by pip when you install pybio over `pip install pybio`.

Optional dependencies include [STAR](https://github.com/alexdobin/STAR) and [salmon](https://combine-lab.github.io/salmon/getting_started/) if you would like to build genome/transcriptome indices and align reads.

### File formats

Supported file formats. Work in progress.

### Genomic Coordinates

All genomic coordinates we operate with inside pybio are 0-based left+right inclusive. This means, when we say for example 100-103, this would include coordinates 100, 101, 102 and 103. The first coordinate is 0.

**Important**

Refseq and Ensembl GTF files are 1-indexed. When we read files from refseq/ensembl, we substract 1 on all coordinates to keep this in line with other coordinate structures inside pybio (which are all 0-indexed).

### Authors

[pybio](https://github.com/grexor/pybio) is developed and supported by [Gregor Rot](https://grexor.github.io).

### Issues and Suggestions

Use the [issues page](https://github.com/grexor/pybio/issues) to report issues and leave suggestions.
