## pybio: basic genomics toolset

pybio is a Python 3 framework for basic genomics/transcriptomics operations/annotations. The package is a dependency to [apa](https://github.com/grexor/apa) (alternative polyadenylation) and [RNAmotifs2](https://github.com/grexor/rnamotifs2) (motif cluster analysis). The pybio package provides:

Features include automagical genome assemblies+annotation download and indexing from Ensembl and STAR, support of Fasta, Fastq, bedGraph and other file formats, motif sequence searches, and specific APA feautes like alternative polyadenylation site-pair classification (same-exon, skipped-exon, composite-exon) and more.

## Contents

* [Installation](#installation)
* [Quick Start](#quickstart)
* [Documentation](#documentation)
  * [Downloading Ensembl genomes](#downloading-Ensembl-genomes)
  * [Retrieving genomic sequence](#retrieving-genomic-sequence)
  * [Annotate genomic position](#annotate-genomic-position)
  * [Importing genome annotation](#importing-genome-annotation)
  * [File formats](#file-formats)
    * [bedGraph](#bedgraph)
* [Authors](#authors)
* [Reporting problems](#reporting-problems)

### Installation

20221218: pip package in preparation, [https://pypi.org/project/pybio](https://pypi.org/project/pybio).

Current installation instructions below:

#### Clone the GitHub repository

Temporarily the most direct way of installing pybio is to clone the repository:

```
git clone https://github.com/grexor/pybio.git
```

Enter the cloned pybio folder and copy the config file from the template:

```
cp pybio.config.template pybio.config
```

Finally, change your PATH and PYTHONPATH environment variables: for example, if you cloned to `/home/user/pybio`, you would add the paths like so:

```
export PYTHONPATH=$PYTHONPATH:/home/user # to import pybio in python with "import pybio"
export PATH=$PATH:/home/user/pybio       # to run pybio on the command line with "pybio"
```

Voila.

#### Dependencies

There are a few software tools pybio depends on. For a quick start, you don't need to have any of these dependencies installed.

* [STAR aligner](https://github.com/alexdobin/STAR), `sudo apt-get install rna-star`
* [pysam](https://pysam.readthedocs.io/en/latest/api.html), `sudo apt-get install python-pysam`
* [numpy](https://numpy.org/), `sudo apt-get install python-numpy`
* [Salmon](https://combine-lab.github.io/salmon/getting_started/), download and install from Salmon webpage
* [samtools](http://www.htslib.org), `sudo apt-get install samtools`

### Quick Start

In progress

### Documentation

Here we provide basic `pybio` usage examples.

#### Downloading Ensembl genomes

Downloading Ensembl genome assemblies and annotation is easy. Simply run these commands on the command line:

```
pybio ensembl homo_sapiens      # downloads the latest version of Ensembl homo_sapiens assembly and annotation
pybio ensembl homo_sapiens 109  # downloads a specific version of Ensembl homo_sapiens assembly and annotation
```

This will download the FASTA sequence, GTF if you have STAR and salmon installed, will also build an index of the genome for both.

Data will be stored in the folder specified in the file `pybio.config`. The genomes folder structure is as follows:

```
homo_sapiens.assembly.ensembl109             # FASTA files of the genome, each chromosome in a separate file
homo_sapiens.annotation.ensembl109           # Annotation in GTF and TAB format
homo_sapiens.assembly.ensembl109.star        # STAR index, GTF annotation aware
homo_sapiens.transcripts.ensembl109          # transcriptome, this is the Ensembl "cDNA" file in FASTA format
homo_sapiens.transcripts.ensembl109.salmon   # Salmon index of the transcriptome
```

#### Retrieving genomic sequences

To retrieve stretches of genomic sequence, we use the seq(genome, chr, strand, position, upstream, downstream) method:

```
pybio.genomes.seq("hg38", "1", "+", 450000, -20, 20) # returns 'TACCCTGATTCTGAAACGAAAAAGCTTTACAAAATCCAAGA' for hg38, Ensembl v98
```

The above command fetches the chr1 sequence from 450000-20..450000+20, the resulting sequence is of length 41.

#### Annotating genomic positions

Given a genomic position, we can quickly retrieve the gene, transcript, exon and utr5/3 information at the given position. If there are several features (genes, transcripts, exons, UTR regions) at the specified position, they are all returned.

```
# annotate position
genes, transcripts, exons, utr5, utr3 = pybio.genomes.annotate("hg38", "1", "+", 11012344)

# print all genes that cover the position
for gene in genes:
   print(gene.gene_id, gene.gene_name)
   print(gene.start, gene.stop)
```

The above command would return:

```
[pybio] loading genome annotation for homo_sapiens with Ensembl version 109
ENSG00000120948, TARDBP
11012343, 11030527
```

You can also easily access all transcripts of each gene with `gene.transcripts` and all exons of each transcript with `transcript.exons`.

#### File formats

##### bedGraph

```
b = pybio.data.bedGraph()
```

### Authors

[pybio](https://github.com/grexor/pybio) is developed and supported by [Gregor Rot](https://grexor.github.io).

### Reporting problems

Use the [issues page](https://github.com/grexor/pybio/issues) to report issues and leave suggestions.
