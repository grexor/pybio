<picture><img src="media/pybio.png" height="30"/></picture>
### pybio: basic genomics toolset

*pybio* is a comprehensive Python framework designed to streamline genomics operations. It offers a direct interface to Ensembl genome assemblies and annotations, while also accommodating custom genomes via FASTA/GTF inputs. The primary objective of *pybio* is to simplify genome management. It achieves this by providing automatic download of Ensembl genome assemblies and annotation, provides Python genomic feature search and sequence retrieval from the managed genomes, STAR indexing and mapping and more.

### Quick Start

Install via pip and download + prepare human genome:

```
pip install pybio
pybio genome homo_sapiens
```

Search genome features (exons, transcripts, genes) from Python:

```
import pybio
result = pybio.core.genomes.annotate("homo_sapiens", "1", "+", 11012344)
genes, transcripts, exons, UTR5, UTR3 = result
```

Retrieve genomic sequences from Python:

```
import pybio
seq = pybio.core.genomes.seq("homo_sapiens", "1", "+", 450000, -20, 20)
```

Check documentation for more examples.

### Documentation

* [PDF reference manual](https://github.com/grexor/pybio/raw/master/docs/pybio_docs.pdf)
* [Google docs](https://docs.google.com/document/d/12KJvdsl78ujXaE3vTdGBK4vDgRRpHHh3RJg9npSVlZ4/edit?usp=sharing) of the above PDF (comment if you like)

### Authors

[pybio](https://github.com/grexor/pybio) is developed and supported by [Gregor Rot](https://grexor.github.io).

### Issues and Suggestions

Use the [issues page](https://github.com/grexor/pybio/issues) to report issues and leave suggestions.

### Change log

**v0.3.12**: November 2023
* updated docs