<picture><img src="media/pybio.png" height="30"/></picture>
### pybio: basic genomics toolset

*pybio* is a comprehensive Python framework designed to streamline genomics operations. It offers a direct interface to Ensembl genome assemblies and annotations, while also accommodating custom genomes via FASTA/GTF inputs. The primary objective of *pybio* is to simplify genome management. It achieves this by providing automatic download of Ensembl genome assemblies and annotation, provides Python genomic feature search and sequence retrieval from the managed genomes, STAR indexing and mapping and more.

### Quick Start

Install and download + prepare human genome:

```
# Option 1: install over PyPi
pip install pybio

# Option 2: install from this repository
pip install git+https://github.com/grexor/pybio.git@master

# Option 3: use over singularity / apptainer / Docker (only if you don't need python imports)
singularity run docker://ghcr.io/grexor/pybio:master pybio

# Download and process homo sapiens genome
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

**0.7**: February 2025
* alignIntronMax for STAR
* other small fixes

**v0.6.3**: December 2024
* updated setup.py to use an entry point instead of a script
* removed `pybio` scripts

**v0.6**: November 2024
* updated Ensembl search and genome versioning offline
* updated custom genome interface

<details>
<summary><b>v0.5</b>: May 2024</summary>

* refreshed Ensembl (112) and Ensembl Genomes (58) database
</details>

<details>
<summary><b>v0.4</b>: April 2024</summary>

* refreshed Ensembl (111) and Ensembl Genomes (58) database
</details>

<details>
<summary><b>v0.3.12</b>: released in November 2023</summary>

* updated docs
</details>

### Citation

If you are using pybio in your research, please cite:

Rot, G., Wehling, A., Schmucki, R., Berntenis, N., Zhang, J. D., & Ebeling, M. (2024)<br>
[splicekit : an integrative toolkit for splicing analysis from short-read RNA-seq](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae121/7735317)<br>
Bioinformatics Advances, 4(1). https://doi.org/10.1093/bioadv/vbae121
