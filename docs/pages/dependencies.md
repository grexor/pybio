# Dependencies

[pysam](https://pysam.readthedocs.io/en/latest/api.html), [numpy](https://numpy.org/) and [samtools](http://www.htslib.org/) are required and installed automatically by `pip install pybio`.

[STAR](https://github.com/alexdobin/STAR) and [salmon](https://combine-lab.github.io/salmon/getting_started/) are optional — install them yourself (and make sure they're on your `PATH`) if you want pybio to build genome/transcriptome indices and align reads. See [Read mapping](mapping.md).

[bedGraphToBigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) (from UCSC's genome browser tools) is only needed if you convert bedGraph tracks to BigWig with `pybio.data.bedgraph_bigwig()` — see [CLIP & interval data](bedgraph.md#converting-to-bigwig).
