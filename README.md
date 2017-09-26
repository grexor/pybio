# pybio: basic genomics toolset

pybio is a Python framework for basic genomics operations. The package is a dependency to [apa](https://github.com/grexor/apa) (alternative polyadenylation) and [RNAmotifs2](https://github.com/grexor/rnamotifs2) (motif cluster analysis). Historically pybio was a side project with the intention of servicing other project requirements. Nevertheless, it can be used independently and provides:

+ download of genome assemblies from Ensembl together with automated index building (STAR short-read aligner),
+ download of genome annotations from Ensembl GTF with fast-search indexing,
+ Fasta, Fastq, bedGraph and other file format handling,
+ motif sequence searches,
+ alternative polyadenylation site-pair classification (same-exon, skipped-exon, composite-exon),
+ and other.

## Authors

[pybio](https://github.com/grexor/pybio) is developed and supported by [Gregor Rot](http://rotlab.info).

## Reporting problems

Use the [issues page](https://github.com/grexor/pybio/issues) to report issues and leave suggestions.
