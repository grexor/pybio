# Genomic coordinates

All genomic coordinates inside pybio are **0-based, left+right inclusive**. So a range written `100-103` includes coordinates 100, 101, 102 and 103. The first coordinate of any sequence is 0.

!!! important
    RefSeq and Ensembl GTF files are 1-based. When pybio reads annotation from RefSeq or Ensembl, it subtracts 1 from every coordinate so that everything inside pybio — genes, transcripts, exons, UTRs, and the coordinates you pass to `seq()` and `annotate()` — stays consistently 0-based.

This matters whenever you compare a position from pybio against one from a 1-based source (a GTF file, a VCF, a genome browser coordinate you typed by hand): subtract 1 from the 1-based position before passing it to pybio, or add 1 back when displaying a pybio position to a user.
