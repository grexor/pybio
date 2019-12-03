"""
pybio3
====

A compact python module developed to handle common bioinformatics file formats, especially
in the next-generation sequencing (NGS) field.
"""

import pybio3.path
import pybio3.data
import pybio3.map
import pybio3.utils
import pybio3.expression
import pybio3.genomes
import pybio3.maths
import pybio3.config
import pybio3.sequence
import os

# initialize path module
pybio3.path.init()
pybio3.genomes.init()

def gff3_from_fasta(fasta_file):
    f = pybio3.data.Fasta(fasta_file)
    while f.read():
        row = [f.id, "pybio", "chromosome", "1", len(f.sequence), ".", ".", ".", "ID=%s;Name=%s" % (f.id, f.id)]
        print("\t".join(str(x) for x in row))
