"""
pybio
====

A compact python module developed to handle common bioinformatics file formats, especially
in the next-generation sequencing (NGS) field.
"""

import pybio.path
import pybio.data
import pybio.map
import pybio.utils
import pybio.expression
import pybio.genomes
import pybio.maths
import pybio.config
import pybio.sequence
#import pybio.core
import os

# initialize path module
pybio.config.init()
pybio.path.init()
pybio.genomes.init()

def gff3_from_fasta(fasta_file):
    f = pybio.data.Fasta(fasta_file)
    while f.read():
        row = [f.id, "pybio", "chromosome", "1", len(f.sequence), ".", ".", ".", "ID=%s;Name=%s" % (f.id, f.id)]
        print("\t".join(str(x) for x in row))
