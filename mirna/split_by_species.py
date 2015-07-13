import pybio
import sys

f = pybio.data.Fasta("mirna.fasta")
while f.read():
    id = f.id
    species = id.split("-")[0]
    fout = open("mirna.%s.fasta" % species, "a")
    fout.write(">%s\n%s\n" % (f.id, f.sequence))
    fout.close()
