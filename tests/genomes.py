import pybio

# AC
assert(pybio.genomes.seq("hg19", "1", "+", 100000, 100002)=="ACT")

