import pybio
import sys
import os

pybio_folder = os.path.dirname(pybio.__file__)
genomes_folder = os.path.join(pybio_folder, "genomes", "genomes_data")
test_annotation_folder = os.path.join(genomes_folder, "hg38.annotation.test")

os.system(f"cp *.gtf {test_annotation_folder}")
pybio.genomes.prepare_gtf("hg38", "test")

for i in range(0, 1000+1):
    print(i, pybio.genomes.annotate("hg38", 1, "+", i, version="test"))
