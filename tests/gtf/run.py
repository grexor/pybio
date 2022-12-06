import pybio
import sys
import os

pybio_folder = os.path.dirname(pybio.__file__)
genomes_folder = os.path.join(pybio_folder, "genomes")
genomes_data_folder = os.path.join(pybio_folder, "genomes", "genomes_data")
detected_genomes_data_folder = False
if os.path.exists(genomes_data_folder):
    detected_genomes_data_folder = True
if not detected_genomes_data_folder:
    genomes_data_folder = os.path.join(pybio_folder, "genomes")
    if os.path.exists(genomes_data_folder):
        detected_genomes_data_folder = True

if not detected_genomes_data_folder:
    print("can't detect the location of the genomes folder")
    sys.exit(1)

print("genomes data folder = ", genomes_data_folder)

# download hg38 assembly for testing
test_assembly_folder = os.path.join(genomes_data_folder, "hg38.assembly.test")
if not os.path.exists(test_assembly_folder):
    os.system(f"cd {genomes_folder}; ./hg38.download.test.sh {genomes_data_folder}")

# prepare pybio test annotation from provided test gtf
test_annotation_folder = os.path.join(genomes_data_folder, "hg38.annotation.test")
os.system(f"mkdir {test_annotation_folder}")
os.system(f"cp *.gtf {test_annotation_folder}")
pybio.genomes.prepare_gtf("hg38", "test")

for i in range(0, 1000+1):
    print(i, pybio.genomes.annotate("hg38", 1, "+", i, version="test"))
