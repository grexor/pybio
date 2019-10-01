rm -r tt.assembly.custom
mkdir tt.assembly.custom
cd tt.assembly.custom
wget http://www.ciliate.org/system/downloads/T_thermophila_June2014_assembly.fasta -O tt.fasta
printf 'import pybio\npybio.data.Fasta("tt.fasta").split()' | python
cd ..

mkdir tt.annotation.custom
cd tt.annotation.custom
wget http://www.ciliate.org/system/downloads/T_thermophila_June2014.gff3 -O tt.gff3

printf 'import pybio\npybio.genomes.prepare("tt")' | python
