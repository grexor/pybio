rm -r tt.assembly.custom
mkdir tt.assembly.custom
cd tt.assembly.custom
wget http://www.ciliate.org/system/downloads/T_thermophila_June2014_assembly.fasta -O tt.fasta
printf 'import pybio\npybio.data.Fasta("tt.fasta").split()' | python
cd ..

mkdir tt.annotation.custom
cd tt.annotation.custom
wget http://www.ciliate.org/system/downloads/T_thermophila_June2014.gff3 -O tt.gff

cd ..

python tt_gff_gtf.py

gzip -f tt.annotation.custom/tt.gff

printf 'import pybio\npybio.genomes.prepare("tt")' | python

rm -r tt.assembly.custom.star
mkdir tt.assembly.custom.star
STAR --runMode genomeGenerate --genomeDir tt.assembly.custom.star --genomeFastaFiles tt.assembly.custom/tt.fasta --genomeSAindexNbases 12 --runThreadN 4
