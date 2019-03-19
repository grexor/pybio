# http://marchantia.info/download/

rm -r mar3.assembly
mkdir mar3.assembly
cd mar3.assembly
wget http://marchantia.info/download/download/JGI_3.1.fasta.gz -O mar3.fasta.gz
gunzip -f mar3.fasta.gz
printf 'import pybio\npybio.data.Fasta("mar3.fasta").split()' | python
cd ..

rm -r mar3.assembly.star
mkdir mar3.assembly.star
STAR --runMode genomeGenerate --genomeDir mar3.assembly.star --genomeFastaFiles mar3.assembly/mar3.fasta --runThreadN 4

mkdir mar3.annotation
cd mar3.annotation
wget http://marchantia.info/download/download/Mpolymorphav3.1.allTrs.gff3.gz -O mar3.gff.gz
#export BM=`sed ':a;N;$!ba;s/\n/ /g' ../at.biomart.ensembl39.xml`
#wget -O at.annotation.ensembl39.tab "http://plants.ensembl.org/biomart/martservice?query=$BM"
#printf 'import pybio\npybio.genomes.prepare("at")' | python
