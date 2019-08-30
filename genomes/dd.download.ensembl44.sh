rm -r dd.assembly.ensembl44
mkdir dd.assembly.ensembl44
cd dd.assembly.ensembl44
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-44/fasta/dictyostelium_discoideum/dna/Dictyostelium_discoideum.dicty_2.7.dna.toplevel.fa.gz -O dd.fasta.gz
gunzip -f dd.fasta.gz
printf 'import pybio\npybio.data.Fasta("dd.fasta").split()' | python
cd ..

mkdir dd.annotation.ensembl44
cd dd.annotation.ensembl44
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-44/gtf/dictyostelium_discoideum/Dictyostelium_discoideum.dicty_2.7.44.gtf.gz
gunzip Dictyostelium_discoideum.dicty_2.7.44.gtf.gz
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../dd.biomart.ensembl44.xml`
wget -O dd.annotation.ensembl44.tab "http://protists.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dd")' | python

cd ..
rm -r dd.assembly.ensembl44.star
mkdir dd.assembly.ensembl44.star
STAR --runMode genomeGenerate --genomeDir dd.assembly.ensembl44.star --genomeFastaFiles dd.assembly.ensembl44/dd.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile dd.annotation.ensembl44/Dictyostelium_discoideum.dicty_2.7.44.gtf
gzip dd.annotation.ensembl44/Dictyostelium_discoideum.dicty_2.7.44.gtf


# another alternative, not Ensembl
# http://marchantia.info/download/

#rm -r mar3.assembly
#mkdir mar3.assembly
#cd mar3.assembly
#wget http://marchantia.info/download/download/JGI_3.1.fasta.gz -O mar3.fasta.gz
#gunzip -f mar3.fasta.gz
#printf 'import pybio\npybio.data.Fasta("mar3.fasta").split()' | python
#cd ..

#rm -r mar3.assembly.star
#mkdir mar3.assembly.star
#STAR --runMode genomeGenerate --genomeDir mar3.assembly.star --genomeFastaFiles mar3.assembly/mar3.fasta --runThreadN 4

#mkdir mar3.annotation
#cd mar3.annotation
#wget http://marchantia.info/download/download/Mpolymorphav3.1.allTrs.gff3.gz -O mar3.gff.gz
#export BM=`sed ':a;N;$!ba;s/\n/ /g' ../at.biomart.ensembl39.xml`
#wget -O at.annotation.ensembl39.tab "http://plants.ensembl.org/biomart/martservice?query=$BM"
#printf 'import pybio\npybio.genomes.prepare("at")' | python
