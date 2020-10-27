rm -r mar5.assembly
mkdir mar5.assembly
cd mar5.assembly
wget https://marchantia.info/download/tak1v5.1/MpTak1_with_female-chromosome.zip -O mar5.fasta.zip
unzip mar5.fasta.zip
mv tak1v5.1_with_female_chr.fasta mar5.fasta
printf 'import pybio\npybio.data.Fasta("mar5.fasta").split()' | python
cd ..

mkdir mar5.annotation
cd mar5.annotation
mv ../mar5.assembly/tak1v5.1r1_with_female_chr.gtf mar5.gtf
#printf 'import pybio\npybio.genomes.prepare("mar3", version="ensembl'${eversion}'")' | python
cd ..

rm -r mar5.assembly.star
mkdir mar5.assembly.star
STAR --runMode genomeGenerate --genomeDir mar5.assembly.star --genomeFastaFiles mar5.assembly/mar5.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile mar5.annotation/mar5.gtf
cd ..

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
