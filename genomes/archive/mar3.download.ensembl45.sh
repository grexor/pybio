eversion="47"

rm -r mar3.assembly.ensembl${eversion}
mkdir mar3.assembly.ensembl${eversion}
cd mar3.assembly.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${eversion}/fasta/marchantia_polymorpha/dna/Marchantia_polymorpha.Marchanta_polymorpha_v1.dna.toplevel.fa.gz -O mar3.fasta.gz
gunzip -f mar3.fasta.gz
printf 'import pybio\npybio.data.Fasta("mar3.fasta").split()' | python
cd ..

mkdir mar3.annotation.ensembl${eversion}
cd mar3.annotation.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${eversion}/gtf/marchantia_polymorpha/Marchantia_polymorpha.Marchanta_polymorpha_v1.${eversion}.gtf.gz
gunzip Marchantia_polymorpha.Marchanta_polymorpha_v1.${eversion}.gtf.gz
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../mar3.biomart.ensembl${eversion}.xml`
wget -O mar3.annotation.ensembl${eversion}.tab "http://plants.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("mar3", version="ensembl'${eversion}'")' | python

cd ..
rm -r mar3.assembly.ensembl${eversion}.star
mkdir mar3.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir mar3.assembly.ensembl${eversion}.star --genomeFastaFiles mar3.assembly.ensembl${eversion}/mar3.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile mar3.annotation.ensembl${eversion}/Marchantia_polymorpha.Marchanta_polymorpha_v1.${eversion}.gtf
gzip mar3.annotation.ensembl${eversion}/Marchantia_polymorpha.Marchanta_polymorpha_v1.${eversion}.gtf

rm -r mar3.transcripts.ensembl${eversion}
mkdir mar3.transcripts.ensembl${eversion}
cd mar3.transcripts.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${eversion}/fasta/marchantia_polymorpha/cdna/Marchantia_polymorpha.Marchanta_polymorpha_v1.cdna.all.fa.gz
cd ..
salmon index -t mar3.transcripts.ensembl${eversion}/Marchantia_polymorpha.Marchanta_polymorpha_v1.cdna.all.fa.gz -i mar3.transcripts.ensembl${eversion}.salmon

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
