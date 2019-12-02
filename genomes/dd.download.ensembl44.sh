eversion="44"

rm -r dd.assembly.ensembl${eversion}
mkdir dd.assembly.ensembl${eversion}
cd dd.assembly.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-${eversion}/fasta/dictyostelium_discoideum/dna/Dictyostelium_discoideum.dicty_2.7.dna.toplevel.fa.gz -O dd.fasta.gz
gunzip -f dd.fasta.gz
printf 'import pybio\npybio.data.Fasta("dd.fasta").split()' | python
cd ..

mkdir dd.annotation.ensembl${eversion}
cd dd.annotation.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-${eversion}/gtf/dictyostelium_discoideum/Dictyostelium_discoideum.dicty_2.7.44.gtf.gz
gunzip Dictyostelium_discoideum.dicty_2.7.44.gtf.gz
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../dd.biomart.ensembl${eversion}.xml`
wget -O dd.annotation.ensembl${eversion}.tab "http://protists.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dd", version="ensembl'${eversion}'")' | python

cd ..
rm -r dd.assembly.ensembl${eversion}.star
mkdir dd.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir dd.assembly.ensembl${eversion}.star --genomeFastaFiles dd.assembly.ensembl${eversion}/dd.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile dd.annotation.ensembl${eversion}/Dictyostelium_discoideum.dicty_2.7.44.gtf
gzip dd.annotation.ensembl${eversion}/Dictyostelium_discoideum.dicty_2.7.44.gtf

rm -r dd.transcripts.ensembl${eversion}
mkdir dd.transcripts.ensembl${eversion}
cd dd.transcripts.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-44/fasta/dictyostelium_discoideum/cdna/Dictyostelium_discoideum.dicty_2.7.cdna.all.fa.gz
cd ..
salmon index -t dd.transcripts.ensembl${eversion}/Dictyostelium_discoideum.dicty_2.7.cdna.all.fa.gz -i dd.transcripts.ensembl${eversion}.salmon

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
