rm -r hg19.assembly
mkdir hg19.assembly
cd hg19.assembly
wget ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.74.dna.primary_assembly.fa.gz -O hg19.fasta.gz
gunzip -f hg19.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg19.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg19 > hg19.chr.ucsc.ensembl
cd ..

rm -r hg19.assembly.star
mkdir hg19.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir hg19.assembly.star --genomeFastaFiles hg19.assembly/hg19.fasta --runThreadN 4

mkdir hg19.annotation.ensembl74
cd hg19.annotation.ensembl74
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../hg19.biomart.xml`
wget -O hg19.annotation.ensembl74.tab "http://dec2013.archive.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg19")' | python
