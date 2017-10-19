rm -r hg19.assembly.ensembl75
mkdir hg19.assembly.ensembl75
cd hg19.assembly.ensembl75
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz -O hg19.fasta.gz
gunzip -f hg19.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg19.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg19 > hg19.chr.ucsc.ensembl
cd ..

rm -r hg19.assembly.ensembl75.star
mkdir hg19.assembly.ensembl75.star
STAR --runMode genomeGenerate --genomeDir hg19.assembly.ensembl75.star --genomeFastaFiles hg19.assembly.ensembl75/hg19.fasta --runThreadN 4

mkdir hg19.annotation.ensembl75
cd hg19.annotation.ensembl75
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../hg19.biomart.ensembl75.xml`
wget -O hg19.annotation.ensembl75.tab "http://feb2014.archive.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg19", version="ensembl75")' | python
