rm -r hg38.assembly.ensembl90
mkdir hg38.assembly.ensembl90
cd hg38.assembly.ensembl90
wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O hg38.fasta.gz
gunzip -f hg38.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg38.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg38 > hg38.chr.ucsc.ensembl
cd ..

rm -r hg38.assembly.ensembl90.star
mkdir hg38.assembly.ensembl90.star
STAR --runMode genomeGenerate --genomeDir hg38.assembly.ensembl90.star --genomeFastaFiles hg38.assembly.ensembl90/hg38.fasta --runThreadN 4

mkdir hg38.annotation.ensembl90
cd hg38.annotation.ensembl90
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../hg38.biomart.ensembl90.xml`
wget -O hg38.annotation.ensembl90.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg38", version="ensembl90")' | python
