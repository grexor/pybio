rm -r hg38.assembly.ensembl96
mkdir hg38.assembly.ensembl96
cd hg38.assembly.ensembl96
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O hg38.fasta.gz
gunzip -f hg38.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg38.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg38 > hg38.chr.ucsc.ensembl
cd ..

mkdir hg38.annotation.ensembl96
cd hg38.annotation.ensembl96
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../hg38.biomart.ensembl96.xml`
wget -O hg38.annotation.ensembl96.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg38", version="ensembl96")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gtf.gz -O Homo_sapiens.GRCh38.96.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.96.chr.gtf.gz # file must be unzipped for STAR to consider it

cd ..
rm -r hg38.assembly.ensembl96.star
mkdir hg38.assembly.ensembl96.star
STAR --runMode genomeGenerate --genomeDir hg38.assembly.ensembl96.star --genomeFastaFiles hg38.assembly.ensembl96/hg38.fasta --runThreadN 4 --sjdbGTFfile hg38.annotation.ensembl96/Homo_sapiens.GRCh38.96.chr.gtf

gzip hg38.annotation.ensembl96/Homo_sapiens.GRCh38.96.chr.gtf # to save space
