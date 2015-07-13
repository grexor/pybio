rm -r xenTro2.assembly
mkdir xenTro2.assembly
cd xenTro2.assembly
wget ftp://ftp.ensembl.org/pub/release-59/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI4.1.59.dna_rm.toplevel.fa.gz -O xenTro2.fasta.gz
gunzip -f xenTro2.fasta.gz
printf 'import pybio\npybio.data.fastasplit("xenTro2.fasta")' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" xenTro2 > xenTro2.chr.ucsc.ensembl
cd ..

rm -r xenTro2.assembly.star
mkdir xenTro2.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir xenTro2.assembly.star --genomeFastaFiles xenTro2.assembly/xenTro2.fasta --runThreadN 4

mkdir xenTro2.annotation.ensembl59
cd xenTro2.annotation.ensembl59
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../xenTro2.biomart.xml`
wget -O xenTro2.annotation.ensembl59.tab "http://aug2010.archive.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("xenTro2", version="ensembl59")' | python
