rm -r anoCar2.assembly
mkdir anoCar2.assembly
cd anoCar2.assembly
wget ftp://ftp.ensembl.org/pub/release-76/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna_rm.toplevel.fa.gz -O anoCar2.fasta.gz
gunzip -f anoCar2.fasta.gz
printf 'import pybio\npybio.data.Fasta("anoCar2.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" anoCar2 > anoCar2.chr.ucsc.ensembl
cd ..

rm -r anoCar2.assembly.star
mkdir anoCar2.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir anoCar2.assembly.star --genomeFastaFiles anoCar2.assembly/anoCar2.fasta --runThreadN 4

mkdir anoCar2.annotation.ensembl76
cd anoCar2.annotation.ensembl76
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../anoCar2.biomart.xml`
wget -O anoCar2.annotation.ensembl76.tab "http://www.biomart.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("anoCar2", version="ensembl76")' | python
