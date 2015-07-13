rm -r panTro2.assembly
mkdir panTro2.assembly
cd panTro2.assembly
wget ftp://ftp.ensembl.org/pub/release-64/fasta/pan_troglodytes/dna/Pan_troglodytes.CHIMP2.1.64.dna_rm.toplevel.fa.gz -O panTro2.fasta.gz
gunzip -f panTro2.fasta.gz
printf 'import pybio\npybio.data.fastasplit("panTro2.fasta")' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" panTro2 > panTro2.chr.ucsc.ensembl
cd ..

rm -r panTro2.assembly.star
mkdir panTro2.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir panTro2.assembly.star --genomeFastaFiles panTro2.assembly/panTro2.fasta --runThreadN 4

mkdir panTro2.annotation.ensembl64
cd panTro2.annotation.ensembl64
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../panTro2.biomart.xml`
wget -O panTro2.annotation.ensembl64.tab "http://sep2011.archive.ensembl.org//biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("panTro2", version="ensembl64")' | python
