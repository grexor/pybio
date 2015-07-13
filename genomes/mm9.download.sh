rm -r mm9.assembly
mkdir mm9.assembly
cd mm9.assembly
wget ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz -O mm9.fasta.gz
gunzip -f mm9.fasta.gz
printf 'import pybio\npybio.data.fastasplit("mm9.fasta")' | python
cd ..

rm -r mm9.assembly.star
mkdir mm9.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir mm9.assembly.star --genomeFastaFiles mm9.assembly/mm9.fasta --runThreadN 4

mkdir mm9.annotation.ensembl67
cd mm9.annotation.ensembl67
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../mm9.biomart.xml`
wget -O mm9.annotation.ensembl67.tab "http://may2012.archive.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("panTro2", version="ensembl67")' | python
