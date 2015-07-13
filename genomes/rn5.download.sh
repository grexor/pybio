rm -r rn5.assembly
mkdir rn5.assembly
cd rn5.assembly
wget ftp://ftp.ensembl.org/pub/release-74/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_5.0.74.dna_rm.toplevel.fa.gz -O rn5.fasta.gz
gunzip -f rn5.fasta.gz
printf 'import pybio\npybio.data.fastasplit("rn5.fasta")' | python
cd ..

rm -r rn5.assembly.star
mkdir rn5.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir rn5.assembly.star --genomeFastaFiles rn5.assembly/rn5.fasta --runThreadN 4

mkdir rn5.annotation.ensembl74
cd rn5.annotation.ensembl74
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../rn5.biomart.xml`
wget -O rn5.annotation.ensembl74.tab "http://www.biomart.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("rn5")' | python
