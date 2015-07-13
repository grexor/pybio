rm -r zv9.assembly
mkdir zv9.assembly
cd zv9.assembly
wget ftp://ftp.ensembl.org/pub/release-74/fasta/danio_rerio/dna/Danio_rerio.Zv9.74.dna_rm.toplevel.fa.gz -O zv9.fasta.gz
gunzip -f zv9.fasta.gz
printf 'import pybio\npybio.data.fastasplit("zv9.fasta")' | python
cd ..

rm -r zv9.assembly.star
mkdir zv9.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir zv9.assembly.star --genomeFastaFiles zv9.assembly/zv9.fasta --runThreadN 4

mkdir zv9.annotation.ensembl74
cd zv9.annotation.ensembl74
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../zv9.biomart.xml`
wget -O zv9.annotation.ensembl74.tab "http://www.biomart.org/biomart/martservice?query=$BM"
