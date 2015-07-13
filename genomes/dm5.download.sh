rm -r dm5.assembly
mkdir dm5.assembly
cd dm5.assembly
wget ftp://ftp.ensembl.org/pub/release-74/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.74.dna_rm.toplevel.fa.gz -O dm5.fasta.gz
gunzip -f dm5.fasta.gz
printf 'import pybio\npybio.data.fastasplit("dm5.fasta")' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" dm5 > dm5.chr.ucsc.ensembl
cd ..

rm -r dm5.assembly.star
mkdir dm5.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir dm5.assembly.star --genomeFastaFiles dm5.assembly/dm5.fasta --runThreadN 4

mkdir dm5.annotation.ensembl74
cd dm5.annotation.ensembl74
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../dm5.biomart.xml`
wget -O dm5.annotation.ensembl74.tab "http://www.biomart.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dm5")' | python
