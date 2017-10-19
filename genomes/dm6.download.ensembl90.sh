rm -r dm6.assembly.ensembl90
mkdir dm6.assembly.ensembl90
cd dm6.assembly
wget ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz -O dm6.fasta.gz
gunzip -f dm6.fasta.gz
printf 'import pybio\npybio.data.Fasta("dm6.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" dm6 > dm6.chr.ucsc.ensembl
cd ..

rm -r dm6.assembly.star
mkdir dm6.assembly.star
STAR --runMode genomeGenerate --genomeDir dm6.assembly.ensembl90.star --genomeFastaFiles dm6.assembly.ensembl90/dm6.fasta --runThreadN 4

mkdir dm6.annotation.ensembl90
cd dm6.annotation.ensembl90
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../dm5.biomart.xml`
wget -O dm6.annotation.ensembl90.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dm6")' | python
