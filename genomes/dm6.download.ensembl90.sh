rm -r dm6.assembly.ensembl90
mkdir dm6.assembly.ensembl90
cd dm6.assembly.ensembl90
wget ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz -O dm6.fasta.gz
gunzip -f dm6.fasta.gz
printf 'import pybio\npybio.data.Fasta("dm6.fasta").split()' | python
# there is no data for dm6 chr name mappings
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" dm6 > dm6.chr.ucsc.ensembl
touch dm6.chr.ucsc.ensembl
cd ..

mkdir dm6.annotation.ensembl90
cd dm6.annotation.ensembl90
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../dm6.biomart.ensembl90.xml`
wget -O dm6.annotation.ensembl90.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dm6")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-90/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.90.chr.gtf.gz -O Drosophila_melanogaster.BDGP6.90.chr.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.90.chr.gtf.gz # file must be unzipped for STAR to consider it
cd ..

rm -r dm6.assembly.ensembl90.star
mkdir dm6.assembly.ensembl90.star
STAR --runMode genomeGenerate --genomeDir dm6.assembly.ensembl90.star --genomeFastaFiles dm6.assembly.ensembl90/dm6.fasta --runThreadN 8 --sjdbGTFfile dm6.annotation.ensembl90/Drosophila_melanogaster.BDGP6.90.chr.gtf

gzip dm6.annotation.ensembl90/Drosophila_melanogaster.BDGP6.90.chr.gtf
