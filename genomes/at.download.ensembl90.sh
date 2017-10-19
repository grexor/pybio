rm -r at.assembly.ensembl90
mkdir at.assembly.ensembl90
cd at.assembly.ensembl90
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-37/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O at.fasta.gz
gunzip -f at.fasta.gz
printf 'import pybio\npybio.data.Fasta("at.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" tair10 > dm6.chr.ucsc.ensembl
cd ..

rm -r at.assembly.ensembl90.star
mkdir at.assembly.ensembl90.star
STAR --runMode genomeGenerate --genomeDir at.assembly.ensembl90.star --genomeFastaFiles at.assembly.ensembl90/at.fasta --runThreadN 4

mkdir at.annotation.ensembl90
cd at.annotation.ensembl90
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../at.biomart.ensembl90.xml`
wget -O at.annotation.ensembl90.tab "http://plants.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("at")' | python
