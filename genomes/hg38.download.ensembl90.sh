#rm -r hg38.assembly
#mkdir hg38.assembly
#cd hg38.assembly
#wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O hg38.fasta.gz
#gunzip -f hg38.fasta.gz
#printf 'import pybio\npybio.data.Fasta("hg38.fasta").split()' | python
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg38 > hg38.chr.ucsc.ensembl
#cd ..

rm -r hg38.assembly.star
mkdir hg38.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir hg38.assembly.star --genomeFastaFiles hg38.assembly/hg38.fasta --runThreadN 4

#mkdir hg19.annotation.ensembl74
#cd hg19.annotation.ensembl74
# download annotation
#export BM=`sed ':a;N;$!ba;s/\n/ /g' ../hg19.biomart.xml`
#wget -O hg19.annotation.ensembl74.tab "http://dec2013.archive.ensembl.org/biomart/martservice?query=$BM"
#printf 'import pybio\npybio.genomes.prepare("hg19")' | python
