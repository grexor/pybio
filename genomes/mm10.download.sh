rm -r mm10.assembly
mkdir mm10.assembly
cd mm10.assembly
wget ftp://ftp.ensembl.org/pub/release-74/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.primary_assembly.fa.gz -O mm10.fasta.gz
gunzip -f mm10.fasta.gz
printf 'import pybio\npybio.data.fastasplit("mm10.fasta")' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" mm10 > mm10.chr.ucsc.ensembl
cd ..

rm -r mm10.assembly.star
mkdir mm10.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir mm10.assembly.star --genomeFastaFiles mm10.assembly/mm10.fasta --runThreadN 4

mkdir mm10.annotation.ensembl74
cd mm10.annotation.ensembl74
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../mm10.biomart.xml`
wget -O mm10.annotation.ensembl74.tab "http://www.biomart.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("mm10")' | python
