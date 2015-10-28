rm -r mm10.assembly
mkdir mm10.assembly
cd mm10.assembly
wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -O mm10.fasta.gz
gunzip -f mm10.fasta.gz
printf 'import pybio\npybio.data.Fasta("mm10.fasta").split()' | python
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" mm10 > mm10.chr.ucsc.ensembl
cd ..

rm -r mm10.assembly.star
mkdir mm10.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir mm10.assembly.star --genomeFastaFiles mm10.assembly/mm10.fasta --runThreadN 40

mkdir mm10.annotation.ensembl82
cd mm10.annotation.ensembl82
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../mm10.biomart.xml`
wget -O mm10.annotation.ensembl82.tab "http://www.biomart.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("mm10")' | python
