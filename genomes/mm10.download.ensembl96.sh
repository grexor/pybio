rm -r mm10.assembly.ensembl96
mkdir mm10.assembly.ensembl96
cd mm10.assembly.ensembl96
wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -O mm10.fasta.gz
gunzip -f mm10.fasta.gz
printf 'import pybio\npybio.data.Fasta("mm10.fasta").split()' | python
cd ..

mkdir mm10.annotation.ensembl96
cd mm10.annotation.ensembl96
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../mm10.biomart.ensembl96.xml`
wget -O mm10.annotation.ensembl96.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("mm10", version="ensembl96")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.chr.gtf.gz -O Mus_musculus.GRCm38.96.chr.gtf.gz
gunzip Mus_musculus.GRCm38.96.chr.gtf.gz # file must be unzipped for STAR to consider it

cd ..
rm -r mm10.assembly.ensembl96.star
mkdir mm10.assembly.ensembl96.star
STAR --runMode genomeGenerate --genomeDir mm10.assembly.ensembl96.star --genomeFastaFiles mm10.assembly.ensembl96/mm10.fasta --runThreadN 8 --sjdbGTFfile mm10.annotation.ensembl96/Mus_musculus.GRCm38.96.chr.gtf

gzip mm10.annotation.ensembl96/Mus_musculus.GRCm38.96.chr.gtf # to save space
