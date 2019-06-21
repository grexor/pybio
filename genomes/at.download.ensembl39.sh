rm -r at.assembly.ensembl39
mkdir at.assembly.ensembl39
cd at.assembly.ensembl39
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O at.fasta.gz
gunzip -f at.fasta.gz
printf 'import pybio\npybio.data.Fasta("at.fasta").split()' | python
cd ..

mkdir at.annotation.ensembl39
cd at.annotation.ensembl39
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.39.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.39.gtf.gz
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../at.biomart.ensembl39.xml`
wget -O at.annotation.ensembl39.tab "http://plants.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("at")' | python

cd ..
rm -r at.assembly.ensembl39.star
mkdir at.assembly.ensembl39.star
STAR --runMode genomeGenerate --genomeDir at.assembly.ensembl39.star --genomeFastaFiles at.assembly.ensembl39/at.fasta --runThreadN 4 --sjdbGTFfile at.annotation.ensembl39/Arabidopsis_thaliana.TAIR10.39.gtf --alignIntronMin 60 --alignIntronMax 6000
gzip at.annotation.ensembl39/Arabidopsis_thaliana.TAIR10.39.gtf
