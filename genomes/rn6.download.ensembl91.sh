rm -r rn6.assembly.ensembl91
mkdir rn6.assembly.ensembl91
cd rn6.assembly.ensembl91
wget ftp://ftp.ensembl.org/pub/release-91/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz -O rn6.fasta.gz
gunzip -f rn6.fasta.gz
printf 'import pybio\npybio.data.Fasta("rn6.fasta").split()' | python
cd ..

mkdir rn6.annotation.ensembl91
cd rn6.annotation.ensembl91
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../rn6.biomart.ensembl91.xml`
wget -O rn6.annotation.ensembl91.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("rn6", version="ensembl91")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-91/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.91.chr.gtf.gz -O Rattus_norvegicus.Rnor_6.0.91.chr.gtf.gz
gunzip Rattus_norvegicus.Rnor_6.0.91.chr.gtf.gz # file must be unzipped for STAR to consider it

cd ..
rm -r rn6.assembly.ensembl91.star
mkdir rn6.assembly.ensembl91.star
STAR --runMode genomeGenerate --genomeDir rn6.assembly.ensembl91.star --genomeFastaFiles rn6.assembly.ensembl91/rn6.fasta --runThreadN 8 --sjdbGTFfile rn6.annotation.ensembl91/Rattus_norvegicus.Rnor_6.0.91.chr.gtf
gzip Rattus_norvegicus.Rnor_6.0.91.chr.gtf # to save space
