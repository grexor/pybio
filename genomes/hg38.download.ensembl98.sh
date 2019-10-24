eversion="98"

rm -r hg38.assembly.ensembl${eversion}
mkdir hg38.assembly.ensembl${eversion}
cd hg38.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O hg38.fasta.gz
gunzip -f hg38.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg38.fasta").split()' | python
cd ..

mkdir hg38.annotation.ensembl${eversion}
cd hg38.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../hg38.biomart.ensembl${eversion}.xml`
wget -O hg38.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
# check this
printf 'import pybio\npybio.genomes.prepare("hg38", version="ensembl'${eversion}'")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${eversion}.chr.gtf.gz -O Homo_sapiens.GRCh38.${eversion}.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd ..
rm -r hg38.assembly.ensembl${eversion}.star
mkdir hg38.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir hg38.assembly.ensembl${eversion}.star --genomeFastaFiles hg38.assembly.ensembl${eversion}/hg38.fasta --runThreadN 4 --sjdbGTFfile hg38.annotation.ensembl${eversion}/Homo_sapiens.GRCh38.${eversion}.chr.gtf

gzip hg38.annotation.ensembl${eversion}/Homo_sapiens.GRCh38.${eversion}.chr.gtf # to save space

rm -r hg38.transcripts.ensembl${eversion}
mkdir hg38.transcripts.ensembl${eversion}
cd hg38.transcripts.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
cd ..
salmon index -t hg38.transcripts.ensembl${eversion}/Homo_sapiens.GRCh38.cdna.all.fa.gz -i hg38.transcripts.ensembl${eversion}.salmon
