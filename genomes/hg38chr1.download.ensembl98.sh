eversion="98"
sdir=`realpath .`
gdir=$1
if [ -z "$gdir" ]
then
  gdir=`realpath .`
else
  mkdir -p $gdir
  gdir=`realpath $1`
fi

cd $gdir
rm -r hg38chr1.assembly.ensembl${eversion}
mkdir hg38chr1.assembly.ensembl${eversion}
cd hg38chr1.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz -O hg38chr1.fasta.gz
gunzip -f hg38chr1.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg38chr1.fasta").split()' | python3

cd $gdir
mkdir hg38chr1.annotation.ensembl${eversion}
cd hg38chr1.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/hg38chr1.biomart.ensembl${eversion}.xml`
wget -O hg38chr1.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg38chr1", version="ensembl'${eversion}'")' | python3

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${eversion}.chr.gtf.gz -O Homo_sapiens.GRCh38.${eversion}.chr.gtf.gz
gunzip -f Homo_sapiens.GRCh38.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r hg38chr1.assembly.ensembl${eversion}.star
mkdir hg38chr1.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir hg38chr1.assembly.ensembl${eversion}.star --genomeFastaFiles hg38chr1.assembly.ensembl${eversion}/hg38chr1.fasta --runThreadN 4 --sjdbGTFfile hg38chr1.annotation.ensembl${eversion}/Homo_sapiens.GRCh38.${eversion}.chr.gtf

gzip -f hg38chr1.annotation.ensembl${eversion}/Homo_sapiens.GRCh38.${eversion}.chr.gtf # to save space

rm -r hg38chr1.transcripts.ensembl${eversion}
mkdir hg38chr1.transcripts.ensembl${eversion}
cd hg38chr1.transcripts.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

cd $gdir
salmon index -t hg38chr1.transcripts.ensembl${eversion}/Homo_sapiens.GRCh38.cdna.all.fa.gz -i hg38chr1.transcripts.ensembl${eversion}.salmon
