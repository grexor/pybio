eversion="75"
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
rm -r hg19.assembly.ensembl${eversion}
mkdir hg19.assembly.ensembl${eversion}
cd hg19.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.${eversion}.dna.primary_assembly.fa.gz -O hg19.fasta.gz
gunzip -f hg19.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg19.fasta").split()' | python3
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg19 > hg19.chr.ucsc.ensembl

cd $gdir
mkdir hg19.annotation.ensembl${eversion}
cd hg19.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/hg19.biomart.ensembl${eversion}.xml`
wget -O hg19.annotation.ensembl${eversion}.tab "http://feb2014.archive.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg19", version="ensembl${eversion}")' | python3

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/homo_sapiens/Homo_sapiens.GRCh37.${eversion}.gtf.gz -O Homo_sapiens.GRCh37.${eversion}.gtf.gz
gunzip Homo_sapiens.GRCh37.${eversion}.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r hg19.assembly.ensembl${eversion}.star
mkdir hg19.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir hg19.assembly.ensembl${eversion}.star --genomeFastaFiles hg19.assembly.ensembl${eversion}/hg19.fasta --runThreadN 8 --sjdbGTFfile hg19.annotation.ensembl${eversion}/Homo_sapiens.GRCh37.${eversion}.gtf
gzip hg19.annotation.ensembl${eversion}/Homo_sapiens.GRCh37.${eversion}.gtf # to save space
