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
rm -r hg38chr22.assembly.ensembl${eversion}
mkdir hg38chr22.assembly.ensembl${eversion}
cd hg38chr22.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz -O hg38chr22.fasta.gz
gunzip -f hg38chr22.fasta.gz
printf 'import pybio\npybio.data.Fasta("hg38chr22.fasta").split()' | python3

cd $gdir
mkdir hg38chr22.annotation.ensembl${eversion}
cd hg38chr22.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/hg38chr22.biomart.ensembl${eversion}.xml`
wget -O hg38chr22.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("hg38chr22", version="ensembl'${eversion}'")' | python3

cd $gdir
rm -r hg38chr22.assembly.ensembl${eversion}.star
mkdir hg38chr22.assembly.ensembl${eversion}.star
STAR --genomeSAindexNbases 11 --runMode genomeGenerate --genomeDir hg38chr22.assembly.ensembl${eversion}.star --genomeFastaFiles hg38chr22.assembly.ensembl${eversion}/hg38chr22.fasta --runThreadN 4
