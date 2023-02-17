#!/usr/bin/env bash
species=${1}
eversion=${2}

declare -A ensembl_longA=( ["hg38"]="homo_sapiens")
declare -A ensembl_longB=( ["hg38"]="Homo_sapiens")
declare -A ensembl_longC=( ["hg38"]="GRCh38")

sdir=`realpath .`
gdir=${3}
echo $gdir
if [ -z "$gdir" ]
then
  gdir=`realpath .`
else
  mkdir -p $gdir
  gdir=`realpath ${3}`
fi

cd $gdir
rm -r ${species}.assembly.ensembl${eversion}
mkdir ${species}.assembly.ensembl${eversion}
cd ${species}.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/${ensembl_longA[hg38]}/dna/${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.dna.primary_assembly.fa.gz -O ${species}.fasta.gz
echo "[pybio] unzipping ${species}.fasta.gz"
gunzip -f ${species}.fasta.gz
python3 -c "import pybio; pybio.data.Fasta('${species}.fasta').split()"

cd $gdir
mkdir ${species}.annotation.ensembl${eversion}
cd ${species}.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/ensembl_download.xml`
wget -O ${species}.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
python3 -c "import pybio; pybio.genomes.prepare('${species}', version='ensembl${eversion}')"

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/${ensembl_longA[hg38]}/${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.${eversion}.chr.gtf.gz -O ${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.${eversion}.chr.gtf.gz
gunzip -f ${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r ${species}.assembly.ensembl${eversion}.star
mkdir ${species}.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir ${species}.assembly.ensembl${eversion}.star --genomeFastaFiles ${species}.assembly.ensembl${eversion}/${species}.fasta --runThreadN 4 --sjdbGTFfile ${species}.annotation.ensembl${eversion}/${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.${eversion}.chr.gtf

gzip -f ${species}.annotation.ensembl${eversion}/${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.${eversion}.chr.gtf # to save space

rm -r ${species}.transcripts.ensembl${eversion}
mkdir ${species}.transcripts.ensembl${eversion}
cd ${species}.transcripts.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/${ensembl_longA[hg38]}/cdna/${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.cdna.all.fa.gz

cd $gdir
salmon index -t ${species}.transcripts.ensembl${eversion}/${ensembl_longB[hg38]}.${ensembl_longC[hg38]}.cdna.all.fa.gz -i ${species}.transcripts.ensembl${eversion}.salmon
