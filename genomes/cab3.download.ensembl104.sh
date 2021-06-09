eversion="104"
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
rm -r cab3.assembly.ensembl${eversion}
mkdir cab3.assembly.ensembl${eversion}
cd cab3.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/equus_caballus/dna/Equus_caballus.EquCab3.0.dna.toplevel.fa.gz -O cab3.fasta.gz
gunzip -f cab3.fasta.gz
printf 'import pybio\npybio.data.Fasta("cab3.fasta").split()' | python3

cd $gdir
mkdir cab3.annotation.ensembl${eversion}
cd cab3.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/cab3.biomart.ensembl${eversion}.xml`
wget -O cab3.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("cab3", version="ensembl'${eversion}'")' | python3

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/equus_caballus/Equus_caballus.EquCab3.0.${eversion}.gtf.gz -O cab3.${eversion}.chr.gtf.gz
gunzip -f cab3.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r cab3.assembly.ensembl${eversion}.star
mkdir cab3.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir cab3.assembly.ensembl${eversion}.star --genomeFastaFiles cab3.assembly.ensembl${eversion}/cab3.fasta --runThreadN 4 --sjdbGTFfile cab3.annotation.ensembl${eversion}/cab3.${eversion}.chr.gtf

gzip -f cab3.annotation.ensembl${eversion}/cab3.${eversion}.chr.gtf # to save space

rm -r cab3.transcripts.ensembl${eversion}
mkdir cab3.transcripts.ensembl${eversion}
cd cab3.transcripts.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/equus_caballus/cdna/Equus_caballus.EquCab3.0.cdna.all.fa.gz

cd $gdir
salmon index -t cab3.transcripts.ensembl${eversion}/Equus_caballus.EquCab3.0.cdna.all.fa.gz -i cab3.transcripts.ensembl${eversion}.salmon
