eversion="91"
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
rm -r rn6.assembly.ensembl${eversion}
mkdir rn6.assembly.ensembl${eversion}
cd rn6.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz -O rn6.fasta.gz
gunzip -f rn6.fasta.gz
printf 'import pybio\npybio.data.Fasta("rn6.fasta").split()' | python

cd $gdir
mkdir rn6.annotation.ensembl${eversion}
cd rn6.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/rn6.biomart.ensembl${eversion}.xml`
wget -O rn6.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("rn6", version="ensembl${eversion}")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.${eversion}.chr.gtf.gz -O Rattus_norvegicus.Rnor_6.0.${eversion}.chr.gtf.gz
gunzip Rattus_norvegicus.Rnor_6.0.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r rn6.assembly.ensembl${eversion}.star
mkdir rn6.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir rn6.assembly.ensembl${eversion}.star --genomeFastaFiles rn6.assembly.ensembl${eversion}/rn6.fasta --runThreadN 8 --sjdbGTFfile rn6.annotation.ensembl${eversion}/Rattus_norvegicus.Rnor_6.0.${eversion}.chr.gtf
gzip Rattus_norvegicus.Rnor_6.0.${eversion}.chr.gtf # to save space
