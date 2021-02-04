eversion="39"
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
rm -r at.assembly.ensembl${eversion}
mkdir at.assembly.ensembl${eversion}
cd at.assembly.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${eversion}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O at.fasta.gz
gunzip -f at.fasta.gz
printf 'import pybio\npybio.data.Fasta("at.fasta").split()' | python3

cd $gdir
mkdir at.annotation.ensembl${eversion}
cd at.annotation.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${eversion}/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${eversion}.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.${eversion}.gtf.gz
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/at.biomart.ensembl${eversion}.xml`
wget -O at.annotation.ensembl${eversion}.tab "http://plants.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("at")' | python3

cd $gdir
rm -r at.assembly.ensembl${eversion}.star
mkdir at.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir at.assembly.ensembl${eversion}.star --genomeFastaFiles at.assembly.ensembl${eversion}/at.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile at.annotation.ensembl${eversion}/Arabidopsis_thaliana.TAIR10.${eversion}.gtf
gzip at.annotation.ensembl${eversion}/Arabidopsis_thaliana.TAIR10.${eversion}.gtf

rm -r at.transcripts.ensembl${eversion}
mkdir at.transcripts.ensembl${eversion}
cd at.transcripts.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${eversion}/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

cd $gdir
salmon index -t at.transcripts.ensembl${eversion}/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -i at.transcripts.ensembl${eversion}.salmon
