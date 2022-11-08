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
rm -r mm10.assembly.ensembl${eversion}
mkdir mm10.assembly.ensembl${eversion}
cd mm10.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -O mm10.fasta.gz
gunzip -f mm10.fasta.gz
printf 'import pybio\npybio.data.Fasta("mm10.fasta").split()' | python

cd $gdir
mkdir mm10.annotation.ensembl${eversion}
cd mm10.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/mm10.biomart.ensembl${eversion}.xml`
wget -O mm10.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("mm10", version="ensembl'${eversion}'")' | python

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/mus_musculus/Mus_musculus.GRCm39.${eversion}.chr.gtf.gz -O Mus_musculus.GRCm39.${eversion}.chr.gtf.gz
gunzip Mus_musculus.GRCm39.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r mm10.assembly.ensembl${eversion}.star
mkdir mm10.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir mm10.assembly.ensembl${eversion}.star --genomeFastaFiles mm10.assembly.ensembl${eversion}/mm10.fasta --runThreadN 4 --sjdbGTFfile mm10.annotation.ensembl${eversion}/Mus_musculus.GRCm39.${eversion}.chr.gtf

gzip mm10.annotation.ensembl${eversion}/Mus_musculus.GRCm39.${eversion}.chr.gtf # to save space

rm -r mm10.transcripts.ensembl${eversion}
mkdir mm10.transcripts.ensembl${eversion}
cd mm10.transcripts.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

cd $gdir
salmon index -t mm10.transcripts.ensembl${eversion}/Mus_musculus.GRCm39.cdna.all.fa.gz -i mm10.transcripts.ensembl${eversion}.salmon