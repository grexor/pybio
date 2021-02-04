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
rm -r dm6.assembly.ensembl${eversion}
mkdir dm6.assembly.ensembl${eversion}
cd dm6.assembly.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz -O dm6.fasta.gz
gunzip -f dm6.fasta.gz
printf 'import pybio\npybio.data.Fasta("dm6.fasta").split()' | python3
touch dm6.chr.ucsc.ensembl

cd $gdir
mkdir dm6.annotation.ensembl${eversion}
cd dm6.annotation.ensembl${eversion}
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/dm6.biomart.ensembl${eversion}.xml`
wget -O dm6.annotation.ensembl${eversion}.tab "http://www.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dm6", version="ensembl'${eversion}'")' | python3

# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-${eversion}/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.${eversion}.chr.gtf.gz -O Drosophila_melanogaster.BDGP6.${eversion}.chr.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.${eversion}.chr.gtf.gz # file must be unzipped for STAR to consider it

cd $gdir
rm -r dm6.assembly.ensembl${eversion}.star
mkdir dm6.assembly.ensembl${eversion}.star
# https://groups.google.com/forum/#!topic/rna-star/6csLuVjxiR0
# dm6 = 140M bases, log2(140M)/2-1 = 12.5, we set --genomeSAindexNbases to 10
STAR --runMode genomeGenerate --genomeSAindexNbases 10 --genomeDir dm6.assembly.ensembl${eversion}.star --genomeFastaFiles dm6.assembly.ensembl${eversion}/dm6.fasta --runThreadN 8 --sjdbGTFfile dm6.annotation.ensembl${eversion}/Drosophila_melanogaster.BDGP6.${eversion}.chr.gtf

gzip dm6.annotation.ensembl${eversion}/Drosophila_melanogaster.BDGP6.${eversion}.chr.gtf

rm -r dm6.transcripts.ensembl${eversion}
mkdir dm6.transcripts.ensembl${eversion}
cd dm6.transcripts.ensembl${eversion}
wget ftp://ftp.ensembl.org/pub/release-${eversion}/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz

cd $gdir
salmon index -t dm6.transcripts.ensembl${eversion}/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz -i dm6.transcripts.ensembl${eversion}.salmon
