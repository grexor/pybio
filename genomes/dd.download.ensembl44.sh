eversion="44"
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
rm -r dd.assembly.ensembl${eversion}
mkdir dd.assembly.ensembl${eversion}
cd dd.assembly.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-${eversion}/fasta/dictyostelium_discoideum/dna/Dictyostelium_discoideum.dicty_2.7.dna.toplevel.fa.gz -O dd.fasta.gz
gunzip -f dd.fasta.gz
printf 'import pybio\npybio.data.Fasta("dd.fasta").split()' | python3

exit

cd $gdir
mkdir dd.annotation.ensembl${eversion}
cd dd.annotation.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-${eversion}/gtf/dictyostelium_discoideum/Dictyostelium_discoideum.dicty_2.7.44.gtf.gz
gunzip Dictyostelium_discoideum.dicty_2.7.44.gtf.gz
export BM=`sed ':a;N;$!ba;s/\n/ /g' $sdir/dd.biomart.ensembl${eversion}.xml`
wget -O dd.annotation.ensembl${eversion}.tab "http://protists.ensembl.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("dd", version="ensembl'${eversion}'")' | python3

cd $gdir
rm -r dd.assembly.ensembl${eversion}.star
mkdir dd.assembly.ensembl${eversion}.star
STAR --runMode genomeGenerate --genomeDir dd.assembly.ensembl${eversion}.star --genomeFastaFiles dd.assembly.ensembl${eversion}/dd.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile dd.annotation.ensembl${eversion}/Dictyostelium_discoideum.dicty_2.7.44.gtf
gzip dd.annotation.ensembl${eversion}/Dictyostelium_discoideum.dicty_2.7.44.gtf

rm -r dd.transcripts.ensembl${eversion}
mkdir dd.transcripts.ensembl${eversion}
cd dd.transcripts.ensembl${eversion}
wget ftp://ftp.ensemblgenomes.org/pub/protists/release-44/fasta/dictyostelium_discoideum/cdna/Dictyostelium_discoideum.dicty_2.7.cdna.all.fa.gz
cd ..
salmon index -t dd.transcripts.ensembl${eversion}/Dictyostelium_discoideum.dicty_2.7.cdna.all.fa.gz -i dd.transcripts.ensembl${eversion}.salmon
