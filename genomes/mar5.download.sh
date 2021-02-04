eversion="custom"
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
rm -r mar5.assembly.${eversion}
mkdir mar5.assembly.${eversion}
cd mar5.assembly.${eversion}
wget https://marchantia.info/download/tak1v5.1/MpTak1_with_female-chromosome.zip -O mar5.fasta.zip
unzip mar5.fasta.zip
mv tak1v5.1_with_female_chr.fasta mar5.fasta
printf 'import pybio\npybio.data.Fasta("mar5.fasta").split()' | python

cd $gdir
mkdir mar5.annotation.${eversion}
cd mar5.annotation.${eversion}
mv ../mar5.assembly.${eversion}/tak1v5.1r1_with_female_chr.gtf mar5.gtf
wget https://raw.githubusercontent.com/AleGirFon/Mytest/master/DB/Mapoly_2_MPGs.dic.txt -O mar5_mar3_gene_names.tab

cd $gdir
rm -r mar5.assembly.${eversion}.star
mkdir mar5.assembly.${eversion}.star
STAR --runMode genomeGenerate --genomeDir mar5.assembly.${eversion}.star --genomeFastaFiles mar5.assembly.${eversion}/mar5.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile mar5.annotation.${eversion}/mar5.gtf

rm -r mar5.transcripts.${eversion}
mkdir mar5.transcripts.${eversion}
cd mar5.transcripts.${eversion}
wget https://marchantia.info/download/tak1v5.1/MpTak1v5.1_r1.mrna.fasta
gzip MpTak1v5.1_r1.mrna.fasta

cd $gdir
salmon index -t mar5.transcripts.${eversion}/MpTak1v5.1_r1.mrna.fasta.gz -i mar5.transcripts.${eversion}.salmon
