eversion="custom"

rm -r mar5.assembly.${eversion}
mkdir mar5.assembly.${eversion}
cd mar5.assembly.${eversion}
wget https://marchantia.info/download/tak1v5.1/MpTak1_with_female-chromosome.zip -O mar5.fasta.zip
unzip mar5.fasta.zip
mv tak1v5.1_with_female_chr.fasta mar5.fasta
printf 'import pybio\npybio.data.Fasta("mar5.fasta").split()' | python
cd ..

mkdir mar5.annotation.${eversion}
cd mar5.annotation.${eversion}
mv ../mar5.assembly.${eversion}/tak1v5.1r1_with_female_chr.gtf mar5.gtf
wget https://raw.githubusercontent.com/AleGirFon/Mytest/master/DB/Mapoly_2_MPGs.dic.txt -O mar5_mar3_gene_names.tab
#printf 'import pybio\npybio.genomes.prepare("mar3", version="ensembl'${eversion}'")' | python

cd ..
rm -r mar5.assembly.${eversion}.star
mkdir mar5.assembly.${eversion}.star
STAR --runMode genomeGenerate --genomeDir mar5.assembly.${eversion}.star --genomeFastaFiles mar5.assembly.${eversion}/mar5.fasta --genomeSAindexNbases 12 --runThreadN 4 --sjdbGTFfile mar5.annotation.${eversion}/mar5.gtf
