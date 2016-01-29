rm -r wbcel235.assembly
mkdir wbcel235.assembly
cd wbcel235.assembly
wget ftp://ftp.ensembl.org/pub/release-77/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_rm.toplevel.fa.gz -O wbcel235.fasta.gz
gunzip -f wbcel235.fasta.gz
printf 'import pybio\npybio.data.Fasta("wbcel235.fasta").split()' | python
cd ..

rm -r wbcel235.assembly.star
mkdir wbcel235.assembly.star
# STAR index
STAR --runMode genomeGenerate --genomeDir wbcel235.assembly.star --genomeFastaFiles wbcel235.assembly/wbcel235.fasta --runThreadN 4

mkdir wbcel235.annotation.ensembl77
cd wbcel235.annotation.ensembl77
# download annotation
export BM=`sed ':a;N;$!ba;s/\n/ /g' ../wbcel235.biomart.xml`
wget -O wbcel235.annotation.ensembl77.tab "http://www.biomart.org/biomart/martservice?query=$BM"
printf 'import pybio\npybio.genomes.prepare("wbcel235", version="ensembl77")' | python
