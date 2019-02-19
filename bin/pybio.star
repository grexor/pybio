#!/bin/sh

output_folder=$1
genome_dir=$2
fastq=$3
name=$4
cpu=$5
minlen=$6

mkdir $output_folder
cd $output_folder

# https://groups.google.com/forum/#!topic/rna-star/7RwKkvNLmI4
# --genomeLoad LoadAndRemove
echo "STAR --readFilesCommand bzip2 -c --outFilterMultimapNmax 1 --genomeDir $2 --readFilesIn $3 --outReadsUnmapped Fastx --runThreadN ${cpu} --outFilterMatchNminOverLread ${minlen} --outFilterScoreMinOverLread ${minlen}"
STAR --outFilterMultimapNmax 1 --genomeDir $2 --readFilesIn $3 --outReadsUnmapped Fastx --runThreadN ${cpu} --outFilterMatchNminOverLread ${minlen} --outFilterScoreMinOverLread ${minlen} --readFilesCommand bunzip2 -c

mv Aligned.out.sam ${name}.sam
mv Log.final.out ${name}.stats.txt
mv Log.out ${name}.log.txt
mv Unmapped.out.mate1 ${name}.unmapped.fastq
mv Log.progress.out ${name}.progress.txt
pybio.sam2bam ${name}.sam
rm ${name}.sam
gzip -f ${name}.unmapped.fastq
mv SJ.out.tab ${name}.splice.tab.out
