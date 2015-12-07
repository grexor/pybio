mkdir $1
cd $1
STAR --outFilterMultimapNmax 1 --genomeDir $2 --readFilesIn $3 --outReadsUnmapped Fastx --readFilesCommand zcat
# https://groups.google.com/forum/#!topic/rna-star/7RwKkvNLmI4
# STAR --outFilterMultimapNmax 1 --genomeDir $2 --readFilesIn $3 --outReadsUnmapped Fastx --readFilesCommand zcat --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --runThreadN $5
mv Aligned.out.sam $4.sam
mv Log.final.out $4.stats.txt
mv Log.out $4.log.txt
mv Unmapped.out.mate1 $4.unmapped.fastq
mv Log.progress.out $4.progress.txt
pybio.sam2bam $4.sam
rm $4.sam
gzip $4.unmapped.fastq
mv SJ.out.tab $4.splice.tab.out
