output_folder=$1
genome_dir=$2
fastq=$3
name=$4
cpu=$5

mkdir $output_folder
cd $output_folder

STAR --outFilterMultimapNmax 1 --genomeDir $2 --readFilesIn $3 --outReadsUnmapped Fastx --readFilesCommand zcat

# https://groups.google.com/forum/#!topic/rna-star/7RwKkvNLmI4
# STAR --outFilterMultimapNmax 1 --genomeDir ${genome_dir} --readFilesIn ${fastq} --outReadsUnmapped Fastx --readFilesCommand zcat --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --runThreadN ${cpu}
mv Aligned.out.sam ${name}.sam
mv Log.final.out ${name}.stats.txt
mv Log.out ${name}.log.txt
mv Unmapped.out.mate1 ${name}.unmapped.fastq
mv Log.progress.out ${name}.progress.txt
pybio.sam2bam ${name}.sam
rm ${name}.sam
gzip ${name}.unmapped.fastq
mv SJ.out.tab ${name}.splice.tab.out
