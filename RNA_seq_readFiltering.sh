#!/usr/bin/bash

# Make sure that the shell script stops running if it hits an error
set -e

# Read sample file path if not already set in a shell parameter
 [ -z "$1" ] && read -p "Provide sample file path --> " sampleFilePath || sampleFilePath=$1

# Read sample name if not already set in a shell parameter
 [ -z "$2" ] && read -p "Provide sample name --> " sampleName || sampleName=$2

echo "Running $sampleName ($sampleFilePath)"

# Read trimming
fastp -i "$sampleFilePath"/"$sampleName"_1.fq.gz -I "$sampleFilePath"/"$sampleName"_2.fq.gz -o "$sampleFilePath"/"$sampleName"_1_fastp.fq.gz -O "$sampleFilePath"/"$sampleName"_2_fastp.fq.gz --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim_poly_x

fastqc "$sampleFilePath"/"$sampleName"_1_fastp.fq.gz
fastqc "$sampleFilePath"/"$sampleName"_2_fastp.fq.gz

rm "$sampleFilePath"/"$sampleName"_1.fq.gz
rm "$sampleFilePath"/"$sampleName"_2.fq.gz

# Read alignment
./STAR-2.7.11b/source/STAR --genomeDir ColCEN/STAR_index --readFilesIn "$sampleFilePath"/"$sampleName"_1_fastp.fq.gz "$sampleFilePath"/"$sampleName"_2_fastp.fq.gz --outSAMtype BAM Unsorted --runThreadN 8 --readFilesCommand zcat --outFileNamePrefix "$sampleFilePath"/"$sampleName" --outTmpDir ~/STAR_temp

samtools sort -n "$sampleFilePath"/"$sampleName"Aligned.out.bam -o "$sampleFilePath"/"$sampleName"_alignedSorted.out.bam

rm "$sampleFilePath"/"$sampleName"Aligned.out.bam

# Read counting
htseq-count -s no -f bam "$sampleFilePath"/"$sampleName"_alignedSorted.out.bam ColCEN/ColCEN_GENES.gtf > "$sampleFilePath"/"$sampleName"_readCounts.txt

echo "Done with $sampleName."