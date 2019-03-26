#!/bin/sh

module load bowtie2/2.2.9-foss-2017a
module load samtools/1.9-foss-2017a

idx=$1 ; shift
dir=$1 ; shift
pref=$1 ; shift
paired=$1 ; shift
threads=$1 ; shift
params=( "$@" )


if [ $paired == 'false' ]; then
	# single end
	bowtie2 -x $idx  -U ${dir}/${pref}.fastq.gz -p $threads $params | samtools view -Sb | samtools sort -n -@ $threads -O BAM -o ${dir}/${pref}.bam
else
	# paired end
	bowtie2 -x $idx  -1 ${dir}/${pref}_1.fastq.gz -2 ${dir}/${pref}_2.fastq.gz -p $threads $params | samtools view -Sb | samtools sort -n -@ $threads -O BAM -o ${dir}/${pref}.bam
fi

samtools sort -@ $threads -o ${dir}/${pref}_sorted.bam ${dir}/${pref}.bam
samtools index -b ${dir}/${pref}_sorted.bam

