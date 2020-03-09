#!/usr/bin/env sh

# module load bowtie2/2.2.9-foss-2017a
# module load samtools/1.9-foss-2017a

idx=$1 ; shift
indir=$1 ; shift
outdir=$1 ; shift
pref=$1 ; shift
paired=$1 ; shift
threads=$1 ; shift
params=( "$@" )

if [ ! -d "$outdir" ]; then
	mkdir -p $outdir
fi

if [ $paired == 'false' ]; then
	# single end
	bowtie2 -x $idx  -U ${indir}/${pref}.fastq.gz -p $threads $(expr $params / 2) | samtools view -Sb | samtools sort -n -@ $(expr $threads / 2) -O BAM -o ${outdir}/${pref}.bam
else
	# paired end
	bowtie2 -x $idx  -1 ${indir}/${pref}_1.fastq.gz -2 ${indir}/${pref}_2.fastq.gz -p $(expr $threads / 2) $params | samtools view -Sb | samtools sort -n -@ $(expr $threads / 2) -O BAM -o ${outdir}/${pref}.bam
fi

samtools sort -@ $threads -o ${outdir}/${pref}_sorted.bam ${outdir}/${pref}.bam
samtools index -b ${outdir}/${pref}_sorted.bam
