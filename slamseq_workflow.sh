#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH -J slmSeq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20000
#SBATCH --output=slmSeq.out
#SBATCH --error=slmSeq.err
#SBATCH --qos=medium
#SBATCH --time=1-0:00:00       # one day

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i BAM -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR -x XREF"
    exit 1
}
# Parse options.
minq=27
readlen=100
while getopts 'i:r:b:a:o:q:l:x:' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# folder of bams
    r) ref="${OPTARG}" ;;			# whole reference fasta for alignment
    b) bed="${OPTARG}" ;;			# bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;			# GTF annotation (for counting)
    o) outdir="${OPTARG}" ;;		# output directory
    q) minq="${OPTARG}" ;;			# minimum alignement quality (10)
    l) readlen="${OPTARG}" ;;		# maximum read length
    x) xref="${OPTARG}" ;;		    # "Name" to other IDs
    *) usage ;;
  esac
done

wait_for_jobs(){
  echo "waiting"
  sleep 60  # seconds, give time to the scheduler to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep kimon.fr | grep -c $1) -eq 0 ]; then
      break
    else
      echo "sleep another" $((sleeptime / 60)) "minutes..."
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}

mkdir -p ${outdir}/fastqc_pre ${outdir}/fastqc_post ${outdir}/fastq ${outdir}/fastq_trimmed ${outdir}/dunk ${outdir}/alleyoop ${outdir}/multiqc_pre ${outdir}/multiqc_post ${outdir}/multiqc_alleyoop

# echo "$bam - BAM to FASTQ"
# fileutilities.py T $bam/*bam --loop sbatch ,-J bam2fq ,--get-user-env ,--wrap "'samtools bam2fq -0 ${outdir}/fastq/{bas}.fastq {abs}'"
# wait_for_jobs bam2fq
# fileutilities.py T ${outdir}/fastq/*.fastq --loop sbatch ,-J fqzip ,--get-user-env ,--wrap "'gzip -f {abs}'"
# wait_for_jobs fqzip
#
# echo "$bam - QC"
# fileutilities.py T ${outdir}/fastq/*.fastq.gz --loop sbatch ,-J fastQC ,--get-user-env ,--wrap "'fastqc -o ${outdir}/fastqc_pre {abs}'"
# wait_for_jobs fastQC
#
# echo "$bam - trim, filter, QC"
# fileutilities.py T ${outdir}/fastq/*.fastq.gz --loop sbatch ,--ntasks=1 ,-c 4 ,-J trimgalo ,--get-user-env trim_galore ,-j 4 ,-o ${outdir}/fastq_trimmed ,--stringency 3 ,--gzip ,--fastqc ,--fastqc_args "',--outdir ${outdir}/fastqc_post/'" {abs}
# wait_for_jobs trimgalo


# # source ~/miniconda3/bin/deactivate
# # source ~/miniconda3/bin/activate slamdunk


echo "$bam - slamdunk"
fileutilities.py T ${outdir}/fastq_trimmed/*_trimmed.fq.gz --loop sbatch ,-c 6 ,--mem=24G ,-J slamdunk ,--wrap "'source ~/miniconda3/bin/activate slamdunk && slamdunk all -r $ref -b $bed -o ${outdir}/dunk -5 12 -n 100 -c 2 -mv 0.2 -t 4 -rl $readlen -mbq $minq -m -ss {abs}'"
wait_for_jobs slamdunk

echo "$bam - alleyoop"
mkdir -p ${outdir}/alleyoop/summary/ ${outdir}/alleyoop/rates ${outdir}/alleyoop/utrrates ${outdir}/alleyoop/tcperreadpos ${outdir}/alleyoop/tcperutrpos

sbatch --qos short --time=0-0:10:00 --ntasks=1 -J alleyoop --wrap "source ~/miniconda3/bin/activate slamdunk && alleyoop summary -t ${outdir}/dunk/count/ -o ${outdir}/alleyoop/summary/summary.txt ${outdir}/dunk/filter/*bam"
fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-c 4 ,-J alleyoop ,--wrap "'source ~/miniconda3/bin/activate slamdunk && alleyoop rates -o ${outdir}/alleyoop/rates -r $ref -t 4 -mq $minq {abs}'"
fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-c 4 ,-J alleyoop ,--wrap "'source ~/miniconda3/bin/activate slamdunk && alleyoop utrrates -o ${outdir}/alleyoop/utrrates -r $ref -t 4 -b $bed -l $readlen {abs}'"
fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-c 4 ,-J alleyoop ,--wrap "'source ~/miniconda3/bin/activate slamdunk && alleyoop tcperreadpos -o ${outdir}/alleyoop/tcperreadpos -r $ref -s ${outdir}/dunk/snp -l $readlen -mq $minq -t 4 {abs}'"
fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-c 4 ,-J alleyoop ,--wrap "'source ~/miniconda3/bin/activate slamdunk && alleyoop tcperutrpos -o ${outdir}/alleyoop/tcperutrpos -r $ref -b $bed -s ${outdir}/dunk/snp -l $readlen -mq $minq -t 4 ${outdir}/dunk/filter/*bam'"
wait_for_jobs alleyoop
# # echo "# slamdunk summary v0.3.4" > ${outdir}/alleyoop/summary.txt
# # echo "FileName	SampleName	SampleType	SampleTime	Sequenced	Mapped	Deduplicated	MQ-Filtered	Identity-Filtered	NM-Filtered	Multimap-Filtered	Retained	Counted	Annotation" >> ${outdir}/alleyoop/summary.txt
# # cat ${outdir}/alleyoop/summary/*summary.txt | tail -n 1 >> ${outdir}/alleyoop/summary.txt

echo "$bam - MultiQC"
srun --qos short --ntasks=1 --mpi=none multiqc -f -o ${outdir}/multiqc_pre ${outdir}/fastqc_pre
srun --qos short --ntasks=1 --mpi=none multiqc -f -o ${outdir}/multiqc_post ${outdir}/fastqc_post
srun --qos short --ntasks=1 --mpi=none multiqc -f -o ${outdir}/multiqc_alleyoop ${outdir}/alleyoop/summary ${outdir}/alleyoop/rates ${outdir}/alleyoop/tcperreadpos ${outdir}/alleyoop/tcperutrpos ${outdir}/alleyoop/utrrates

if [ -e "${xref}" ] ; then
    echo "$bam - Adding crossreferencing IDs"
    xref.R ${outdir}/dunk/count $xref 4 1 '.tsv'
fi



# # Prepare file headers for pasting
# # mkdir ${outdir}/dunk/count_renamed
# # fileutilities.py T ${outdir}/dunk/count --dir tcount.tsv | perl -e 'while(<>){~s/(\t\d+).*$/$1/;print}' |  ???

# # Extract and collate the counts, lose the commented command.
# # fileutilities.py T ${outdir}/dunk/count/*txt -i -l --cols 6 | perl -e 'while(<>){~s/_\|6//g;print}' > ${outdir}/counts/all_counts.tsv


# Clean up file pollution
# rm ./slurm*
