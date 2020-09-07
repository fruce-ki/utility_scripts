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
    echo "      $0 -i BAM -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR -l RLEN [-x XREF]"
    exit 1
}
# Parse options.
minq=27
readlen=100
pre=1
dunk=1
post=1
while getopts 'i:r:b:a:o:q:l:x:123' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# folder of bams
    r) ref="${OPTARG}" ;;			# whole reference fasta for alignment
    b) bed="${OPTARG}" ;;			# bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;			# GTF annotation (for counting)
    o) outdir="${OPTARG}" ;;		# output directory
    q) minq="${OPTARG}" ;;			# minimum alignement quality (10)
    l) readlen="${OPTARG}" ;;		# maximum read length
    x) xref="${OPTARG}" ;;		    # "Name" to other IDs
    1) pre=0 ;;         		    # Skip bam 2 FQ and trim
    2) dunk=0 ;;        		    # Skip slamdunk
    3) post=0 ;;		            # Skip alleyoop and mutliqc
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

if [ "$pre" -eq 1 ]; then
    echo "$bam - BAM to FASTQ"
    fileutilities.py T $bam/*bam --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bam2fq ,--get-user-env ,--wrap "'samtools bam2fq -0 ${outdir}/fastq/{bas}.fastq {abs}'"
    wait_for_jobs bam2fq
    fileutilities.py T ${outdir}/fastq/*.fastq --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqzip ,--get-user-env ,--wrap "'gzip -f {abs}'"
    wait_for_jobs fqzip

    echo "$bam - QC"
    fileutilities.py T ${outdir}/fastq/*.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fastQC ,--get-user-env ,--wrap "'fastqc -o ${outdir}/fastqc_pre {abs}'"
    # wait_for_jobs fastQC

    echo "$bam - trim, filter, QC"
    fileutilities.py T ${outdir}/fastq/*.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,-J trimgalo ,--get-user-env trim_galore ,-o ${outdir}/fastq_trimmed ,--stringency 3 ,--gzip ,--fastqc ,--fastqc_args "',--outdir ${outdir}/fastqc_post/'" {abs}
    wait_for_jobs trimgalo
    wait_for_jobs fastQC

    echo "$bam - MultiQC"
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_pre ${outdir}/fastqc_pre
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_post ${outdir}/fastqc_post ${outdir}/fastq_trimmed
fi

fileutilities.py T ${outdir}/fastq_trimmed/ --dir 'fastq.gz$|fq.gz$' | perl -e 'while(<>){~s/\.fastq$\.fq$/\tdummy\t0/; print}' > ${outdir}/dunk/samplesheet.tsv
threads=$( wc -l ${outdir}/dunk/samplesheet.tsv | perl -e 'while(<>){~/^(\d+)/; print $1}' )
if [ "$threads" -gt 24 ]; then
    threads=24
fi
# memory=$((4 * $threads))
memory='40G'
echo "$bam - Threads: $threads Memory: $memory"

if [ "$dunk" -eq 1 ]; then
    echo "$bam - slamdunk"
    sbatch --qos=medium -c $threads --mem=$memory -J slamdunk -o /dev/null -e /dev/null --wrap "slamdunk all -r $ref -b $bed -o ${outdir}/dunk -5 12 -n 100 -c 2 -mv 0.2 -t $threads -rl $readlen -mbq $minq -m -ss ${outdir}/dunk/samplesheet.tsv"
    wait_for_jobs slamdunk
fi

if [ "$post" -eq 1 ]; then
    echo "$bam - alleyoop"
    mkdir -p ${outdir}/alleyoop/summary/ ${outdir}/alleyoop/rates ${outdir}/alleyoop/utrrates ${outdir}/alleyoop/tcperreadpos ${outdir}/alleyoop/tcperutrpos

    sbatch -J alleyoop -o /dev/null -e /dev/null --wrap "alleyoop summary -t ${outdir}/dunk/count/ -o ${outdir}/alleyoop/summary/summary.txt ${outdir}/dunk/filter/*bam"
    # source ~/miniconda3/bin/activate slamdunk && <...>
    fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-o /dev/null ,-e /dev/null ,-c $threads  ,-J alleyoop ,--wrap "'alleyoop rates -o ${outdir}/alleyoop/rates -r $ref -t $threads  -mq $minq {abs}'"
    fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-o /dev/null ,-e /dev/null ,-c $threads  ,-J alleyoop ,--wrap "'alleyoop utrrates -o ${outdir}/alleyoop/utrrates -r $ref -t $threads  -b $bed -l $readlen {abs}'"
    fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-o /dev/null ,-e /dev/null ,-c $threads  ,-J alleyoop ,--wrap "'alleyoop tcperreadpos -o ${outdir}/alleyoop/tcperreadpos -r $ref -s ${outdir}/dunk/snp -l $readlen -mq $minq -t $threads  {abs}'"
    fileutilities.py T ${outdir}/dunk/filter/*bam --loop sbatch ,-o /dev/null ,-e /dev/null ,-c $threads  ,-J alleyoop ,--wrap "'alleyoop tcperutrpos -o ${outdir}/alleyoop/tcperutrpos -r $ref -b $bed -s ${outdir}/dunk/snp -l $readlen -mq $minq -t $threads  ${outdir}/dunk/filter/*bam'"
    wait_for_jobs alleyoop
    # echo "# slamdunk summary v0.3.4" > ${outdir}/alleyoop/summary.txt
    # echo "FileName	SampleName	SampleType	SampleTime	Sequenced	Mapped	Deduplicated	MQ-Filtered	Identity-Filtered	NM-Filtered	Multimap-Filtered	Retained	Counted	Annotation" >> ${outdir}/alleyoop/summary.txt
    # cat ${outdir}/alleyoop/summary/*summary.txt | tail -n 1 >> ${outdir}/alleyoop/summary.txt

    echo "$bam - MultiQC"
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_alleyoop ${outdir}/fastqc_post ${outdir}/alleyoop/summary ${outdir}/alleyoop/rates ${outdir}/alleyoop/tcperreadpos ${outdir}/alleyoop/tcperutrpos ${outdir}/alleyoop/utrrates

    echo "$bam - RPMu"
    srun slamseq_rpmu.R ${outdir}/dunk/count/*tcount.tsv

    echo "$bam - Merge"
    fileutilities.py T ${outdir}/dunk/count/*rpmu.txt -i -r --appnd > ${outdir}/all_rpmu.txt

    dups=$(head -n1 ${outdir}/all_rpmu.txt | fileutilities.py D --swap "\n" | perl -e '$i=0; while($field = <STDIN>){print "$i " if $field=~/Name/; $i++} print "\n";')
    srun --mem=10G --ntasks=1 fileutilities.py T ${outdir}/all_rpmu.txt -r --mrgdups $dups > ${outdir}/all_rpmu_dedup.txt
    nc=$(perl -e '$ARGV[0] =~/^(\d+)/; print $1' $(fileutilities.py T ${outdir}/all_rpmu_dedup.txt --cntcols))
    srun --mem=10G --ntasks=1 fileutilities.py T ${outdir}/all_rpmu_dedup.txt -r --cols $(expr $nc - 1) 0 1:$(expr $nc - 2) > ${outdir}/all_rpmu_dedup_reord.txt
    rm ${outdir}/all_rpmu_dedup.txt
    mv ${outdir}/all_rpmu_dedup_reord.txt ${outdir}/all_rpmu.txt

    if [ -e "${xref}" ] ; then
      # Unlike quantseq, the output counts of slamdunk are per UTR, not per gene. So Entrez are not unique and cannot be used as keys in pandas/fileutilities.py.
      echo "$bam - Adding crossreferencing IDs"
      slamseq_xref.R ${outdir} 'all_rpmu.txt' $xref 1 1
    fi
fi

# Clean up file pollution
# rm ./slurm*

echo "$bam - SLAMseq workflow finished!"
