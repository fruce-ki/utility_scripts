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


## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i BAM -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR -l RLEN [-x XREF] [-k SPIKEIN_CHROM_LIST]"
    exit 1
}
# Parse options.
minq=27
readlen=100
pre=1
dunk=1
post=1
spikes=0
while getopts 'i:r:b:a:o:q:l:x:k123' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# folder of bams
    r) ref="${OPTARG}" ;;			# whole reference fasta for alignment
    b) bed="${OPTARG}" ;;			# bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;			# GTF annotation (for counting)
    o) outdir="${OPTARG}" ;;		# output directory
    q) minq="${OPTARG}" ;;			# minimum alignement quality (10)
    l) readlen="${OPTARG}" ;;		# maximum read length
    x) xref="${OPTARG}" ;;		    # xref table from "Name" to other IDs
    k) spikes=1 ;;		# Drosophila spike-in
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
    # ,-o /dev/null ,-e /dev/null
    fileutilities.py T $bam --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bam2fq ,--get-user-env ,--wrap "'samtools bam2fq -0 ${outdir}/fastq/{bas}.fastq {abs}'"
    wait_for_jobs bam2fq
    fileutilities.py T ${outdir}/fastq --dir 'fastq$|fq$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqzip ,--get-user-env ,--wrap "'gzip -f {abs}'"
    wait_for_jobs fqzip

    echo "$bam - QC"
    fileutilities.py T ${outdir}/fastq --dir 'fastq.gz$|fq.gz$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fastQC ,--get-user-env ,--wrap "'fastqc -o ${outdir}/fastqc_pre {abs}'"
    # wait_for_jobs fastQC

    echo "$bam - trim, filter, QC"
    fileutilities.py T ${outdir}/fastq --dir 'fastq.gz$|fq.gz$' | fileutilities.py P  --loop sbatch ,-o /dev/null ,-e /dev/null ,-J trimgalo ,--get-user-env trim_galore ,-o ${outdir}/fastq_trimmed ,--stringency 3 ,--gzip ,--fastqc ,--fastqc_args "',--outdir ${outdir}/fastqc_post/'" {abs}
    wait_for_jobs trimgalo
    wait_for_jobs fastQC

    echo "$bam - MultiQC"
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_pre ${outdir}/fastqc_pre
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_post ${outdir}/fastqc_post ${outdir}/fastq_trimmed
fi

fileutilities.py T ${outdir}/fastq_trimmed/ --dir 'fastq.gz$|fq.gz$' | perl -e 'while(<>){~s/\.fastq$|\.fq$/\tdummy\t0/; print}' > ${outdir}/dunk/samplesheet.tsv
threads=$( wc -l ${outdir}/dunk/samplesheet.tsv | perl -e 'while(<>){~/^(\d+)/; print $1}' )
if [ "$threads" -gt 24 ]; then
    threads=24
fi
# memory=$((4 * $threads))
memory='40G'
echo "$bam - Threads: $threads Memory: $memory"

if [ "$dunk" -eq 1 ]; then
    if ! [ -x "$(command -v slamdunk)" ]; then
        echo "Error: slamdunk could not be found. Did you load the module?" >&2
        exit 3
    fi

    echo "$bam - slamdunk"
    sbatch -o /dev/null -e /dev/null --qos=medium -c $threads --mem=$memory -J slamdunk --wrap "slamdunk all -r $ref -b $bed -o ${outdir}/dunk -5 12 -n 100 -c 2 -mv 0.2 -t $threads -rl $readlen -mbq $minq -m -ss ${outdir}/dunk/samplesheet.tsv"
    wait_for_jobs slamdunk

    if [ "$spikes" -eq 1 ]; then
        echo "$bam - separate spike-in from subject"
        # Restructure
        rm -r ${outdir}/dunk/filter_presplit
        mv ${outdir}/dunk/filter ${outdir}/dunk/filter_presplit
        mkdir -p ${outdir}/dunk/filter
        mkdir -p ${outdir}/dunk/filter_split
        rm ${outdir}/dunk/count/* ${outdir}/dunk/filter_split/*
        # De-spike
        fileutilities.py T ${outdir}/dunk/filter_presplit --dir 'bam$' | fileutilities.py P --loop sbatch ,--qos=medium ,-o /dev/null ,-e /dev/null ,--mem=$memory ,-J slamspik Split_BAM_Spikein.py ,-a {abs} ,-o ${outdir}/dunk/filter_split
        wait_for_jobs slamspik
        # Plug subject reads back into slamdunk
        fileutilities.py T ${outdir}/dunk/filter_split --dir subject | perl -e 'while(<>){~s/_subject$//;print}' | fileutilities.py P --link ${outdir}/dunk/filter/

        echo "$bam - requantify"
        # Quantify features.
        # Yes I could split apart `slamdunk all` into its constituent steps and insert the de-spiking in between and quantify only once on the correct files,
        # but that would add failure opportunitiess that I don't want to deal with at the moment.
        fileutilities.py T ${outdir}/dunk/filter --dir 'bam$' | fileutilities.py P --loop sbatch ,-J samidx ,-o /dev/null ,-e /dev/null ,--wrap "'samtools index {abs}'"
        wait_for_jobs samidx
        fileutilities.py T ${outdir}/dunk/filter --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J slamcnt slamdunk count ,-o ${outdir}/dunk/count ,-s ${outdir}/dunk/snp ,-r $ref ,-b $bed ,-l $readlen ,-q $minq ,-t 1 {abs}
        # Quantify library partitions. I can use samtools because slamdunk filters out multimappers. (therefore ambiguous should be empty)
        fileutilities.py T ${outdir}/dunk/filter_split --dir 'spikein.bam$' | fileutilities.py P --loop printf '"%s\t"' '"{cor}"' '>>' ${outdir}/dunk/filter_split/spike_counts.txt \&\& samtools view ,-F 4 ,-c {abs} '>>' ${outdir}/dunk/filter_split/spike_counts.txt
        fileutilities.py T ${outdir}/dunk/filter_split --dir 'ambiguous.bam$' | fileutilities.py P --loop printf '"%s\t"' '"{cor}"' '>>' ${outdir}/dunk/filter_split/ambiguous_counts.txt \&\& samtools view ,-F 4 ,-c {abs} '>>' ${outdir}/dunk/filter_split/ambiguous_counts.txt
        fileutilities.py T ${outdir}/dunk/filter_split --dir 'subject.bam$' | fileutilities.py P --loop printf '"%s\t"' '"{cor}"' '>>' ${outdir}/dunk/filter_split/subject_counts.txt \&\& samtools view ,-F 4 ,-c {abs} '>>' ${outdir}/dunk/filter_split/subject_counts.txt
        fileutilities.py T ${outdir}/dunk/filter_split/*counts.txt --appnd -i | perl -e 'while(<>){~s/_\|1//g;print}' > ${outdir}/spike_summary_counts.txt
        wait_for_jobs slamcnt
        # rm ${outdir}/dunk/filter_split/*ambiguous.bam
    fi
fi

if [ "$post" -eq 1 ]; then
    echo "$bam - alleyoop"
    mkdir -p ${outdir}/alleyoop/summary/ ${outdir}/alleyoop/rates ${outdir}/alleyoop/utrrates ${outdir}/alleyoop/tcperreadpos ${outdir}/alleyoop/tcperutrpos
    fileutilities.py T ${outdir}/dunk/count --dir 'tcount.tsv$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J alleyoop ,--wrap "'alleyoop collapse -o ${outdir}/dunk/count {abs}'"
    fileutilities.py T ${outdir}/dunk/filter --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J alleyoop ,--wrap "'alleyoop rates -o ${outdir}/alleyoop/rates -r $ref -mq $minq {abs}'"
    fileutilities.py T ${outdir}/dunk/filter --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J alleyoop ,--wrap "'alleyoop utrrates -o ${outdir}/alleyoop/utrrates -r $ref -b $bed -l $readlen {abs}'"
    fileutilities.py T ${outdir}/dunk/filter --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J alleyoop ,--wrap "'alleyoop tcperreadpos -o ${outdir}/alleyoop/tcperreadpos -r $ref -s ${outdir}/dunk/snp -l $readlen -mq $minq {abs}'"
    fileutilities.py T ${outdir}/dunk/filter --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J alleyoop ,--wrap "'alleyoop tcperutrpos -o ${outdir}/alleyoop/tcperutrpos -r $ref -b $bed -s ${outdir}/dunk/snp -l $readlen -mq $minq {abs}'"
    # wait_for_jobs alleyoop
    sbatch -J alleyoop -o /dev/null -e /dev/null --wrap "alleyoop summary -t ${outdir}/dunk/count/ -o ${outdir}/alleyoop/summary/summary.txt ${outdir}/dunk/filter/*bam"
    wait_for_jobs alleyoop

    echo "$bam - MultiQC"
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_alleyoop ${outdir}/fastqc_post ${outdir}/alleyoop/summary ${outdir}/alleyoop/rates ${outdir}/alleyoop/tcperreadpos ${outdir}/alleyoop/tcperutrpos ${outdir}/alleyoop/utrrates

    echo "$bam - RPMu"
    slamseq_rpmu.R collapsed ${outdir}/dunk/count/*tcount_collapsed.csv
    slamseq_rpmu.R utr ${outdir}/dunk/count/*tcount.tsv

    echo "$bam - Merge"
    fileutilities.py T ${outdir}/dunk/count --dir 'filtered_tcount_collapsed.rpmu.txt$' | fileutilities.py P -i -r --appnd > ${outdir}/all_collapsed_rpmu.txt
    fileutilities.py T ${outdir}/dunk/count --dir 'filtered_tcount.rpmu.txt$' | fileutilities.py P -i -r --appnd > ${outdir}/all_utr_rpmu.txt
    # After merging on rowID, the Name field needs to be deduplicated and then brought first to merge with the xref.
    dedup_table_field.R ${outdir}/all_utr_rpmu.txt Name
    fileutilities.py T ${outdir}/all_utr_rpmu_dedup.tsv -r --cols 1 0 2:$(fileutilities.py T ${outdir}/all_utr_rpmu_dedup.tsv --cntcols | cut -f 1) > ${outdir}/all_utr_rpmu.txt
    rm ${outdir}/all_utr_rpmu_dedup.tsv

    if [ -e "${xref}" ] ; then
      # Merge xref on gene identifier. But genes are not unique row identifiers for the UTR results, so I can't use my usual fileutilities/pandas.
      echo "$bam - Adding crossreferencing IDs"
      slamseq_xref.R ${outdir} 'all_utr_rpmu.txt' $xref 1 1
      slamseq_xref.R ${outdir} 'all_collapsed_rpmu.txt' $xref 1 1     # could have use fileutilities here, but consistency is simpler.
    fi
fi

# Clean up file pollution
# rm ./slurm*

echo "$bam - SLAMseq workflow finished!"
