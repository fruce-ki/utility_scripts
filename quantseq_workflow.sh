#!/bin/bash

set -x

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i BAM_DIR -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR -l READ_LEN [-x XREF] [-u UMI_LENGTH] [-q MINALQ] [-T TRIM_5] [-t TRIM_3]"
    exit 1
}
# Parse options.
minq=1    # featurecounts : min alignment quality
minq3=20   # cutadapt : min base quality at 3'
minlen=20  # cutadapt : min read length after adaptor trimming
adaptor='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'  # quantFWD read 1, single end 5'adaptor
pre=1
dunk=1
post=1
umilen=0
mm_umi=0
debam=1
while getopts 'i:r:b:a:o:q:t:T:x:u:U:l:F123' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;	            # folder of BAMs
    F) debam=0;;                            # files are already fastq
    r) ref="${OPTARG}" ;;	            # whole reference fasta for alignment
    b) bed="${OPTARG}" ;;		    # bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;		    # GTF annotation (for counting)
    o) outdir="${OPTARG}" ;;		    # output directory
    q) minq="${OPTARG}" ;;	            # minimum alignement quality (1)
    t) trim3="${OPTARG}" ;;                 # cutadapt trim this many bases from end of reads
    T) trim5="${OPTARG}" ;;                 # cutadapt trim this many bases from start of reads. Applied *after* UMI extraction, if any.
    x) xref="${OPTARG}" ;;                  # table with additional IDs
    u) umilen="${OPTARG}" ;;                # Length of UMI at 5' of read. If non-zero, will be trimmed and used to deduplicate
    U) mm_umi="${OPTARG}" ;;                # Mismatch allowance in UMI
    l) rlen="${OPTARG}" ;;                  # Read length
    1) pre=0 ;;         		    # Skip bam 2 FQ and trim
    2) dunk=0 ;;        		    # Skip slamdunk
    3) post=0 ;;		            # Skip post-alignemnt stuff
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

if [ ! -z "$trim3" ]; then
    trim3=",-u -${trim3}"
fi
if [ ! -z "$trim5" ]; then
    trim5=",-u +${trim5}"
fi


prefix=$(basename $bam)
prefix=${prefix/.bam/}

mkdir -p ${outdir}/fastqc_pre ${outdir}/fastqc_post ${outdir}/fastq ${outdir}/fastq_trimmed ${outdir}/dunk/map ${outdir}/dunk/filter ${outdir}/counts ${outdir}/multiqc ${outdir}/multiqc_pre ${outdir}/multiqc_post


if [ "$pre" -eq 1 ]; then
    if ["$debam" -eq 1 ]; then
        echo "$bam - BAM2FQ"
        fileutilities.py T $bam --dir 'bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bam2fq ,--get-user-env ,--wrap "'samtools bam2fq -0 ${outdir}/fastq/{bas}.fastq {abs}'"
        wait_for_jobs bam2fq
        fileutilities.py T ${outdir}/fastq --dir 'fastq$|fq$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqzip ,--get-user-env ,--wrap "'gzip -f {abs}'"
        wait_for_jobs fqzip
    else
        fileutilities.py T $bam --dir 'fastq$|fq$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqzip ,--get-user-env ,--wrap "'gzip -f {abs}'"
        wait_for_jobs fqzip
        fileutilities.py T $bam --dir 'fastq.gz$|fq.gz$' | fileutilities.py P --link "${outdir}/fastq/"
    fi

    fqdir="${outdir}/fastq/"

    echo "$bam - UMI"
    if [ "$umilen" -gt 0 ]; then
        mkdir -p ${outdir}/dedup

        echo "$bam - Retrieving UMIs"
        repl() {        # Create a UMI pattern of appropriate length.
            printf "N"'%.s' $(seq 1 $1);
        }
        fileutilities.py T ${fqdir} --dir 'fastq.gz$|fq.gz$' | fileutilities.py P --loop sbatch ,--qos=medium ,-J umitrim ,-o /dev/null ,-e /dev/null umi_tools extract ,-I {abs} ,-S ${outdir}/dedup/{cor}_umi-clipped.fastq.gz ,-p $(repl $umilen) ,--extract-method=string ,--quality-encoding=phred33
        wait_for_jobs umitrim

        fqdir="${outdir}/dedup/"
    fi

    echo "$bam - TRIM"
    fileutilities.py T ${fqdir} --dir 'fastq.gz$|fq.gz$' | fileutilities.py P --loop sbatch ,-J cutadapt ,-o ${outdir}/fastq_trimmed/{cor}.cutadapt.log ,-e ${outdir}/fastq_trimmed/{cor}.cutadapt.log ,--get-user-env cutadapt ,-a A{18} ,-g T{18} $trim5 $trim3 ,-q $minq3 ,-m $minlen ,-o ${outdir}/fastq_trimmed/{cor}_trimmed.fastq.gz {abs}
    wait_for_jobs cutadapt

    echo "$bam - QC"
    fileutilities.py T ${fqdir} --dir 'fastq.gz$|fq.gz$' | fileutilities.py P --loop sbatch ,-J fastQC ,-o /dev/null ,-e /dev/null ,--get-user-env ,--wrap "'fastqc -o ${outdir}/fastqc_pre {abs}'"
    fileutilities.py T ${fqdir} --dir 'fastq.gz$|fq.gz$' | fileutilities.py P --loop sbatch ,-J fastQC ,-o /dev/null ,-e /dev/null ,--get-user-env ,--wrap "'fastqc -o ${outdir}/fastqc_post {abs}'"
    wait_for_jobs fastQC
    sbatch -J multiqc -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_pre ${outdir}/fastqc_pre
    sbatch -J multiqc -o /dev/null -e /dev/null multiqc -f -o ${outdir}/multiqc_post ${outdir}/fastqc_post
fi

bamdir="${outdir}/dunk/filter"  # need it accessible outside the block
threads=4                       # also
memory='40G'

if [ "$dunk" -eq 1 ]; then
    if ! [ -x "$(command -v slamdunk)" ]; then
        echo "Error: slamdunk could not be found. Did you load the module?" >&2
        exit 3
    fi

    fqdir="${outdir}/fastq_trimmed"

    fileutilities.py T $fqdir --dir 'fastq.gz$|fq.gz$' | perl -e 'while(<>){~s/\.fastq$|\.fq$/\tdummy\t0/; print}' > ${outdir}/dunk/samplesheet.tsv
    threads=$( wc -l ${outdir}/dunk/samplesheet.tsv | perl -e 'while(<>){~/^(\d+)/; print $1}' )
    # memory=$((4 * $threads))
    if [ "$threads" -gt 24 ]; then
      threads=24
    fi

    # echo "$bam - MAP"
    # sbatch -J slamaln -o /dev/null -e /dev/null --qos=medium -c $threads --mem=$memory --wrap "slamdunk map -t $threads -r $ref -o ${outdir}/dunk/map -5 0 -n 100 --quantseq -ss ${outdir}/dunk/samplesheet.tsv"
    # wait_for_jobs slamaln

    # echo "$bam - FILTER"
    # sbatch -J slamfltr -c $threads --mem=$memory -o /dev/null -e /dev/null --wrap "slamdunk filter -t $threads -o ${outdir}/dunk/filter -b $bed ${outdir}/dunk/map/*.bam"
    # wait_for_jobs slamfltr

    if [ "$umilen" -gt 0 ]; then
        mkdir -p ${outdir}/dedup ${outdir}/dedup/tmp ${outdir}/dedup/map_sorted

        # echo "$bam - Sort and index"
        # fileutilities.py T $bamdir/*bam --loop sbatch ,-J samsort ,-o /dev/null ,-e /dev/null ,--wrap "'samtools sort -o ${outdir}/dedup/map_sorted/{cor}.bam {abs}'"
        # wait_for_jobs samsort
        # fileutilities.py T ${outdir}/dedup/map_sorted/*bam --loop sbatch ,-J samidx ,-o /dev/null ,-e /dev/null ,--wrap "'samtools index {abs}'"
        # wait_for_jobs samidx

        echo ""
        echo "$bam - Removing duplicates"
        fileutilities.py T ${outdir}/dedup/map_sorted/*.bam --loop sbatch ,--qos=medium ,--mem=$memory ,-J umidedup umi_tools dedup ,-I {abs} ,-S ${outdir}/dedup/{cor}.deduped.bam ,-L ${outdir}/dedup/{cor}_dedup.log  ,--edit-distance-threshold=${mm_umi}
        wait_for_jobs umidedup
        #  ,-o /dev/null ,-e /dev/null
    fi
fi

if [ "$umilen" -gt 0 ]; then    # So that it points correctly even if I'm repeating only the 3rd step
  bamdir="${outdir}/dedup"
fi

if [ "$post" -eq 1 ]; then
    threads=$( wc -l ${outdir}/dunk/samplesheet.tsv | perl -e 'while(<>){~/^(\d+)/; print $1}' )
    if [ "$threads" -gt 24 ]; then
      threads=24
    fi

    echo "$bam - QUANTIFY alignments with qual >$minq"
    fileutilities.py T ${bamdir}/*.bam --loop sbatch ,-J ftrcnt ,-o /dev/null ,-e /dev/null ,-c 4 ,--wrap "'featureCounts -T 4 -M --primary -Q $minq -a $annot -o ${outdir}/counts/{cor}.txt {abs}'"
    wait_for_jobs ftrcnt

    echo "$bam - MultiQC"
    if [ "$umilen" -gt 0 ]; then
        sbatch -o /dev/null -e /dev/null -J multiqc multiqc -f -o ${outdir}/multiqc ${outdir}/fastqc_pre ${outdir}/fastqc_post ${outdir}/dedup ${outdir}/fastq_trimmed ${outdir}/counts
    else
        sbatch -o /dev/null -e /dev/null -J multiqc multiqc -f -o ${outdir}/multiqc ${outdir}/fastqc_pre ${outdir}/fastqc_post ${outdir}/fastq_trimmed ${outdir}/counts
    fi

    echo "$bam - Collect all"
    fileutilities.py T ${outdir}/counts/*trimmed.txt -i -l --cols 6 | perl -e 'while(<>){~s/(_umi-clipped)?(_trimmed)?(\.fa?s?t?q_slamdunk_mapped_filtered)?_\|6//g;print}' > ${outdir}/all_counts.tsv
    rpm.R -f ${outdir}/all_counts.tsv -i 1

    if [ ! -z "$xref" ]; then
      # Unlike slamseq, the output of featureCounts is per gene/Entrez, not per UTR, so here Entrez are unique per row. No need for the special xref.R .
      echo "$bam - Cross-reference"
      fileutilities.py T $xref ${outdir}/all_counts.tsv -i -r --appnd > ${outdir}/all_counts_xref.txt
      fileutilities.py T $xref ${outdir}/all_counts_rpm.txt -i -r --appnd > ${outdir}/all_counts_rpm_xref.txt
    fi
fi


echo "$bam - Quantseq workflow finished!"
