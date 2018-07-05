#!/bin/bash

# Send array job for STAR. 
# SGE_TASK_ID <-> replicate number


# Set defaults.
log=false
threads=1
maxram=1000000000
tool=/homes/mgierlinski/software/STAR_2.5.0a/STAR
project=$(pwd)
cond=''
fastqdir='fastq'
indexdir='starindex'
index=''
annotdir='genome'
annot=''
samdir='sam'
scriptsdir="${HOME}/scripts"

function usage() {
    echo "Usage:"
    echo "      $0 -c PREFIX -a ANNOTATION -n INDEX  [-A ANNOTDIR] [-i FASTQDIR] [-o SAMDIR] [-N INDEXDIR] [-t THREADS] [-m MAXMEMORY] [-X path/star] [-D PROJECTDIR] [-L]"
    echo "Defaults:"
    echo "      -A $annotdir -i $fastqdir -o $samdir -N ${indexdir} -t ${threads} -m $maxram -X $tool -D $project -L $log"  
    echo "All directories must be relative to PROJECTDIR. -L requires ${scriptsdir}/mylogs.py"
    exit 1
}




# Parse options.
while getopts 'c:i:o:n:N:a:A:t:m:X:D:L' flag; do
  case "${flag}" in
    c) cond="${OPTARG}" ;;
    i) fastqdir="${OPTARG}" ;;
    o) samdir="${OPTARG}" ;;
    n) index="${OPTARG}" ;;
    N) indexdir="${OPTARG}" ;;
    a) annot="${OPTARG}" ;;
    A) annotdir="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    m) maxram=="${OPTARG}" ;;
    X) tool="${OPTARG}" ;;
    D) project="${OPTARG}" ;;
    L) log=true ;;
    *) usage ;;
  esac
done

# No point going further without at least these 3.
[[ -z "$cond" ]] && usage
[[ -z "$annot" ]] && usage
[[ -z "$index" ]] && usage

# Tidy up the arguments.
myparams="--runThreadN ${threads} --limitBAMsortRAM ${maxram} --outFileNamePrefix ${project}/${samdir}/${cond}_${SGE_TASK_ID}_ --genomeDir ${project}/${indexdir}/${index} --sjdbGTFfile ${project}/${annotdir}/${annot}"
miscparams="--quantMode GeneCounts --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory --readFilesCommand zcat --outReadsUnmapped Fastx --outSAMmode Full --outSAMstrandField intronMotif --outSAMunmapped Within --outSJfilterIntronMaxVsReadN 5000 10000 15000 20000 --outWigType None --outWigNorm None --outFilterType BySJout --outFilterMultimapNmax 2 --outFilterMismatchNmax 5"

# Paired files or singles?
input=''
r1="${project}/${fastqdir}/${cond}_${SGE_TASK_ID}_R1.fastq.gz"
r2="${project}/${fastqdir}/${cond}_${SGE_TASK_ID}_R2.fastq.gz"
if [ -f $r2 ] ;
    then
        input="--readFilesIn $r1 $r2"
    else
        input="--readFilesIn $r1"
fi


# Log and Execute.

command="$tool $input $myparams $miscparams"
if [ "$log" = true ] ; then  
	python "${scriptsdir}/mylogs.py" "$command"
fi

$command
 
 
 
 
 

  
 

