#!/usr/bin/env sh

# Process the samples of one Run, from .sra to _barcode_report.pdf

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -d BASEDIR -b BATCHDIR -i BOWTIE2_IDX [-D DATADIR] [-p PROCESSING_DIR] [-r RESULTS_DIR] [-a AUX_DIR] [-o OFFSET] [[-s] | [-f]]"
		exit 1
}
# Defaults
data='data'
process='process'
results='results'
aux='aux'
offsets='-:0:0'
fastasuffix='fa'
XDIR='/groups/pavri/Kimon/mutpe'
five='TCTCCACAGGTGTCCACTCC'
three='GTGAGTCCTTACAACCTCTC'
adapterr=0.2
# Parse options.
while getopts 'd:D:b:p:r:a:i:o:l:F:sf' flag; do
  case "${flag}" in
    d) base="${OPTARG}" ;;        # Base dir
	D) data="${OPTARG}" ;;     		# Data dir in which the batch is located, relative to base
	b) run="${OPTARG}" ;;         # Batch dir, relative to data dir
	p) process="${OPTARG}" ;;     # Dir in which intermediate steps output files are created, relative to base
	r) results="${OPTARG}" ;;     # Dir for the final files, relative to base
	a) aux="${OPTARG}" ;;     		# Dir where the bowtie index is located, relative to base
    i) bowtie2idx="${OPTARG}" ;;  # Bowtie2 index prefix
    o) offsets="${OPTARG}" ;;     # Reference deletions: "REF:START:LENGTH" ie. "HDR2:280:6,HDR2:280:6"
    s) issra=1 ;;                 # Input is .sra format
	f) renamefq=1 ;;              # Files are uncompressed .fq instead of compressed .fastq.gz
    F) fastasuffix="${OPTARG}";;  # Reference fasta suffix. (assumes the prefix is the same as the bowtie2 index prefix)
	*) usage ;;
  esac
done

oldPath=$PATH
oldPyPath=$PYTHONPATH
export PATH="${XDIR}:${PATH}"
export PYTHONPATH="${XDIR}:${PYTHONPATH}"


wait_for_jobs(){
  sleep 60  # seconds, give time to the schedular to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep kimon.fr | grep -c $1) -eq 0 ]; then
      break
    else
      echo sleep another $((sleeptime / 60)) minutes...
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


module load pandoc/2.0.5
##module load samtools/1.9-foss-2017a
# using samtools from my conda (because it's there anyway).

if [[ ! -d "${base}/${process}/${run}" ]]; then
	mkdir -p "${base}/${process}/${run}"
  echo "Created destination: ${base}/${process}/${run}"
fi
if [[ ! -d "${base}/${results}/${run}" ]]; then
	mkdir -p "${base}/${results}/${run}"
  echo "Created destination: ${base}/${results}/${run}"
fi


if [[ $renamefq ]]; then
	echo "Renaming & compressing ${data}/${run}/*.fq to *.fastq.gz"
	fileutilities.py T ${base}/${data}/${run} --dir fq | fileutilities.py P --loop mv {abs} {dir}/{bas}.fastq \; srun gzip {dir}/{bas}.fastq
fi
if [[ $issra ]]; then
	echo "Converting ${data}/${run}/*.sra to *.fastq.gz"
	fileutilities.py T ${base}/${data}/${run} --dir .sra | fileutilities.py P --loop srun ~/utility_scripts/sra2fastq.sh ${base}/${data}/${run} ${base}/${data}/${run}/{val}
fi
if [[ ! $renamefq ]] && [[ ! $issra ]]; then
	echo "Extracting ${data}/${run}/*.bam to *_1/2.fastq.gz"
	fileutilities.py T ${base}/${data}/${run} --dir .bam | fileutilities.py P --loop srun ,-J bam2fq ,--qos short samtools fastq ,-c 1 ,-1 {dir}/{ali}_1.fastq.gz ,-2 {dir}/{ali}_2.fastq.gz {abs} \&
	wait_for_jobs bam2fq
fi

echo "Merging read pairs into ${process}/${run}/*.fastq.gz"
mkdir -p ${base}/${process}/${run}/flash
module load flash/1.2.11-foss-2017a
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop srun flash ,-t 1 ,-M 150 ,-x 0.25 ,-d ${base}/${process}/${run} ,-o {ali} ,-z ${base}/${data}/${run}/{ali}_1.fastq.gz ${base}/${data}/${run}/{ali}_2.fastq.gz 2\>\&1 \| tee ${base}/${process}/${run}/flash/{ali}.log \&
wait_for_jobs flash
module unload flash

echo "Running FastQC for ${run}/*.extendedFrags.fastq.gz"
mkdir -p ${base}/${process}/${run}/fastqc
module load fastqc/0.11.5-java-1.8.0_121
srun fastqc -o ${base}/${process}/${run}/fastqc ${base}/${process}/${run}/*.extendedFrags.fastq.gz
module unload fastqc

echo "Trimming overhangs for ${run}/*.extendedFrags.fastq.gz"
mkdir -p ${base}/${process}/${run}/cutadapt
##module load cutadapt/1.16-foss-2017a-python-2.7.13
## fileutilities won't work with python 2.7, so I'm using cutadapt from my conda where it's built for python 3.6
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop srun cutadapt ,-g $five ,-e $adapterr ,-o ${base}/${process}/${run}/{ali}.trim5.fastq ${base}/${process}/${run}/{ali}.extendedFrags.fastq.gz \> ${base}/${process}/${run}/cutadapt/{ali}_5.log \&
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop srun cutadapt ,-a $three ,-e $adapterr ,-o ${base}/${process}/${run}/{ali}.trim53.fastq.gz ${base}/${process}/${run}/{ali}.trim5.fastq \> ${base}/${process}/${run}/cutadapt/{ali}_3.log \&
wait_for_jobs cutadapt
rm ${base}/${process}/${run}/*.trim5.fastq
##module unload cutadapt

echo "Running FastQC for ${run}/*.trim53.fastq.gz"
##mkdir -p ${base}/${process}/${run}/fastqc
module load fastqc/0.11.5-java-1.8.0_121
srun fastqc -o ${base}/${process}/${run}/fastqc ${base}/${process}/${run}/*.trim53.fastq.gz
module unload fastqc

echo "Aligning ${run}/*.trim53.fastq.gz to ${aux}/${bowtie2idx}"
mkdir -p ${base}/${process}/${run}/bowtie2
module load bowtie2/2.2.9-foss-2017a
# With the overhangs of the reads trimmed off, global alignment is appropriate.
# Also apply a mapping quality filter to get rid of weirdos with lots of indels.
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop srun ,-J bowtie2 ,-c 4 bowtie2 ,-x ${base}/${aux}/${bowtie2idx} ,-U ${base}/${process}/${run}/{ali}.trim53.fastq.gz ,-p 2 ,--very-sensitive 2\> ${base}/${process}/${run}/bowtie2/{ali}.trim53.log \| samtools view ,-Sbq 20 \| samtools sort ,-n ,-@ 2 ,-O BAM ,-o ${base}/${process}/${run}/{ali}.aln.bam \&
wait_for_jobs bowtie2
module unload bowtie2

echo "Collecting MultiQC for ${run}"
mkdir -p ${base}/${process}/${run}
module load multiqc/1.3-foss-2017a-python-2.7.13
srun multiqc -f -o ${base}/${process}/${run}/multiqc ${base}/${process}/${run} ${base}/${process}/${run}/flash ${base}/${process}/${run}/fastqc ${base}/${process}/${run}/bowtie2 ${base}/${process}/${run}/cutadapt
module unload multiqc

echo "Stratifying ${run}/*.aln.bam by number of mismatches, and sorting by coordinates"
fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop srun ,-J stratBAM ,--qos short ${XDIR}/split_BAM_by_numMismatches.sh {abs} 1 2 \&
fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop srun ,-J stratBAM ,--qos short ${XDIR}/split_BAM_by_numMismatches.sh {abs} 3 10 \&
fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop srun ,-J stratBAM ,--qos short ${XDIR}/split_BAM_by_numMismatches.sh {abs} 11 30 \&
wait_for_jobs stratBAM

echo "Piling up ${run}/*\d+-\d+.bam onto ${aux}/${bowtie2idx}.${fastasuffix}"
fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.bam$' | fileutilities.py P --loop srun ,-J pileup ,--qos short samtools mpileup ,-A ,-B ,-f ${base}/${aux}/${bowtie2idx}.${fastasuffix} ,-d 1000000 {abs} \> {dir}/{bas}.pileup \&
wait_for_jobs pileup

echo "Counting ${run}/*\d+-\d+.pileup"
fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.pileup$' | fileutilities.py P --loop srun ,-J summpil ,--qos short ${XDIR}/summarize_mpileup.py {abs} \> {dir}/{bas}.stats \&
wait_for_jobs summpil

echo "Visualising ${run}/*.stats into ${results}/${run}/${run//\//_}*.html/pdf"
srun --qos short ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results}/${run} ${run/\//_} NULL no $offsets ${base}/${process}/${run}/*.stats

echo "Finished ${run}"
##module unload samtools

export PATH="${oldPath}"
export PYTHONPATH="${oldPyPath}"
