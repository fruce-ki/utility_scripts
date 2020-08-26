#!/usr/bin/env sh

# Process the samples of one Run, from .sra to _barcode_report.pdf

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -d BASEDIR -b BATCHDIR -i BOWTIE2_IDX [-o OFFSET] [-D DATADIR -p PROCESSING_DIR -r RESULTS_DIR -a AUX_DIR] [-x MM_RATE -m MIN_OVERLAP -M MAX_OVERLAP] [[-s] | [-f]] [-X SCRIPTS_DIR]"
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
refLen=361                      # fragment size
mergelen=30                     # +- nt from refLen for merge filtering
flash=1                         # merge read pairs
ovlpmm=0.10						# FLASh : Maximum allowed mismatch rate when merging overlapping mates
maxOver=250                     # FLASh : max allowed overlap of mates
minOver=80                      # FLASh : min allowed overlap of mates
five='TCTCCACAGGTGTCCACTCC'     # cutadapt : read1 5'
THREE='GGAGTGGACACCTGTGGAGA'    # cutadapt : read2 3', revcomp of $five
three='GTGAGTCCTTACAACCTCTC'    # cutadapt : read1 3'
FIVE='GAGAGGTTGTAAGGACTCAC'     # cutadapt : read2 5', revcomp of $three
adapterr=0.1                    # cutadapt : mismatches in adapter identification
minL1=200                       # cutadapt : minimum allowed read1 length post-trim
minL2=100                       # cutadapt : minimum allowed read2 length post-trim
# clip1=-5                        # cutadapt : trim from read1
# clip2=-150                      # cutadapt : trim from read2
orient="--fr"                   # bowtie2 : paired-end orientation ("" / "--fr" / "--rf" / "--ff")
# minfrag=370                     # bowtie2 : paired-end max fragment size
# maxfrag=350                     # bowtie2 : paired-end min fragment size
minalq=13						# samtools : Minimum allowed read mapping quality for contributing to the pileup
minbasq=30						# samtools : Minimum allowed base quality for contributing to mutation stats

# Parse options.
while getopts 'd:D:b:p:r:a:i:o:F:m:M:L:l:x:X:H:R:f:S:Psf' flag; do
  case "${flag}" in
    d) base="${OPTARG}" ;;        # Base dir
	D) data="${OPTARG}" ;;     		# Data dir in which the batch is located, relative to base
	b) run="${OPTARG}" ;;         # Batch dir, relative to data dir
	p) process="${OPTARG}" ;;     # Dir in which intermediate steps output files are created, relative to base
	r) results="${OPTARG}" ;;     # Dir for the final files, relative to base
	a) aux="${OPTARG}" ;;     		# Dir where the bowtie index is located, relative to base
    i) bowtie2idx="${OPTARG}" ;;  # Bowtie2 index prefix
    o) offsets="${OPTARG}" ;;     # Reference deletions: "REF:START:LENGTH" ie. "HDR2:280:6,HDR2:280:6"
    F) fastasuffix="${OPTARG}";;  # Reference fasta suffix. (assumes the prefix is the same as the bowtie2 index prefix)
    R) refLen="${OPTARG}";;
    f) mergelen="${OPTARG}";;
    m) minOver="${OPTARG}" ;;
    M) maxOver="${OPTARG}" ;;
    x) ovlpmm="${OPTARG}" ;;
    X) adapterr="${OPTARG}" ;;
    L) minL1="${OPTARG}" ;;
    l) minL2="${OPTARG}" ;;
    # C) clip1="${OPTARG}" ;;
    # c) clip2="${OPTARG}" ;;
    H) XDIR="${OPTARG}" ;;        # Path to the mutpe scripts directory (ie. the local repo clone).
    s) issra=1 ;;                 # Input is .sra format
	f) renamefq=1 ;;              # Files are uncompressed .fq instead of compressed .fastq.gz
    P) flash=0;;
    S) subs="${OPTARG}" ;;       # file specifying the substraction patterns (from_prefix \t minus_prefix)
	*) usage ;;
  esac
done

# Read count tracking
metrics="${base}/${process}/${run}/metrics.txt"

oldPath=$PATH
oldPyPath=$PYTHONPATH
export PATH="${XDIR}:${PATH}"
export PYTHONPATH="${XDIR}:${PYTHONPATH}"
echo $(which fileutilities.py)

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



# echo ""
# if [[ ! -d "${base}/${process}/${run}" ]]; then
# 	mkdir -p "${base}/${process}/${run}"
#   echo "Created destination: ${base}/${process}/${run}"
# fi
# if [[ ! -d "${base}/${results}/${run}" ]]; then
# 	mkdir -p "${base}/${results}/${run}"
#   echo "Created destination: ${base}/${results}/${run}"
# fi
#
# if [[ $renamefq ]]; then
#     echo ""
#     echo "Renaming & compressing ${data}/${run}/*.fq to *.fastq.gz"
# 	fileutilities.py T ${base}/${data}/${run} --dir fq | fileutilities.py P --loop mv {abs} {dir}/{bas}.fastq \; srun gzip {dir}/{bas}.fastq
# fi
# if [[ $issra ]]; then
#     echo ""
#     echo "Converting ${data}/${run}/*.sra to *.fastq.gz"
# 	fileutilities.py T ${base}/${data}/${run} --dir .sra | fileutilities.py P --loop srun ~/utility_scripts/sra2fastq.sh ${base}/${data}/${run} ${base}/${data}/${run}/{val}
# fi
#
#
# if [[ ! $renamefq ]] && [[ ! $issra ]]; then
#     echo ""
#     echo "Extracting ${data}/${run}/*.bam to *_1/2.fastq.gz"
# 	fileutilities.py T ${base}/${data}/${run} --dir .bam | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J bam2fq ,--qos short ,--wrap "'samtools fastq -c 1 -1 ${base}/${data}/${run}/{ali}_1.fastq.gz -2 ${base}/${data}/${run}/{ali}_2.fastq.gz {abs}'"
# 	wait_for_jobs bam2fq
# fi
#
# echo ""
# echo "Running FastQC for ${run}/*_1/2.fastq.gz"
# mkdir -p ${base}/${process}/${run}/fastqc_raw
# fileutilities.py T ${base}/${data}/${run}/*_1.fastq.gz ${base}/${data}/${run}/*_2.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/${run}/fastqc_raw {abs}'"
# wait_for_jobs fastqc
#
#
# fileutilities.py T ${base}/${data}/${run}/*fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
# wait_for_jobs fgzcnt
#
#
# echo ""
# echo "Trimming ${run}/*_1/2.fastq.gz"
# mkdir -p ${base}/${process}/${run}/trimmomatic
# mkdir -p ${base}/${process}/${run}/cutadapt
# fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J trimomma ,-c 4 ,--wrap "'trimmomatic PE -threads 4 -phred33 -summary ${base}/${process}/${run}/trimmomatic/{ali}.log ${base}/${data}/${run}/{ali}_1.fastq.gz ${base}/${data}/${run}/{ali}_2.fastq.gz ${base}/${process}/${run}/{ali}_1.trim3.fastq /dev/null ${base}/${process}/${run}/{ali}_2.trim3.fastq /dev/null SLIDINGWINDOW:5:25 TRAILING:25 MINLEN:$(( minL1 < minL2 ? minL1 : minL2 ))'"
# wait_for_jobs trimmoma
# fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J cutadapt ,-c 4 ,--wrap "'cutadapt --cores 4 -q 25 -g $five -a $three -G $FIVE -A $THREE -e $adapterr -m ${minL1}:${minL2} -o ${base}/${process}/${run}/{ali}_1.trim53.fastq.gz -p ${base}/${process}/${run}/{ali}_2.trim53.fastq.gz ${base}/${process}/${run}/{ali}_1.trim3.fastq ${base}/${process}/${run}/{ali}_2.trim3.fastq > ${base}/${process}/${run}/cutadapt/{ali}.log'"
# wait_for_jobs cutadapt
#
#
# fileutilities.py T ${base}/${process}/${run}/*trim3.fastq --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqcnt ,--wrap "'countfq.sh {val} {abs} ${metrics}_{bas}.txt'"
# wait_for_jobs fqcnt
# fileutilities.py T ${base}/${process}/${run}/*trim53.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fgzcnt ,--wrap "'countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
# wait_for_jobs fgzcnt
#
# rm ${base}/${process}/${run}/*.trim3.fastq
#
#
# echo ""
# echo "Running FastQC for ${run}/*_1/2.trim53.fastq.gz"
# mkdir -p ${base}/${process}/${run}/fastqc_posttrim
# fileutilities.py T ${base}/${process}/${run}/*trim53.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/${run}/fastqc_posttrim {abs}'"
# wait_for_jobs fastqc
#
# if [[ $flash -eq 1 ]]; then
#     echo ""
#     echo "Merging read pairs into ${process}/${run}/*.trim53.extendedFrags.fastq.gz"
# 	mkdir -p ${base}/${process}/${run}/flash
# 	fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J flash ,--wrap "'flash -t 1 -M $maxOver -m $minOver -x $ovlpmm -d ${base}/${process}/${run} -o {ali} -z ${base}/${process}/${run}/{ali}_1.trim53.fastq.gz ${base}/${process}/${run}/{ali}_2.trim53.fastq.gz 2>&1 > ${base}/${process}/${run}/flash/{ali}.log'"
# 	wait_for_jobs flash
# 	# I have no use for the uncombined and they get caught up in *_1/2.fastq.gz
# 	rm ${base}/${process}/${run}/*notCombined_1.fastq.gz ${base}/${process}/${run}/*notCombined_2.fastq.gz
#
#
#     fileutilities.py T ${base}/${process}/${run}/*extendedFrags.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
#     wait_for_jobs fgzcnt
#
#
#     echo ""
#     echo "Running FastQC for ${run}/*.extendedFrags.fastq.gz"
#     mkdir -p ${base}/${process}/${run}/fastqc_postmerge
#     fileutilities.py T ${base}/${process}/${run}/*.extendedFrags.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00  ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/${run}/fastqc_postmerge {abs}'"
#     wait_for_jobs fastqc
#
#     echo ""
#     echo "Filtering out unacceptable merged lengths from ${run}/*.extendedFrags.fastq.gz"
#     fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J cutadapt ,-c 4 ,--wrap "'cutadapt --cores 4 -m $((refLen - mergelen)) -M $((refLen + mergelen)) -o ${base}/${process}/${run}/{ali}.extendedFrags.fltr.fastq.gz ${base}/${process}/${run}/{ali}.extendedFrags.fastq.gz > ${base}/${process}/${run}/cutadapt/{ali}_mrg.log'"
#     wait_for_jobs cutadapt
#
#
#     fileutilities.py T ${base}/${process}/${run}/*extendedFrags.fltr.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
#     wait_for_jobs fgzcnt
#
#
#     echo ""
#     echo "Aligning ${run}/*.extendedFrags.fltr.fastq.gz single-end to ${aux}/${bowtie2idx}"
# 	mkdir -p ${base}/${process}/${run}/bowtie2
# 	fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:20:00 ,-J bowtie2 ,-c 4 ,--wrap "'bowtie2 --no-unal -x ${base}/${aux}/${bowtie2idx} -U ${base}/${process}/${run}/{ali}.extendedFrags.fltr.fastq.gz -p 2 --very-sensitive-local 2> ${base}/${process}/${run}/bowtie2/{ali}.trim53.log | samtools view -Sb | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/${run}/{ali}.aln.bam'"
# 	wait_for_jobs bowtie2
# else
#     echo ""
#     echo "Aligning ${run}/*.trim53.fastq.gz paired-end to ${aux}/${bowtie2idx}"
# 	mkdir -p ${base}/${process}/${run}/bowtie2
#     # Paired-end
# 	fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bowtie2 ,-c 4 ,--time=0-00:20:00 ,--wrap "'bowtie2 -p 2 --very-sensitive-local --no-unal $orient --no-discordant --no-contain --no-mixed -x ${base}/${aux}/${bowtie2idx} -1 ${base}/${process}/${run}/{ali}_1.trim53.fastq.gz -2 ${base}/${process}/${run}/{ali}_2.trim53.fastq.gz 2> ${base}/${process}/${run}/bowtie2/{ali}.trim53.log | samtools view -Sb | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/${run}/{ali}.aln.bam'"
#     # # OR Single-end
#     # fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bowtie2 ,-c 4 ,--wrap "'bowtie2 -x ${base}/${aux}/${bowtie2idx} -r ${base}/${process}/${run}/{ali}_1.trim53.fastq.gz,${base}/${process}/${run}/{ali}_2.trim53.fastq.gz -p 2 --very-sensitive-local 2> ${base}/${process}/${run}/bowtie2/{ali}.trim53.log | samtools view -Sb | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/${run}/{ali}.aln.bam'"
# 	wait_for_jobs bowtie2
#
#     # This seems to remove the non-overlapping part of read2 instead of the overlap between read1 and read2. Fail.
#     # echo ""
#     # echo "Clipping pair overlaps in ${run}/*.alnovlp.bam"
#     # fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J clipOver ,-c 4 ,--wrap "'bam clipOverlap --readName --stats --in ${base}/${process}/${run}/{ali}.alnovlp.bam --out ${base}/${process}/${run}/{ali}.aln.bam'"
#     # wait_for_jobs clipOver
#     # rm ${base}/${process}/${run}/*.alnovlp.bam
# fi
#
#
# fileutilities.py T ${base}/${process}/${run}/*.aln.bam --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J bamcnt ,--wrap "'countbam.sh {val} {abs} ${metrics}_{bas}.txt'"
# wait_for_jobs bamcnt
#
#
# echo ""
# echo "Collecting MultiQC for ${run}"
# mkdir -p ${base}/${process}/${run}/multiqc_pre
# mkdir -p ${base}/${process}/${run}/multiqc_posttrim
# srun --ntasks=1 multiqc -f -o ${base}/${process}/${run}/multiqc_pre ${base}/${process}/${run}/fastqc_raw
# if [[ $flash -eq 1 ]]; then
#     mkdir -p ${base}/${process}/${run}/multiqc_postmerge
# 	srun --ntasks=1 multiqc -f -o ${base}/${process}/${run}/multiqc_posttrim ${base}/${process}/${run}/fastqc_posttrim
#     srun --ntasks=1 multiqc -f -o ${base}/${process}/${run}/multiqc_postmerge ${base}/${process}/${run}/fastqc_postmerge ${base}/${process}/${run}/bowtie2
# else
# 	srun --ntasks=1 multiqc -f -o ${base}/${process}/${run}/multiqc_posttrim ${base}/${process}/${run}/fastqc_posttrim ${base}/${process}/${run}/bowtie2
# fi
#
# echo ""
# echo "Stratifying ${run}/*.aln.bam by number of mismatches"
# # # REMEMBER to UPDATE the call to PLOTMETRICS.R , if changing the strata
# # # REMEMBER to UPDATE the SUBSTRACTIONS command below, if changing the strata
#
# # # Original Yeap et al stratification.
# # fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--qos short ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 1 2'"
# # fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--qos short ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 3 10'"
# # fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--qos short ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 11 30'"
# # # More detailed stratification.
# # # Fewer reads have many mutations, so higher strata still need to be pooled.
# fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J stratBAM ,--qos short ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 1 1'"
# fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J stratBAM ,--qos short ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 2 2'"
# fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J stratBAM ,--qos short ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 3 3'"
# fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J stratBAM ,--qos short ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 4 4'"
# fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J stratBAM ,--qos short ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 5 15'"
# fileutilities.py T ${base}/${process}/${run} --dir 'aln.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J stratBAM ,--qos short ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 16 30'"
#
#
# fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-0:10:00 ,-J bamcnt ,--wrap "'countbam.sh {val} {abs} ${metrics}_{bas}.txt'"
# wait_for_jobs bamcnt
#
#
# echo ""
# echo "Piling up ${run}/*\d+-\d+.bam onto ${aux}/${bowtie2idx}.${fastasuffix}"
# # Apply map and base quality thresholds to reduce noise.
# # Do not penalize alignment qualities for high mismatches.
# # Do not allow unpaired reads (for unmerged use-case).
# # This has caught me out several times already, so check explicitly:
# if [ ! -e "${base}/${aux}/${bowtie2idx}.${fastasuffix}" ]; then
#     echo "${base}/${aux}/${bowtie2idx}.${fastasuffix} not found! Did you forget to link it?"
#     exit 1
# fi
# fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.bam$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J pileup ,--qos short ,--wrap "'samtools mpileup -f ${base}/${aux}/${bowtie2idx}.${fastasuffix} -d 500000000 -q $minalq -Q $minbasq {abs} > {dir}/{bas}.pileup'"
# wait_for_jobs pileup
#
# echo ""
# echo "Counting ${run}/*\d+-\d+.pileup"
# fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.pileup$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J summpil ,--qos short ,--wrap "'${XDIR}/summarize_mpileup.py -p {abs} -m {dir}/{bas}.point.stats -i {dir}/{bas}.indel.stats'"
# wait_for_jobs summpil

echo ""
echo "Visualising ${run}/*point.stats into ${results}/${run}/${run//\//_}*.html/pdf"
srun --qos short --ntasks=1 --mpi=none ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results}/${run} ${run/\//_}_point NULL yes no $offsets ${base}/${process}/${run}/*.point.stats
echo "Visualising ${run}/*indel.stats into ${results}/${run}/${run//\//_}*.html/pdf"
srun --qos short --ntasks=1 --mpi=none ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results}/${run} ${run/\//_}_indel NULL yes no $offsets ${base}/${process}/${run}/*.indel.stats
echo "BedGraphs for ${results}/${run}"
fileutilities.py T ${base}/${process}/${run} --dir bedGraph | srun --qos short --ntasks=1 --mpi=none fileutilities.py P --loop mv {abs} ${base}/${process}/${run}/tmp \&\& printf '"track type=bedGraph name=%s\n"' "{ali}" \> {abs} \&\& cat ${base}/${process}/${run}/tmp \>\> {abs}
mv ${base}/${process}/${run}/*point*bedGraph ${base}/${results}/${run}/

# echo ""
# echo "Compiling readcounts from $metrics"
# # Combine and plot all metrics files
# printf '%s\t%s\n' "file" "count" > $metrics
# cat ${metrics}_* >> $metrics
# rm ${metrics}_*
# head -n 1 $metrics > ${base}/${process}/${run}/tmp
# tail -n +2 $metrics | sort >> ${base}/${process}/${run}/tmp
# mv ${base}/${process}/${run}/tmp $metrics
# UPDATE strata HERE
# srun --qos short --ntasks=1 --mpi=none ${XDIR}/plotmetrics.R $metrics ${base}/${results}/${run}/${run/\//_}_readcounds.pdf "1-1" "2-2" "3-3" "4-4" "5-15" "16-30"


if [ ! -z "$subs" ] ; then
    echo ""
    echo "Creating substractions for ${results}/${run}"
    mkdir -p ${base}/${process}/${run}/substractions
    mkdir -p ${base}/${results}/${run}/substractions
    # UPDATE strata HERE
    srun --qos short --ntasks=1 --mpi=none fileutilities.py L $subs -V --loop fileutilities.py T 1-1 2-2 3-3 4-4 5-15 16-30 ,--loop ${XDIR}/mutPE_mutation-stats_substract.R ${base}/${process}/${run}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).{{val}}.point.stats {abs}.aln_{{val}}.point.stats {ali}.aln_{{val}}.point.stats
    echo "Visualising substractions for ${results}/${run}"
    srun --qos short --ntasks=1 --mpi=none ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results}/${run}/substractions ${run/\//_}_substracted NULL yes no $offsets ${base}/${process}/${run}/substractions/*.stats
    echo "BedGraphs for ${results}/${run}/substractions"
    fileutilities.py T ${base}/${process}/${run}/substractions --dir bedGraph | srun --qos short --ntasks=1 --mpi=none fileutilities.py P --loop mv {abs} ${base}/${process}/${run}/substractions/tmp \&\& printf '"track type=bedGraph name=%s\n"' "{ali}" \> {abs} \&\& cat ${base}/${process}/${run}/substractions/tmp \>\> {abs}
    mv ${base}/${process}/${run}/substractions/*point*bedGraph ${base}/${results}/${run}/substractions/
fi


echo ""
echo "Cleaning up file-pollution"
# # Trash (most) intermediate files
# rm ${base}/${process}/${run}/*hist
# rm ${base}/${process}/${run}/*histogram
# rm ${base}/${process}/${run}/*pileup
# fileutilities.py T ${base}/${process}/${run} --dir '\d+\-\d+.bam' | fileutilities.py P --loop rm {abs}
# if [[ $flash -eq 1 ]]; then
#     rm ${base}/${process}/${run}/*trim53.fastq.gz
#     rm ${base}/${process}/${run}/*extendedFrags*.fastq.gz
# fi
# # rm -r ${base}/${results}/${run}/*.stats_pileup_MutTypes_files
rm -r ${base}/${results}/${run}/*coverages_files
rm -r ${base}/${results}/${run}/*pileups_files
if [ ! -z "$subs" ] ; then
    # # rm -r ${base}/${results}/${run}/substractions/*.stats_pileup_MutTypes_files
    rm -r ${base}/${results}/${run}/substractions/*coverages_files
    rm -r ${base}/${results}/${run}/substractions/*pileups_files
fi
# rm ./slurm*

echo ""
echo "Finished ${run}"

export PATH="${oldPath}"
export PYTHONPATH="${oldPyPath}"
