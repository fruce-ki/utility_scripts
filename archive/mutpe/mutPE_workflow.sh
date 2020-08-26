#!/usr/bin/env sh

echo "MutPE Workflow"
# Process the samples of one Run, from .bam/.fq

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
unpaired=0                      # ignore pairs and align single-end
ovlpmm=0.10						# FLASh : Maximum allowed mismatch rate when merging overlapping mates
maxOver=250                     # FLASh : max allowed overlap of mates
minOver=80                      # FLASh : min allowed overlap of mates
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
                                          # Reads are trimmed for base quality at 25 (hard-coded), so that's kind of a floor for $minbasq
dopre=1
doal=1
dopile=1
doviz=1
dopost=1
allowoutties=0
vdj="B18"       # Default
allelecutoff=1

# Parse options.
while getopts 'd:D:b:p:r:a:i:o:F:m:M:L:l:x:X:H:R:f:S:V:A:OPUsf12345' flag; do
  case "${flag}" in
    d) base="${OPTARG}" ;;        # Base dir
    b) run="${OPTARG}" ;;         # Batch subdir. If defined, it will be appended to $data, $process and $results
	D) data="${OPTARG}" ;;        # BAM Data dir relative to base
	p) process="${OPTARG}" ;;     # Dir in which intermediate steps output files are created, relative to base
	r) results="${OPTARG}" ;;     # Dir for the final files, relative to base
	a) aux="${OPTARG}" ;;     		# Dir where the bowtie index is located, relative to base
    i) bowtie2idx="${OPTARG}" ;;  # Bowtie2 index prefix
    o) offsets="${OPTARG}" ;;     # Reference deletions: "REF:START:LENGTH" ie. "HDR2:2301:6,HDR4:2301:6"
    F) fastasuffix="${OPTARG}";;  # Reference fasta suffix. (assumes the prefix is the same as the bowtie2 index prefix)
    R) refLen="${OPTARG}";;
    A) allelecutoff="${OPTARG}";; # Ignore mutations with frequencies higher than this, as likely clonal SNPs.
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
    P) flash=0;;                  # paired-end instead of merged pairs
    O) allowoutties=1;;           # also allow read pairs that overlap on their 5' ends instead of their 3' ends.
    U) unpaired=1;;                 # single-end. Only applicable if also -P
    S) subs="${OPTARG}" ;;       # file specifying the substraction patterns (from_prefix \t minus_prefix)
    1) dopre=0;;                   # skip all steps before alignment, up to and including trimming
    2) doal=0;;                    # Skip alignment, including read pair merge and filtering
    3) dopile=0;;                  # Skip stratification and pileup
    4) doviz=0;;                   # Skip visualisations
    5) dopost=0;;                  # Skip multiqc and collective readcounts
    V) vdj="${OPTARG}" ;;          # Configure adapter trimming and coordinates range for B18, CH12, Ramos1, Ramos2, Ramos3
	*) usage ;;
  esac
done

if [ ! -z "$run" ]; then
    data="${data}/${run}"
    process="${process}/${run}"
    results="${results}/${run}"
fi

# Configure adapters for trimming
if [ "$vdj" == 'B18' ] || [ "$vdj" == 'CH12' ]; then
    # First-round PCR primer 3' ends (also trimming the ends of the amplicon)
    five='TCTCCACAGGTGTCCACTCC'     # cutadapt : read1 5'
    three='GTGAGTCCTTACAACCTCTC'    # cutadapt : read1 3', revcomp of $FIVE
    FIVE='GAGAGGTTGTAAGGACTCAC'     # cutadapt : read2 5'
    THREE='GGAGTGGACACCTGTGGAGA'    # cutadapt : read2 3', revcomp of $five
    # echo "Mouse VDJ adapters"
elif [ "$vdj" == 'B18_1' ]; then    # for round 1, where the amplicons were different
    # Two amplicons in the same sequencing pool/barcode, so two primers to trim
    five='CTTTCTCTCCACAGGTGTCCACTCC -g GGATTGATCCTAATAGTGGTGGTAC'     # cutadapt : read1 5'
    three='GTTCAAGAGCAAGGCCACACTGAC -a GGTGAGTCCTTACAACCTCTCTCTT'    # cutadapt : read1 3', revcomp of $FIVE
    FIVE='GTCAGTGTGGCCTTGCTCTTGAAC -G AAGAGAGAGGTTGTAAGGACTCACC'     # cutadapt : read2 5'
    THREE='GGAGTGGACACCTGTGGAGAGAAAG -A GTACCACCACTATTAGGATCAATCC'    # cutadapt : read2 3', revcomp of $five
elif [ "$vdj" == 'CH12_1' ]; then    # for round 1, where the amplicons were different
    # Two amplicons in the same sequencing pool/barcode, so two primers to trim
    five='CTTTCTCTCCACAGGTGTCCACTCC -g ATCCTAGCAATGGTGGTACTAACTAC'     # cutadapt : read1 5'
    three='GTTCAAGAGCAAGGCCACACTGAC -a GGTGAGTCCTTACAACCTCTCTCTT'    # cutadapt : read1 3', revcomp of $FIVE
    FIVE='GTCAGTGTGGCCTTGCTCTTGAAC -G AAGAGAGAGGTTGTAAGGACTCACC'     # cutadapt : read2 5'
    THREE='GGAGTGGACACCTGTGGAGAGAAAG -A GTAGTTAGTACCACCATTGCTAGGAT'    # cutadapt : read2 3', revcomp of $five
elif [ "$vdj" == 'Ramos1' ] || [ "$vdj" == 'Ramos2' ] || [ "$vdj" == 'Ramos3' ] || [ "$vdj" == 'Ramos' ]; then
    # # Middle of first-round PCR primer (preserve the amplicon-specific 3' ends)
    five='TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
      FIVE='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    three='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    THREE='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA'
    # # Just the shared middle part
    # five='GCTCTTCCGATCT'
    # FIVE='GCTCTTCCGATCT'
    # three='AGATCGGAAGAGC'
    # THREE='AGATCGGAAGAGC'
    # echo "Human VDJ adapters"
else
    echo "$vdj is not among the pre-defined configurations for adapter trimming."
    exit 1
fi


# Read count tracking
metrics="${base}/${process}/metrics.txt"

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



echo ""
if [[ ! -d "${base}/${process}" ]]; then
	mkdir -p "${base}/${process}"
  echo "Created destination: ${base}/${process}"
fi
if [[ ! -d "${base}/${results}" ]]; then
	mkdir -p "${base}/${results}"
  echo "Created destination: ${base}/${results}"
fi



if [[ "$dopre" -eq 1 ]]; then
    if [[ $renamefq ]]; then
        echo ""
        echo "Renaming & compressing ${data}/*.fq to *.fastq.gz"
    	${XDIR}/fileutilities.py T ${base}/${data} --dir fq | ${XDIR}/fileutilities.py P --loop mv {abs} {dir}/{bas}.fastq \; srun gzip {dir}/{bas}.fastq
    fi
    # if [[ $issra ]]; then
    #     echo ""
    #     echo "Converting ${data}/*.sra to *.fastq.gz"
    # 	${XDIR}/fileutilities.py T ${base}/${data} --dir .sra | ${XDIR}/fileutilities.py P --loop srun ${XDIR}/sra2fastq.sh ${base}/${data} ${base}/${data}/{val}
    # fi


    if [[ ! $renamefq ]] && [[ ! $issra ]]; then
        echo ""
        echo "Extracting ${data}/*.bam to *_1/2.fastq.gz"
    	${XDIR}/fileutilities.py T ${base}/${data} --dir .bam | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J bam2fq ,--wrap "'samtools fastq -c 1 -1 ${base}/${data}/{ali}_1.fastq.gz -2 ${base}/${data}/{ali}_2.fastq.gz {abs}'"
    	wait_for_jobs bam2fq
    fi

    echo ""
    echo "Running FastQC for *_1/2.fastq.gz"
    mkdir -p ${base}/${process}/fastqc_raw
    ${XDIR}/fileutilities.py T ${base}/${data}/*_1.fastq.gz ${base}/${data}/*_2.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/fastqc_raw {abs}'"
    wait_for_jobs fastqc


    ${XDIR}/fileutilities.py T ${base}/${data}/*fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs fgzcnt


    echo ""
    echo "Trimming *_1/2.fastq.gz"
    mkdir -p ${base}/${process}/trimmomatic
    mkdir -p ${base}/${process}/cutadapt
    ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J trimomma ,-c 4 ,--wrap "'trimmomatic PE -threads 4 -phred33 -summary ${base}/${process}/trimmomatic/{ali}.log ${base}/${data}/{ali}_1.fastq.gz ${base}/${data}/{ali}_2.fastq.gz ${base}/${process}/{ali}_1.trim3.fastq /dev/null ${base}/${process}/{ali}_2.trim3.fastq /dev/null SLIDINGWINDOW:25:25 TRAILING:25 MINLEN:$(( minL1 < minL2 ? minL1 : minL2 ))'"
    wait_for_jobs trimmoma
    ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J cutadapt ,-c 4 ,--wrap "'cutadapt --cores 4 -q 25 -g $five -a $three -G $FIVE -A $THREE -e $adapterr -m ${minL1}:${minL2} -o ${base}/${process}/{ali}_1.trim53.fastq.gz -p ${base}/${process}/{ali}_2.trim53.fastq.gz ${base}/${process}/{ali}_1.trim3.fastq ${base}/${process}/{ali}_2.trim3.fastq > ${base}/${process}/cutadapt/{ali}.log'"
    wait_for_jobs cutadapt


    ${XDIR}/fileutilities.py T ${base}/${process}/*trim3.fastq --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqcnt ,--wrap "'${XDIR}/countfq.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs fqcnt
    ${XDIR}/fileutilities.py T ${base}/${process}/*trim53.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs fgzcnt

    rm ${base}/${process}/*.trim3.fastq


    echo ""
    echo "Running FastQC for *_1/2.trim53.fastq.gz"
    mkdir -p ${base}/${process}/fastqc_posttrim
    ${XDIR}/fileutilities.py T ${base}/${process}/*trim53.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/fastqc_posttrim {abs}'"
    wait_for_jobs fastqc
fi

if [[ "$doal" -eq 1 ]]; then
    if [[ "$flash" -eq 1 ]]; then
        echo ""
        echo "Merging read pairs into ${process}/*.trim53.extendedFrags.fastq.gz"
        mkdir -p ${base}/${process}/flash
    	if [[ "$allowoutties" -eq 0 ]]; then
        	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J flash ,--wrap "'flash -t 1 -M $maxOver -m $minOver -x $ovlpmm -d ${base}/${process} -o {ali} -z ${base}/${process}/{ali}_1.trim53.fastq.gz ${base}/${process}/{ali}_2.trim53.fastq.gz 2>&1 > ${base}/${process}/flash/{ali}.log'"
        else
            ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J flash ,--wrap "'flash -O -t 1 -M $maxOver -m $minOver -x $ovlpmm -d ${base}/${process} -o {ali} -z ${base}/${process}/{ali}_1.trim53.fastq.gz ${base}/${process}/{ali}_2.trim53.fastq.gz 2>&1 > ${base}/${process}/flash/{ali}.log'"
        fi
        wait_for_jobs flash
    	# I have no use for the uncombined and they get caught up in *_1/2.fastq.gz
    	rm ${base}/${process}/*notCombined_1.fastq.gz ${base}/${process}/*notCombined_2.fastq.gz


        ${XDIR}/fileutilities.py T ${base}/${process}/*extendedFrags.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
        wait_for_jobs fgzcnt


        echo ""
        echo "Running FastQC for *.extendedFrags.fastq.gz"
        mkdir -p ${base}/${process}/fastqc_postmerge
        ${XDIR}/fileutilities.py T ${base}/${process}/*.extendedFrags.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00  ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/fastqc_postmerge {abs}'"
        wait_for_jobs fastqc

        echo ""
        echo "Filtering out unacceptable merged lengths from *.extendedFrags.fastq.gz"
        ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J cutadapt ,-c 4 ,--wrap "'cutadapt --cores 4 -m $((refLen - mergelen)) -M $((refLen + mergelen)) -o ${base}/${process}/{ali}.extendedFrags.fltr.fastq.gz ${base}/${process}/{ali}.extendedFrags.fastq.gz > ${base}/${process}/cutadapt/{ali}_mrg.log'"
        wait_for_jobs cutadapt


        ${XDIR}/fileutilities.py T ${base}/${process}/*extendedFrags.fltr.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
        wait_for_jobs fgzcnt


        echo ""
        echo "Aligning *.extendedFrags.fltr.fastq.gz single-end to ${aux}/${bowtie2idx}"
    	mkdir -p ${base}/${process}/bowtie2
    	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-01:00:00 ,-J bowtie2 ,-c 4 ,--wrap "'bowtie2 --no-unal -x ${base}/${aux}/${bowtie2idx} -U ${base}/${process}/{ali}.extendedFrags.fltr.fastq.gz -p 2 --very-sensitive-local 2> ${base}/${process}/bowtie2/{ali}.trim53.log | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/{ali}.aln.bam'"
    	wait_for_jobs bowtie2
    else
        echo ""
        mkdir -p ${base}/${process}/bowtie2
        if [[ $unpaired -eq 0 ]]; then
            # Paired-end
            echo "Aligning *.trim53.fastq.gz paired-end to ${aux}/${bowtie2idx}"
        	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bowtie2 ,-c 4 ,--time=0-01:00:00 ,--wrap "'bowtie2 -p 2 --very-sensitive-local --no-unal --no-discordant --no-contain --no-mixed -x ${base}/${aux}/${bowtie2idx} -1 ${base}/${process}/{ali}_1.trim53.fastq.gz -2 ${base}/${process}/{ali}_2.trim53.fastq.gz 2> ${base}/${process}/bowtie2/{ali}.trim53.log | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/{ali}.aln.bam'"

            # This seems to remove the non-overlapping part of read2 instead of the overlap between read1 and read2. Fail.
            # echo ""
            # echo "Clipping pair overlaps in *.alnovlp.bam"
            # ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J clipOver ,-c 4 ,--wrap "'bam clipOverlap --readName --stats --in ${base}/${process}/{ali}.alnovlp.bam --out ${base}/${process}/{ali}.aln.bam'"
            # wait_for_jobs clipOver
            # rm ${base}/${process}/*.alnovlp.bam
        else
            # OR Single-end
            echo "Aligning *.trim53.fastq.gz single-end to ${aux}/${bowtie2idx}"
        	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J bowtie2 ,-c 4 ,--time=0-01:00:00 ,--wrap "'bowtie2 -p 2 --very-sensitive-local --no-unal -x ${base}/${aux}/${bowtie2idx} -U ${base}/${process}/{ali}_1.trim53.fastq.gz,${base}/${process}/{ali}_2.trim53.fastq.gz 2> ${base}/${process}/bowtie2/{ali}.trim53.log | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/{ali}.aln.bam'"
        fi
    	wait_for_jobs bowtie2
    fi


    ${XDIR}/fileutilities.py T ${base}/${process}/*.aln.bam --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J bamcnt ,--wrap "'${XDIR}/countbam.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs bamcnt
fi


if [[ "$dopile" -eq 1 ]]; then
    echo ""
    echo "Stratifying *.aln.bam by number of mismatches"
    # # REMEMBER to UPDATE the call to PLOTMETRICS.R , if changing the strata
    # # REMEMBER to UPDATE the SUBSTRACTIONS command below, if changing the strata

    # # Original Yeap et al stratification.
    # ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 1 2'"
    # ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 3 10'"
    # ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 11 30'"
    # # More detailed stratification.
    # # Fewer reads have many mutations, so higher strata still need to be pooled.
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 1 1'"
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 2 2'"
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 3 3'"
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 4 4'"
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 5 15'"
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 16 30'"
    wait_for_jobs stratBAM


    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln(\.\d+-\d+)?.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-0:15:00 ,-J bamcnt ,--wrap "'${XDIR}/countbam.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs bamcnt


    echo ""
    echo "Piling up *\d+-\d+.bam onto ${aux}/${bowtie2idx}.${fastasuffix}"
    # Apply map and base quality thresholds to reduce noise.
    # Do not penalize alignment qualities for high mismatches.
    # Do not allow unpaired reads (for unmerged use-case).
    if [ ! -e "${base}/${aux}/${bowtie2idx}.${fastasuffix}" ]; then
        # This has caught me out several times already, so check explicitly
        echo "${base}/${aux}/${bowtie2idx}.${fastasuffix} not found! Did you forget to link it there?"
        exit 1
    fi
    # mpileup wants them sorted, it doesn't know it's an amplicon. But I don't need both files, and I don't want the .sorted suffix in everything downstream.
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln(\.\d+-\d+)?.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J samsort ,--wrap "'samtools sort {abs} > {dir}/{bas}.sorted.bam'"
    wait_for_jobs samsort
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln(\.\d+-\d+)?.sorted.bam$' | perl -e 'while(<>){~s/.sorted$//;print}' | ${XDIR}/fileutilities.py P --loop mv {abs} {dir}/{ali}.bam
    ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln(\.\d+-\d+)?.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,--mem=10G ,-o /dev/null ,-e /dev/null ,--time=0-00:15:00 ,-J sampile ,--wrap "'samtools mpileup -f ${base}/${aux}/${bowtie2idx}.${fastasuffix} -d 5000000 -q $minalq -Q $minbasq {abs} > {dir}/{bas}.pileup'"
    wait_for_jobs sampile

    echo ""
    echo "Counting *\d+-\d+.pileup"
    ${XDIR}/fileutilities.py T ${base}/${process} --dir '\.pileup$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J summpil ,--wrap "'${XDIR}/summarize_mpileup.py -p {abs} -m {dir}/{bas}.point.stats -i {dir}/{bas}.indel.stats'"
    wait_for_jobs summpil
fi

if [[ "$dopost" -eq 1 ]]; then
    echo ""
    echo "Compiling readcounts from $metrics"
    # Combine and plot all metrics files
    printf '%s\t%s\n' "file" "count" > $metrics
    cat ${metrics}_* >> $metrics
    #rm ${metrics}_*
    head -n 1 $metrics > ${base}/${process}/tmp
    tail -n +2 $metrics | sort >> ${base}/${process}/tmp
    mv ${base}/${process}/tmp $metrics
    # UPDATE strata HERE
    srun --ntasks=1 --mpi=none ${XDIR}/plotmetrics.R $metrics ${base}/${results}/${run/\//_}_readcounds.pdf "1-1" "2-2" "3-3" "4-4" "5-15" "16-30"

    echo ""
    echo "Collecting MultiQC for ${run}"
    mkdir -p ${base}/${process}/multiqc_pre
    mkdir -p ${base}/${process}/multiqc_posttrim
    srun --ntasks=1 multiqc -f -o ${base}/${process}/multiqc_pre ${base}/${process}/fastqc_raw
    if [[ $flash -eq 1 ]]; then
        mkdir -p ${base}/${process}/multiqc_postmerge
        srun --ntasks=1 multiqc -f -o ${base}/${process}/multiqc_posttrim ${base}/${process}/fastqc_posttrim
        srun --ntasks=1 multiqc -f -o ${base}/${process}/multiqc_postmerge ${base}/${process}/fastqc_postmerge ${base}/${process}/bowtie2
    else
        srun --ntasks=1 multiqc -f -o ${base}/${process}/multiqc_posttrim ${base}/${process}/fastqc_posttrim ${base}/${process}/bowtie2
    fi
fi

if [[ "$doviz" -eq 1 ]]; then
    echo ""
    echo "Visualising *point.stats into ${results}/${run//\//_}*.html/pdf"
    srun --ntasks=1 --mpi=none ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results} ${run/\//_}_point NULL yes no $offsets $vdj $allelecutoff ${base}/${process}/*.point.stats
    echo "Visualising *indel.stats into ${results}/${run//\//_}*.html/pdf"
    # Output only the unstratified indels. It doesn't make much sense to stratify indels by point mutation load.
    srun --ntasks=1 --mpi=none ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results} ${run/\//_}_indel NULL yes no $offsets $vdj $allelecutoff ${base}/${process}/*aln.indel.stats
    echo "BedGraphs for ${results}"
    # Add track header
    ${XDIR}/fileutilities.py T ${base}/${process} --dir bedGraph | srun --ntasks=1 --mpi=none ${XDIR}/fileutilities.py P --loop mv {abs} ${base}/${process}/tmp \&\& printf '"track type=bedGraph name=%s\n"' "{ali}" \> {abs} \&\& cat ${base}/${process}/tmp \>\> {abs}
    mv ${base}/${process}/*point*bedGraph ${base}/${results}/

    if [ ! -z "$subs" ] ; then
        echo ""
        echo "Creating substractions for ${results}"
        mkdir -p ${base}/${process}/substractions
        mkdir -p ${base}/${results}/substractions
        # UPDATE strata HERE
        # Looping cannot take just '' or '.' target values, tey get auto-substituted for current directory, which breaks the unstratified.
        # So I have to keep some more of the filename, even though it's repetitive.
        srun --ntasks=1 --mpi=none ${XDIR}/fileutilities.py L $subs -V --loop ${XDIR}/fileutilities.py T point 1-1.point 2-2.point 3-3.point 4-4.point 5-15.point 16-30.point ,--loop ${XDIR}/mutPE_mutation-stats_substract.R ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).{{val}}.stats {abs}.aln.{{val}}.stats {ali}.aln.{{val}}.stats
        echo "Visualising substractions for ${results}"
        srun --ntasks=1 --mpi=none ${XDIR}/mutPE_mutation-stats_viz.R ${base}/${results}/substractions ${run/\//_}_substracted NULL yes no $offsets $vdj $allelecutoff ${base}/${process}/substractions/*.stats
        echo "BedGraphs for ${results}/substractions"
        ${XDIR}/fileutilities.py T ${base}/${process}/substractions --dir bedGraph | srun --ntasks=1 --mpi=none ${XDIR}/fileutilities.py P --loop mv {abs} ${base}/${process}/substractions/tmp \&\& printf '"track type=bedGraph name=%s\n"' "{ali}" \> {abs} \&\& cat ${base}/${process}/substractions/tmp \>\> {abs}
        mv ${base}/${process}/substractions/*point*bedGraph ${base}/${results}/substractions/
    fi
fi

echo ""
echo "Cleaning up file-pollution"
# Trash (most) intermediate files
# if [[ "$dopre" -eq 1 ]]; then
#     rm ${base}/${process}/*trim53.fastq.gz
# fi
# if [[ "$flash" -eq 1 ]] && [[ "$doaln" -eq 1 ]]; then
#     rm ${base}/${process}/*extendedFrags*.fastq.gz
#     if [ "$allowoutties" -eq 1];  then
#         rm ${base}/${process}/*histogram.outie
#         rm ${base}/${process}/*histogram.innie
#         rm ${base}/${process}/*hist.innie
#         rm ${base}/${process}/*hist.outie
#     else
#         rm ${base}/${process}/*histogram
#         rm ${base}/${process}/*hist
#     fi
# fi
# if [[ "$dopile" -eq 1 ]]; then
#     ${XDIR}/fileutilities.py T ${base}/${process} --dir '\d+\-\d+.bam' | ${XDIR}/fileutilities.py P --loop rm {abs}
#     rm ${base}/${process}/*pileup
# fi
if [[ "$doviz" -eq 1 ]]; then
    # rm -r ${base}/${results}/*.stats_pileup_MutTypes_files
    rm -r ${base}/${results}/*coverages_files
    rm -r ${base}/${results}/*pileups_files
    if [ ! -z "$subs" ] ; then
        # rm -r ${base}/${results}/substractions/*.stats_pileup_MutTypes_files
        rm -r ${base}/${results}/substractions/*coverages_files
        rm -r ${base}/${results}/substractions/*pileups_files
    fi
fi

echo ""
echo "Finished ${run}"
