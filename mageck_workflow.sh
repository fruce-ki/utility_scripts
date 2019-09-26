#!/usr/bin/env sh
#
#SBATCH --get-user-env
#SBATCH -J mageck-wf
#SBATCH --mem=50000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=mageck-wf.out
#SBATCH --error=mageck-wf.err


## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i INDIR -l GUIDES_LIB -b DEMUX_TABLE -n COUNTS_OUTDIR -c CONTRASTS_TABLE -m MAGECK_OUTDIR [-1] [-2]"
    echo "      There are many more options and I've lost track. Consult getops in the source code."
    exit 1
}
# Defaults
pad="C"
countthresh=50
spacer='TTCCAGCATAGCTCTTAAAC'
umi=6
demux=4
bcmm=0
bcoffset=-4
smm=2
variable=0
revcomp=0
do_pre=0
do_comparison=0
guideLen=20
legacy=""
mageck_branch=""
# Parse options.
while getopts 'i:l:b:B:n:c:m:r:C:G:z:Z:p:s:u:d:O:g:M:A:E:e:Vr12L' flag; do
  case "${flag}" in
    i) indir="${OPTARG}" ;;           # Input directory with unaligned BAMs
    l) library="${OPTARG}" ;;         # sgRNA library
    b) barcodes="${OPTARG}" ;;        # Demultiplexing table: lane \t sample_name \t barcode
    n) countsdir="${OPTARG}" ;;       # Output directory for FASTQs and counts table
    c) contrasts="${OPTARG}" ;;       # Comparison details for MAGECK nextflow: name \t control \t treatment \t norm_method \t fdr_method \t lfc_method \t cnv_correction \t filter
    m) mageckdir="${OPTARG}" ;;       # Output directory for MAGECK nextflow
    z) countthresh="${OPTARG}" ;;     # Minimum read count per guide (50)
    p) pad="${OPTARG}" ;;             # Padding base ("C")
    s) spacer="${OPTARG}" ;;          # Anchor sequence ('TTCCAGCATAGCTCTTAAAC')
    u) umi="${OPTARG}" ;;             # UMI length (6)
    d) demux="${OPTARG}" ;;           # De-multiplexing barcode length (4)
    O) bcoffset="${OPTARG}" ;;        # De-multiplexing barcode postition relastive to anchor (-4)
    g) guideLen="${OPTARG}" ;;        # Guide length (20). For clipping.
    M) bcmm="${OPTARG}" ;;            # Barcode mismatch allowance (0)
    A) smm="${OPTARG}" ;;             # Anchor mismatch allowance (2)
    Z) refSamps="${OPTARG}" ;;        # Comma-separated sample-names to apply the counts threshold for control guide purposes
    C) ctrlguides="${OPTARG}" ;;      # Comma-separated list of control guide group names.
    G) guidespergene="${OPTARG}" ;;   # Guides per gene.
    V) variable=1 ;;                  # Demultiplex manually, for staggered libraries (no).
    r) revcomp=1 ;;                   # Reverse complement the reads to match the barcodes? (no)
    1) do_pre=1 ;;                    # Execute demultiplexing, alignment and quantification.
    2) do_comparison=1 ;;             # Execute Mageck
    e) guidexref="${OPTARG}" ;;       # Append Entrez and other IDs from this guides-level file
    E) genexref="${OPTARG}" ;;        # Append Entrez and other IDs from this genes-level file
    L) legacy='--legacy' ;;			# Use mageck 0.5.5 instead of latest
    B) mageck_branch="-r ${OPTARG}";;			# non-master branch
    *) usage ;;
  esac
done

# Check we got the minimum set
if [ -z "$indir" ] || [ -z "$library" ] || [ -z "$barcodes" ] || [ -z "$contrasts" ] || [ -z "$countsdir" ] || [ -z "$mageckdir" ]; then
  echo "-i $indir -l $library -b $barodes -c $contrasts -n $countsdir -m $mageckdir"
  usage
  exit 1
fi

if [ $revcomp -eq 1 ]; then
  revcomp="--reverse_complement"
else
  revcomp=""
fi

## Workflow ##

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

set -e

libname=$(basename $library)
libname=${libname/.txt/}

#################### DEMULTIPLEX & COUNT

if [ $do_pre -eq 1 ]; then
  if [ $variable -eq 0 ]; then
    echo ''
    echo "Demultiplex BAM, align and count. Using the standard pipeline."
    nextflow run zuberlab/crispr-process-nf $revcomp --inputDir $indir --library $library --padding_base $pad --spacer_seq $spacer --spacer_length ${#spacer} --barcode_demux_length $demux --barcode_random_length $umi --barcode_demux_mismatches $bcmm --barcodes $barcodes --outputDir $countsdir -profile ii2
  else
    mkdir -p ${countsdir}/fastq ${countsdir}/fastqc ${countsdir}/counts/${libname} ${countsdir}/aligned/${libname}

    echo ''
    echo "Demultiplexing BAM using anchor sequence."
    module load python-levenshtein/0.12.0-foss-2017a-python-2.7.13
    module load pysam/0.14.1-foss-2017a-python-2.7.13
    fileutilities.py T ${indir}/*.bam --loop srun ,--mem=50000 ~/crispr-process-nf/bin/demultiplex_by_anchor-pos.py ,-i {abs} ,-D ${countsdir}/fastq ,-l ${countsdir}/fastq/{bas}.log ,-o $bcoffset ,-s $spacer ,-g $guideLen ,-b $barcodes ,-m $bcmm ,-M $smm ,-q 33 ,-Q \&
    module unload python-levenshtein/0.12.0-foss-2017a-python-2.7.13
    module unload pysam/0.14.1-foss-2017a-python-2.7.13
    wait_for_jobs demultip

    echo ''
    echo "FastQC (in the background)." # and don't wait for it. I don't need its output for a while.
    module load fastqc/0.11.5-java-1.8.0_121
    fileutilities.py T ${countsdir}/fastq/*/*.fqc --loop srun ,--mem-per-cpu=5000 ,--cpus-per-task 4 fastqc ,-q ,-t 4 ,-f fastq ,-o ${countsdir}/fastqc {abs} \&
    module unload fastqc/0.11.5-java-1.8.0_121

    echo ''
    echo "Compressing FASTQ (in the background)."
    fileutilities.py T ${countsdir}/fastq/*/*.fq --loop srun ,--mem=50000 gzip {abs} \&

    echo ''
    echo "Guides library to FASTA."
    cw=$(realpath $(pwd))
    cd $(dirname $library)
    srun ~/crispr-process-nf/bin/process_library.R $(basename $library) C
    cd $cw

    echo ''
    echo "Bowtie2 indexing."
    module load bowtie2/2.2.9-foss-2017a
    srun bowtie2-build ${library/.txt/.fasta} ${library/.txt/}

    echo ''
    echo "... waiting for gzip to catch up."
    wait_for_jobs gzip

    echo ''
    echo "Bowtie2 aligning."
    fileutilities.py T ${countsdir}/fastq/*/*.fq.gz --loop srun ,--mem=10000 ,--cpus-per-task=4 bowtie2 ,-x ${library/.txt/}  ,-U {abs} ,--threads 4 ,-L 20 ,--score-min 'C,0,-1' ,-N 0 ,--seed 42 '2>' ${countsdir}/aligned/${libname}/{bas}.log \> ${countsdir}/aligned/${libname}/{bas}.sam \&
    module unload bowtie2/2.2.9-foss-2017a
    wait_for_jobs bowtie2

    echo ''
    echo "Quantifying with featureCounts."
    module load subread/1.6.4-foss-2017a
    fileutilities.py T ${countsdir}/aligned/${libname}/*.sam --loop srun ,--mem-per-cpu=5000 ,--cpus-per-task=4 featureCounts ,-T 4 ,-a ${library/.txt/.saf} ,-F SAF ,-o ${countsdir}/counts/${libname}/{bas}.txt {abs} \&
    module unload subread/1.6.4-foss-2017a
    wait_for_jobs featureC

    echo ''
    echo "Combining samples into one table."
    srun --mem=5000 ~/crispr-process-nf/bin/combine_counts.R $library ${countsdir}/counts/${libname}/*.fq.txt > ${countsdir}/counts/${libname}/counts_mageck.txt
    # Fix header. Strip path, strip file extension, strip lane
    mv ${countsdir}/counts/${libname}/counts_mageck.txt ${countsdir}/counts/${libname}/_counts_mageck.txt
    head -n 1 ${countsdir}/counts/${libname}/_counts_mageck.txt | perl -e 'while(<STDIN>){~s/\S+\/(\w{9}_\d_)+(\d{8}\w?_\d{8}_)?//g;~s/\.fq//g;print}' > ${countsdir}/counts/${libname}/counts_mageck.txt
    tail -n +2 ${countsdir}/counts/${libname}/_counts_mageck.txt >> ${countsdir}/counts/${libname}/counts_mageck.txt

    echo ''
    echo "MultiQC"
    wait_for_jobs fastqc  # It should be long finished by now, but better ask.
    module load multiqc/1.3-foss-2017a-python-2.7.13
    srun multiqc -f -x *.run -o ${countsdir}/multiqc ${countsdir}/fastqc ${countsdir}/aligned/${libname} ${countsdir}/counts/${libname}
    module unload multiqc/1.3-foss-2017a-python-2.7.13

    echo ''
    echo "Cleaning up intermediate files"
    rm ${countsdir}/fastq/*/*.fqc
    rm -r ${countsdir}/fastqc
    rm ${countsdir}/counts/${libname}/*fq.txt ${countsdir}/counts/${libname}/*fq.txt.summary
  fi

  echo ''
  echo "Pre-processing finished!"
fi

#################### MAGECK

counts="${countsdir}/counts/${libname}/counts_mageck.txt"
if [[ ! -f $counts ]]; then
  exit 1 "Counts file does not exist."
fi

if [ $do_comparison -eq 1 ]; then
  echo ''
  echo "Group control guides into genes."
  if ! [ -z "$ctrlguides" ]; then
    srun mageck_nonTargetGuides2controlGenes.R -c $ctrlguides -f $counts -o ${counts/.txt/_ctrls-grouped.txt} -n $guidespergene -g group -t id -m $countthresh -z $refSamps
    counts="${counts/.txt/_ctrls-grouped.txt}"
  else
    ctrlguides="hakunamatata_dummy" # dummy value that will not match patterns later on
  fi

  echo ''
  echo "MAGECK."
  nextflow run zuberlab/crispr-mageck-nf --contrasts $contrasts --counts $counts --outputDir $mageckdir --min_count $countthresh -profile ii2 $legacy $mageck_branch

  echo ''
  echo "Rename columns."
  renamed="renamed"
  srun fileutilities.py T $mageckdir --dir | fileutilities.py P --loop sh ~/utility_scripts/mageck_rename_columns.sh {abs} {abs}/${renamed}

  ## WARNING:
  ## If appending throws ValueError about the shape not matching the index, it means that there are repeated row keys (probably in the xref files).
  ## Do NOT switch from --appnd (pandas concat) to merge (pandas merge)!!! Merge will do all the possible combinations.

  echo ''
  echo "Merge all gene-level outputs."
  genes="${mageckdir}/genes_all.tsv"
  if [ -f "$genes" ]; then
    rm $genes # clean up previous run, otherwise weird things happen
  fi
  srun --mem=10000 fileutilities.py T ${mageckdir}/*/${renamed}/genes_pos_stats.txt ${mageckdir}/*/${renamed}/genes_neg_stats.txt -r -i --appnd outer > $genes

  echo ''
  echo "Merge all guide-level outputs."
  guides="${mageckdir}/guides_all.tsv"
  if [ -f "$guides" ]; then
    rm $guides # clean up previous run, otherwise weird things happen
  fi
  srun --mem=10000 fileutilities.py T ${mageckdir}/*/${renamed}/guides_stats.txt -r -i --appnd outer > $guides

  if [ ! -z "$guidexref" ]; then
    echo ''
    echo "Add other cross-referencing ID fields to guides"
    srun --mem=10000 fileutilities.py T $guidexref $guides -r -i --appnd outer > ${guides/.tsv/_xref.tsv}
    guides="${guides/.tsv/_xref.tsv}"
  fi

  if [ ! -z "$genexref" ]; then
    echo ''
    echo "Add other cross-referencing ID fields to genes"
    srun --mem=10000 fileutilities.py T $genexref $genes -r -i --appnd outer > ${genes/.tsv/_xref.tsv}
    genes="${genes/.tsv/_xref.tsv}"
  fi

  echo ''
  echo "Deduplicate group columns."
  dups=$(head -n1 $genes | fileutilities.py D --swap "\n" | perl -e '$i=0; while($field = <STDIN>){print "$i " if $field=~/group/; $i++} print "\n";')
  srun --mem=10000 fileutilities.py T $genes -r --mrgdups $dups > ${genes/.tsv/_dedup.tsv}
  genes="${genes/.tsv/_dedup.tsv}"
  nc=$(perl -e '$ARGV[0] =~/^(\d+)/; print $1' $(fileutilities.py T $genes --cntcols))
  srun --mem=10000 fileutilities.py T $genes -r --cols 0 $(expr $nc - 1) 1:$(expr $nc - 2) > ${genes/.tsv/_reord.tsv}
  genes="${genes/.tsv/_reord.tsv}"

  dups=$(head -n1 $guides | fileutilities.py D --swap "\n" | perl -e '$i=0; while($field = <STDIN>){print "$i " if $field=~/group/; $i++} print "\n";')
  srun --mem=10000 fileutilities.py T $guides -r --mrgdups $dups > ${guides/.tsv/_dedup.tsv}
  guides="${guides/.tsv/_dedup.tsv}"
  nc=$(perl -e '$ARGV[0] =~/^(\d+)/; print $1' $(fileutilities.py T $guides --cntcols))
  srun --mem=10000 fileutilities.py T $guides -r --cols 0 $(expr $nc - 1) 1:$(expr $nc - 2) > ${guides/.tsv/_reord.tsv}
  guides="${guides/.tsv/_reord.tsv}"

  echo ''
  echo "Add -log10(p)."
  srun --mem=10000 mageck_add_log10p.R -i $genes -o ${genes/.tsv/_l10p.txt} -r group
  srun --mem=10000 mageck_add_log10p.R -i $guides -o ${guides/.tsv/_l10p.txt}
  genes="${genes/.tsv/_l10p.txt}"
  guides="${guides/.tsv/_l10p.txt}"

  echo ''
  echo "Add good guides ratio to genes."
  srun --mem=10000 mageck_add_ggratio.R -i $genes -o ${genes/.txt/_gg.txt} -f $counts -m $countthresh -z $refSamps
  genes="${genes/.txt/_ggratio.txt}"

  echo ''
  echo "Cleaning up intermediate files."
  rm ${mageckdir}/genes_all.tsv
  rm ${mageckdir}/genes_all_xref.tsv
  rm ${mageckdir}/genes_all_xref_dedup.tsv
  rm ${mageckdir}/genes_all_xref_dedup_reord.tsv
  rm ${mageckdir}/genes_all_xref_dedup_reord_l10p.txt
  rm ${mageckdir}/guides_all.tsv
  rm ${mageckdir}/guides_all_xref.tsv
  rm ${mageckdir}/guides_all_xref_dedup.tsv
  rm ${mageckdir}/guides_all_xref_dedup_reord.tsv
  rm -r ${mageckdir}/*/${renamed}

  echo ''
  echo "Comparisons finished!"
fi

if [ -d "./work" ]; then
	rm -r ./work/
	rm -r ./.nextflow*
	rm timeline*
fi

echo ''
echo "All done!"
