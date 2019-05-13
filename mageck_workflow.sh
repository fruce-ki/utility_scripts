#!/usr/bin/sh
#
#SBATCH --get-user-env
#SBATCH -J mageck-wf
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=mageck-wf.out
#SBATCH --error=mageck-wf.err

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i INDIR -l GUIDES_LIB -b DEMUX_TABLE -n COUNTS_OUTDIR -c CONTRASTS_TABLE -m MAGECK_OUTDIR
                [-r REVCOMP] [-p PAD_BASE] [-s SPACER_LENGTH] [-u UMI_LENGTH] [-d DEMUX_BC_LENGTH] [-M BC_MISMATCHES] [-z MINCOUNT] [-E ENTREZ_FIELD]
                [-Z REF_SAMPLES -C CTRL_GUIDE_GRPS -G GUIDES_PER_GENE]"
		exit 1
}
# Defaults
pad="C"
countthresh=50
spacer=20
umi=6
demux=4
bcmm=1
revcomp="no"
entrezfield=1   # 0-based index
# Parse options.
while getopts 'i:l:b:r:n:c:m:r:C:G:z:Z:p:s:u:d:M:E' flag; do
  case "${flag}" in
    i) indir="${OPTARG}" ;;           # Input directory with unaligned BAMs
    l) library="${OPTARG}" ;;         # sgRNA library
    b) barcodes="${OPTARG}" ;;        # Demultiplexing table: lane \t sample_name \t barcode
    n) countsdir="${OPTARG}" ;;       # Output directory for FASTQs and counts table
    c) contrasts="${OPTARG}" ;;       # Comparison details for MAGECK nextflow: name \t control \t treatment \t norm_method \t fdr_method \t lfc_method \t cnv_correction \t filter
    m) mageckdir="${OPTARG}" ;;       # Output directory for MAGECK nextflow
    r) revcomp="${OPTARG}" ;;         # Reverse complement the reads to match the barcodes? 'yes'\'no' (no)
    C) ctrlguides="${OPTARG}" ;;      # Comma-separated list of control guide group names.
    G) guidespergene="${OPTARG}" ;;   # Guides per gene
    z) countthresh="${OPTARG}" ;;     # Minimum read count per guide (50)
    Z) refSamps="${OPTARG}" ;;        # Comma-separated sample-names to apply the counts threshold for control guide purposes
    p) pad="${OPTARG}" ;;             # Padding base ("C")
    s) spacer="${OPTARG}" ;;          # Spacer length (20)
    u) umi="${OPTARG}" ;;             # UMI length (6)
    d) demux="${OPTARG}" ;;           # De-multiplexing barcode length (4)
    M) bcmm="${OPTARG}" ;;            # Barcode mismatch allowance (1)
    E) entrezfield="${OPTARG}" ;;     # 0-based index position in the undesrcore-separated composite guide IDs that represents the gene Entrez ID (1).
    *) usage ;;
  esac
done

if [ "$revcomp" = "yes" ]; then
  revcomp="--reverse_complement"
else
  revcomp=""
fi

# Check we got the minimum set
if [ -z "$indir" ] || [ -z "$library" ] || [ -z "$barcodes" ] || [ -z "$contrasts" ] || [ -z "$countsdir" ] || [ -z "$mageckdir" ]; then
  usage
  exit 1
fi


## Workflow ##

echo "Pre-process to BAM to FASTQ, align, and count."
nextflow run zuberlab/crispr-process-nf $revcomp --inputDir $indir --library $library --padding_base $pad --spacer_length $spacer --barcode_demux_length $demux --barcode_random_length $umi --barcode_demux_mismatches $bcmm --barcodes $barcodes --outputDir $countsdir -profile ii2
ld=$(basename $library)
counts="${countsdir}/counts/${ld/.txt/}/counts_mageck.txt"

# echo "Merge some columns, add some columns."
# Rscript ~/utility_scripts/sum_cols.R ./process/crispr_nf/counts/library/counts_mageck.txt ${counts/.txt/_merged.txt} 'NoIFNG_d0,NoIFNG_d0_3,NoIFNG_d0_4' NoIFNG_d0 TRUE TRUE
# $counts=${counts/.txt/_merged.txt}
# fileutilities.py T ./process/crispr_nf/counts/library/counts_d0-merged.txt ./data/counts_plasmid-pool.txt -r -i --appnd > ${counts/.txt/_added.txt}
# $counts=${counts/.txt/_added.txt}

echo ''
echo "Group control guides into genes."
if ! [ -z "$ctrlguides" ]; then
  srun Rscript ~/utility_scripts/nonTargetGuides2controlGenes.R -c $ctrlguides -f $counts -o ${counts/.txt/_ctrls-grouped.txt} -n $guidespergene -g group -t id -m $countthresh -z $refSamps
  counts="${counts/.txt/_ctrls-grouped.txt}"
else
  ctrlguides="hakunamatata_dummy" # dummy value that will not match patterns later on
fi

echo ''
echo "MAGECK."
nextflow run zuberlab/crispr-mageck-nf --contrasts $contrasts --counts $counts --outputDir $mageckdir --min_count $countthresh -profile ii2 --legacy

echo ''
echo "Rename columns."
renamed="/renamed"
srun fileutilities.py T $mageckdir --dir | fileutilities.py P --loop S sh ~/utility_scripts/mageck_rename_columns.sh {abs} {abs}${renamed}

echo ''
echo "Extract entrez IDs from sgRNA IDs."
# Also use the re-grouped group field.
lib="$(dirname $library)"
srun cut -f 1 $library | perl -e '$pat = join "|", split ",", $ARGV[0]; while($line=<STDIN>){ if ($line=~/(?<!\w)id(?!\w)/) { print "Entrez\n" } elsif ( $line=~/^($pat)/ ){ print "$1\n" } else { @a = split(/_/, $line); print "$a[$ARGV[1]]\n" }}' $ctrlguides $entrezfield > ${lib}/entrez.txt
srun fileutilities.py T $library ${lib}/entrez.txt -r --appnd outer | cut -f 1,3,4 > ${library}_entrez.txt
srun cut -f 1,2 $counts > tmp.txt
srun fileutilities.py T tmp.txt ${library}_entrez.txt -i -r --appnd outer > ${library}_forguides.txt
head -n 1 ${library}_forguides.txt | cut -f 2,4 > ${library}_forgenes.txt                         # prevent header from getting sorted to another row
srun tail -n +2 ${library}_forguides.txt | cut -f 2,4 | sort | uniq >> ${library}_forgenes.txt    # prevent header from getting sorted to another row
rm tmp.txt

echo ''
echo "Merge all gene-level outputs."
genes="${mageckdir}/genes_all.tsv"
srun --mem=10000 fileutilities.py T ${library}_forgenes.txt ${mageckdir}/*${renamed}/genes_pos_stats.txt ${mageckdir}/*${renamed}/genes_neg_stats.txt -r -i --appnd outer > $genes

echo ''
echo "Merge all guide-level outputs and clean up redundant columns."
guides="${mageckdir}/guides_all.tsv"
srun --mem=10000 fileutilities.py T ${library}_forguides.txt ${mageckdir}/*${renamed}/guides_stats.txt -r -i --appnd outer > $guides
dups=$(head -n1 $guides | fileutilities.py D --swap "\n" | perl -e '$i=0; while($field = <STDIN>){print "$i " if $field=~/group/; $i++} print "\n";')
srun --mem=10000 fileutilities.py T $guides -r --mrgdups $dups > ${guides/.tsv/_dedup.tsv}
guides="${guides/.tsv/_dedup.tsv}"
nc=$(perl -e '$ARGV[0] =~/^(\d+)/; print $1' $(fileutilities.py T $guides --cntcols))
srun --mem=10000 fileutilities.py T $guides -r --cols 0 $(expr $nc - 1) 1:$(expr $nc - 2) > ${guides/.tsv/_reord.tsv}
guides="${guides/.tsv/_reord.tsv}"

echo ''
echo "Add -log10(p)."
srun --mem=10000 add_log10p.R -i $genes -o ${genes/.tsv/_l10p.tsv} -r group
srun --mem=10000 add_log10p.R -i $guides -o ${guides/.tsv/_l10p.tsv}

echo ''
echo "Cleaning up intermediate files."
rm ${library}_entrez.txt
rm ${library}_forguides.txt
rm ${library}_forgenes.txt
rm ${mageckdir}/genes_all.tsv
rm ${mageckdir}/guides_all.tsv
rm ${mageckdir}/guides_all_dedup.tsv
rm ${mageckdir}/guides_all_dedup_reord.tsv
rm -r ${mageckdir}/*${renamed}

echo "All done!"
