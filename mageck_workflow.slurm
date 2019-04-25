#!/usr/bin/sh
#
#SBATCH --get-user-env
#SBATCH --mem=10000
#SBATCH -J mageck-wf
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --output=mageck-wf.out
#SBATCH --error=mageck-wf.err

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i INDIR -l GUIDES_LIB -b DEMUX_TABLE -n COUNTS_OUTDIR -c CONTRASTS_TABLE -m MAGECK_OUTDIR
                [-r REVCOMP] [-p PAD_BASE] [-s SPACER_LENGTH] [-u UMI_LENGTH] [-d DEMUX_BC_LENGTH] [-z MINCOUNT] [-E ENTREZ_FIELD]
                [-Z SAMPLES -C CTRL_GUIDE_GRPS -G GUIDES_PER_GENE]"
		exit 1
}
# Parse options.
while getopts 'i:l:b:r:n:c:m:r:C:G:z:Z:p:s:u:d:E' flag; do
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
    Z) countthresh="${OPTARG}" ;;     # Comma-separated sample-names to apply the counts threshold for control guide purposes
    p) pad="${OPTARG}" ;;             # Padding base ("C")
    s) spacer="${OPTARG}" ;;          # Spacer length (20)
    u) umi="${OPTARG}" ;;             # UMI length (6)
    d) demux="${OPTARG}" ;;           # De-multiplexing barcode length (4)
    E) entrezfield="${OPTARG}" ;;     # 0-based index position in the undesrcore-separated composite guide IDs that represents the gene Entrez ID (2).
    *) usage ;;
  esac
done
# Check we got the minimum set
if [[ -z "$indir" || -z "$library" || -z "$barcodes" || -z "$revcomp" || -z "$contrasts" || -z "$countsdir" || -z "$mageckdir" ]]; then
  usage
  exit 1
fi

# Defaults
if [ "$revcomp" == "yes" ]; then
  revcomp="--reverse_complement"
else
  revcomp=""
fi

if [ -z "$pad" ]; then
  pad="C"
fi

if [ -z "$countthresh" ]; then
  countthresh=50
fi

if [ -z "$spacer" ]; then
  spacer=20
fi

if [ -z "$umi" ]; then
  umi=6
fi
if [ -z "$demux" ]; then
  demux=4
fi

if [ -z "$entrezfield" ]; then
  entrezfield=2
fi

## Workflow ##
empty=""


echo "Pre-process to BAM to FASTQ, align, and count."
nextflow run zuberlab/crispr-process-nf $revcomp --inputDir $indir --library $library --padding_base $pad --spacer_length $spacer --barcode_demux_length $demux --barcode_random_length $umi --barcodes $barcodes --outputDir $countsdir -profile ii2 --resume
ld=$(basename $library)
$counts="${countsdir}/counts/${ld/.txt/$empty}/counts_mageck.txt"

# echo "Merge some columns, add some columns."
# Rscript ~/utility_scripts/sum_cols.R ./process/crispr_nf/counts/library/counts_mageck.txt ${counts/.txt/_merged.txt} 'NoIFNG_d0,NoIFNG_d0_3,NoIFNG_d0_4' NoIFNG_d0 TRUE TRUE
# $counts=${counts/.txt/_merged.txt}
# fileutilities.py T ./process/crispr_nf/counts/library/counts_d0-merged.txt ./data/counts_plasmid-pool.txt -r -i --appnd > ${counts/.txt/_added.txt}
# $counts=${counts/.txt/_added.txt}

echo "Group control guides into genes."
if ! [ -z "$ctrlguides" ]; then
  Rscript ~/utility_scripts/nonTargetGuides2controlGenes.R -c $ctrlguides -f $counts -o ${counts/.txt/_ctrls-grouped.txt} -n $guidespergene -g group -t id -m $countthresh -z 'NoIFNG_d0,imc_surface_plasmid'
  $counts="${counts/.txt/_ctrls-grouped.txt}"
else
  $ctrlguides="hakunamatata_dummy" # dummy value that will not match patterns later on
fi

echo "Extract entrez IDs from sgRNA IDs."
$lib="$(dirname $library)"
cut -f 1 $library | perl -e '$pat = join "|", split "," $ARGV[0]; while($line=<STDIN>){ if ($line=~/(?<!\w)id(?!\w)/) { print "Entrez\n" } elsif ( $line=~/^($pat)/ ){ print "$1\n" } else { @a = split(/_/, $line); print "$a[$ARGV[2]]\n" }}' $ctrlguides $entrezfield > ${lib}/entrez.txt
fileutilities.py T $library ${lib}/entrez.txt --appnd outer -r > ${library}_entrez.txt
library=${library}_entrez.txt

echo "MAGECK."
nextflow run zuberlab/crispr-mageck-nf --contrasts $contrasts --counts $counts --outputDir $mageckdir --min_count $countthresh -profile ii2 --legacy --resume

echo "Rename columns."
renamed="/renamed"
fileutilities.py T $mageckdir --dir | fileutilities.py P --loop S sh ~/utility_scripts/rename_mageck_columns.sh {abs} {abs}${renamed}

echo "Merge all gene-level outputs."
genes="${mageckdir}/genes_all.tsv"
fileutilities.py T ${mageckdir}/*${renamed}/genes_pos_stats.txt ${mageckdir}/*${renamed}/genes_neg_stats.txt -r -i --appnd outer > $genes

echo "Merge all guide-level outputs and clean up redundant columns."
guides="${mageckdir}/guides_all.tsv"
fileutilities.py T $library ${mageckdir}/*${renamed}/guides_stats.txt -r -i --appnd outer > $guides
fileutilities.py T $guides -r --mrgdups $(head -n1 $guides | fileutilities.py D --swap "\n" | perl -e '$i=0; while($field = <STDIN>){print "$i " if $field=~/group/; $i++} print "\n";') > ${guides/.tsv/_dedup.tsv}
guides="${guides/.tsv/_dedup.tsv}"
nc=$(fileutilities.py T $guides --cntcols)
fileutilities.py T $guides -r --cols $(expr $nc - 1) $(expr $nc - 2) 0 $(expr $nc - 3) 1:$(expr $nc - 4) > ${guides/.tsv/_reord.tsv}
guides="${guides/.tsv/_reord.tsv}"

echo "Merge guide-level and gene-level outputs and clean up redundant columns."
$out="${mageckdir}/guides+genes.tsv"
fileutilities.py T $guides $genes -r --merge outer yes no > $out
fileutilities.py T $out -r --mrgdups $(head -n1 $out | fileutilities.py D --swap "\n" | perl -e '$i=0; while($field = <STDIN>){print "$i " if $field=~/group/; $i++}') > ${out/.tsv/_dedup.tsv}
out="${out/.tsv/_dedup.tsv}"
nc=$(fileutilities.py T $out --cntcols)
fileutilities.py T $out -r --cols $(expr $nc - 1) 0:$(expr $nc - 2) > ${out/.tsv/_reord.tsv}
out="${out/.tsv/_reord.tsv}"

echo "Add -log10(p)."
l10="${mageckdir}/tmp_log10p.tsv"
Rscript ~/utility_scripts/p2log10p.R $out $l10
fileutilities T $out $l10 -i -r --appnd outer > ${out/.tsv/_l10p.tsv}
#rm l10

echo "All done!"
