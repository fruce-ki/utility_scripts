#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
	'alpha',          'a', 1, 'numeric',   "Significance level.",
	'de_file',        'd', 1, 'character', "DE results. At least a name column and a value column.",
	'fx_name',        'e', 1, 'character', "Name of value column by which to rank genes (log2FoldChange).",
	'gsdb_file',      'g', 1, 'character', "Gene set database GMT file.",
	'out_dir',        'o', 1, 'character', "Destination directory for output (./results).",
	'p_name',         'p', 1, 'character', "Name of p-value column, if applicable, for gene selection.",
	'p_cut',          'P', 1, 'character', "P-value cutoff (0.05) for gene selection.",
	'reportTemplate', 'T', 1, 'character', "Rmd report template file (~/utility_scripts/GSEA_report_template.Rmd).",
	'xorg_file',      'x', 1, 'character', "CHIP file for gene identifier conversion, if necessary."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(gsdb_file="/Volumes/groups/busslinger/Bioinf/src/cbe/broad/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.mf.v7.5.1.symbols.gmt", xorg_file="/Volumes/groups/busslinger/Bioinf/src/cbe/broad/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.5.1.chip", de_file='/Volumes/groups/busslinger/Kimon/tanja/R13870_RNAseq_timecourse/results/DE/exon_genecounts.CellType_actB.Condition_Ikaros_00_5phIAA_vs_Aiolos_00_5phIAA.deseq2.nolab.spotfire.txt', alpha='0.05', out_dir='/Volumes/groups/busslinger/Kimon/tanja/R13870_RNAseq_timecourse/results/GSEA', fx_name='log2FoldChange.shrink', p_cut=0.05, p_name='padj')

stopifnot (!is.null(opt$de_file))
stopifnot (!is.null(opt$gsdb_file))
if (is.null(opt$out_dir)) opt$out_dir <- "./results"
if (is.null(opt$alpha)) opt$alpha <- 0.05
if (is.null(opt$p_cut)) opt$p_cut <- 0.05
if (is.null(opt$fx_name)) opt$fx_name <- 'log2FoldChange'
if (is.null(opt$reportTemplate)) opt$reportTemplate <- "~/utility_scripts/GSEA_report_template.Rmd"

# dir.create(opt$out_dir, recursive=TRUE)

name <- gsub('nolab\\.|spotfire\\.|_gmt', '', sub('tsv$|txt$|csv$', paste('gsea', gsub('\\.', '_', basename(opt$gsdb_file)), 'html', sep='.'), basename(opt$de_file)))

rmarkdown::render(opt$reportTemplate,
									output_file = name,
									output_dir = opt$out_dir,
									params=list(selfname = name,
															resultsDir = opt$out_dir,
															alpha = opt$alpha,
															gsdb = opt$gsdb_file,
															xorg = opt$xorg_file,
															de = opt$de_file,
															fx = opt$fx_name,
															p = opt$p_name,
															pcut= opt$p_cut
															)
)
