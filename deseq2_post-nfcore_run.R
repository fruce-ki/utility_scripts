#!/usr/bin/env Rscript

###
# Not intended to be given to nf-core/differentialabundance directly to replace the default template.
# Instead, meant to be run aseparately after completion of the pipeline. It ralies on the outputs of the pipeline
###

library(getopt)

spec = matrix(c(
  'baseDir',        'b', 1, "character", "Full path to base directory of the nf-core/differentialabundance results. All other paths relative from here.",
  'countsFile',     'f', 1, "character", "Salmon-formatted gene counts file.",
  'tpmFile',        'F', 1, "character", "Salmon-formatted gene TPMs file.",
  'contrasts',      'c', 1, "character", "File with the contrasts for nf-core/differentialabundance.",
  'genelists',       'g', 1, "character", "Comma-separated files with (vertical) list of genes to highlight.",
  'pathways',       'G', 1, "character", "Comma-separated (significant) pathways from which to plot leading genes.",
  'minCount',       'm', 1, "integer",   "Minimum required mean count in at least one condition (50).",
  'minTPM',         'M', 1, "numeric",   "Minimum required mean TPM in at least one single sample (5).",
  'minLFC',         'l', 1, "numeric",   "Minimum required log2 fold-change (1).",
  'maxQ',           'p', 1, "numeric",   "FDR level (0.05).",
  'nhit',           'n', 1, "integer",   "Number of hits to report (50).",
  'pad',            'z', 1, "logical",   "Disable padding of zeros.",
  'suffix',         'x', 1, "character", "A custom string to add to the report file name.",
  'reportTemplate', 'T', 1, "character", "Full path to template Rmd file (~/utility_scripts/deseq2_post-nfcore_report_template.Rmd).",
#  'useAll',         'a', 0, "logical",   "Do NOT include all the samples, include only the ones marked by 'use' (FALSE).", ## Not meaningful, because deseq2 is already finished before this script
  'help',           'h', 0, "logical",   "Help."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(baseDir = '/SCRATCH/diffex/tess',
#            countsFile = '../master.gene_counts.tsv',
#            tpmFile = '../master.gene_tpm.tsv',
#            samplesFile = '../samplesheet_tess.differentialabundance.csv',
#            contrasts = '../contrasts_tess.differentialabundance.csv',
#            pathways = 'WP_CHOLESTEROL_METABOLISM_WITH_BLOCH_AND_KANDUSCH_RUSSELL_PATHWAYS,WP_CHOLESTEROL_SYNTHESIS_DISORDERS,WP_CHOLESTEROL_BIOSYNTHESIS_PATHWAY'
#            )

# opt <- list(baseDir = 'SOURCE/sandbox/linz',
#            countsFile = '../master.gene_counts.tsv',
#            tpmFile = '../master.gene_tpm.tsv',
#            samplesFile = '../samplesheet_tess.differentialabundance.csv',
#            contrasts = '../contrasts_linz.differentialabundance.csv',
#            reportTemplate = 'SOURCE/sandbox/deseq2_post-nfcore_report_template.Rmd'
#            )

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


if (is.null(opt$minCount))
  opt$minCount <- 50
if (is.null(opt$minTPM))
  opt$minTPM <- 5
if (is.null(opt$minLFC))
  opt$minLFC <- 1
if (is.null(opt$maxQ))
  opt$maxQ <- 0.05
if (is.null(opt$nhit))
  opt$nhit <- 50
if (is.null(opt$reportTemplate))
  opt$reportTemplate <- "~/utility_scripts/deseq2_post-nfcore_report_template.Rmd"
if (is.null(opt$pad))
  opt$pad <- FALSE
# if (is.null(opt$allSamples))
  # opt$useAll <- FALSE


genesets <- list('/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.symbols.gmt',
                # '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c5.all.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c5.hpo.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.symbols.gmt',
                # '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c2.all.v2023.2.Hs.symbols.reduced.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c2.cp.pid.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c2.cp.biocarta.v2023.2.Hs.symbols.gmt',
                # '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c3.all.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c3.tft.gtrd.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c3.tft.tft_legacy.v2023.2.Hs.symbols.gmt',
                '/SCRATCH/REFERENCES/msigdb_v2023.2.Hs_GMTs/c1.all.v2023.2.Hs.symbols.gmt')

# if(!is.null(opt$genelists)) {
#   genelists <- strsplit(opt$genelists, '\\s*,\\s*')[[1]]
# }
opt$genelists <- list('/DROPBOX/Neurolentech Dropbox/NEUROLENTECH/BIOINFORMATICS/REFERENCES/custom_genesets/autism.sfari_long.genelist.txt',
                  '/DROPBOX/Neurolentech Dropbox/NEUROLENTECH/BIOINFORMATICS/REFERENCES/custom_genesets/epilepsy.all_uniq.genelist.txt')

if(!is.null(opt$pathways)) {
  opt$pathways <- strsplit(opt$pathways, '\\s*,\\s*')[[1]]
}


dir.create(file.path(opt$baseDir, 'custom_reports', 'reports'), recursive = TRUE)
dir.create(file.path(opt$baseDir, 'custom_reports', 'plots', 'differential'), recursive = TRUE)
dir.create(file.path(opt$baseDir, 'custom_reports', 'plots', 'gsea'), recursive = TRUE)
dir.create(file.path(opt$baseDir, 'custom_reports', 'tables', 'differential'), recursive = TRUE)
dir.create(file.path(opt$baseDir, 'custom_reports', 'tables', 'gsea'), recursive = TRUE)


contrasts <- read.table(file.path(opt$baseDir, opt$contrasts), sep = ',', header = TRUE)

for (comparison in contrasts$id) {
  # comparison <- contrasts$id[7]

  if( ! file.exists(file.path(opt$baseDir, 'other', 'deseq2', paste0(comparison, '.dds.rld.rds'))) ) {
    print(comparison, ' does not have a DESeq2 results RDS!')
    warning(comparison, ' does not have a DESeq2 results RDS!')
    next
  }

  print(comparison)
  # print(opt)

  rmarkdown::render(opt$reportTemplate,
                    output_file = paste0(comparison, opt$suffix, '.deseq2_fgsea.report.html'),
                    output_dir = file.path(opt$baseDir, 'custom_reports', 'reports'),
                    params=list(comparison = comparison,
                                da_path = opt$baseDir,
                                contrasts = opt$contrasts,
                                counts = opt$countsFile,
                                tpm = opt$tpmFile,
                                highlight_n = opt$nhit,
                                min_count = opt$minCount,
                                min_tpm = opt$minTPM,
                                min_lfc = opt$minLFC,
                                max_q = opt$maxQ,
                                genelists = opt$genelists,
                                pathways = opt$pathways,
                                pad = !opt$pad,
                                # useAll = !opt$useAll,
                                genesets = genesets )
                    )
}
