#!/users/kimon.froussios/miniconda3/envs/bioinfo/bin/Rscript

library(data.table)

# Extract coordinates of natural V,D, and J

args <- commandArgs(trailingOnly = TRUE)
refseq <- args[1]
out <- args[2]

# refseq <- '/Volumes/groups/pavri/Kimon/ursi/PRO/aux/mm10_ncbi/GCF_000001635.26_GRCm38.p6_feature_table.txt'
# out <- '/Volumes/groups/pavri/Kimon/ursi/PRO/aux/mm10_ncbi/GCF_000001635.26_GRCm38.p6_VDJ.bed'

dt <- fread(refseq)

dt <- dt[chromosome=='12', ]
dt <- dt[(grepl('Igh', dt$symbol)), ]
names(dt)[1] <- 'feature'
print(paste( 'IgH:', min(dt$start), '-', max(dt$end) ))
dt <- dt[(grepl('[VDJ]_segment', dt$feature, fixed=FALSE)), ]

v1 <- dt[(grepl('V_segment', dt$feature, fixed=FALSE)) & seq_type == 'chromosome', ]
v2 <- dt[(grepl('V_segment', dt$feature, fixed=FALSE)) & seq_type != 'chromosome', ]
d <- dt[(grepl('D_segment', dt$feature, fixed=FALSE)), ]
j <- dt[(grepl('J_segment', dt$feature, fixed=FALSE)), ]
print(paste( 'IgH main V:', min(v1$start), '-', max(v1$end) ))
print(paste( 'IgH alternative V:', min(v2$start), '-', max(v2$end) ))
print(paste( 'IgH D:', min(d$start), '-', max(d$end) ))
print(paste( 'IgH J:', min(j$start), '-', max(j$end) ))

fwrite(rbindlist(list(v1,d,j))[, .(chromosome, start, end)], file=out, sep="\t", col.names=FALSE)
