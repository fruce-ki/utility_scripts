library(biomaRt)
library(data.table)

DT1 <- fread('/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9179_quantseq/process/counts/all_counts.tsv')
length(DT1$Geneid) == length(unique(DT1$Geneid))
setkey(DT1, Geneid)

listMarts()
# ens = useMart("ENSEMBL_MART_MOUSE")
ens = useMart("ENSEMBL_MART_ENSEMBL")

listDatasets(ens)
# ens = useDataset("mc57bl6nj_gene_ensembl", mart=ens)
# ens = useDataset("mmusculus_gene_ensembl", mart=ens)

ens = useDataset("hsapiens_gene_ensembl", mart=ens)
attributes = listAttributes(ens)
attributes[grepl('[Ee]ntrez', attributes$name, perl=TRUE),]

wanted = c("ensembl_gene_id", "external_gene_name", "entrezgene_id")
xref <- getBM(attributes=wanted, mart = ens)

# fwrite(xref, file="/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/aux/xref_entrez.txt")
fwrite(xref, file="/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9179_quantseq/aux/xref_hs_genes.txt", sep="\t", col.names = TRUE)


