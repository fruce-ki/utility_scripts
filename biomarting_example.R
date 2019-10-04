library(biomaRt)
library(data.table)

# DT1 <- fread('/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/results/condition_all.tsv')
# length(DT1$row_ID) == length(unique(DT1$row_ID))
# setkey(DT1, row_ID)
# 
# DT2 <- fread('/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/process/count/counts_out_all.txt')
# length(DT2$row_ID) == length(unique(DT2$row_ID))
# setkey(DT2, row_ID)

listMarts()
ensmm = useMart("ENSEMBL_MART_MOUSE")
listDatasets(ensmm)
ensmm = useDataset("mc57bl6nj_gene_ensembl", mart=ensmm)
attributes = listAttributes(ensmm)
attributes[grepl('[Ee]ntrez', attributes$name, perl=TRUE),]

wanted = c("external_gene_name", "entrezgene")

xref <- getBM(attributes=wanted, mart = ensmm)
fwrite(xref, file="/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/aux/xref_entrez.txt")


