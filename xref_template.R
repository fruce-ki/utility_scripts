library(data.table)

# Merge all files with a certain pattern form a folder, with a single other table, individually

FF <- dir('/Volumes/groups/zuber/zubarchive/TRANSFER/kimon_to_markus/20191125_slam_mouse', full.names = TRUE)
FF <- FF[grep('tsv', FF)]

X <- fread('/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9262_slamseq/aux/GRCm38/xref.tsv', colClasses="character")
setkey(X, Entrez)
length(unique(X$Entrez)) == length(unique(X$Name))
length(unique(X$Entrez)) == length(unique(X$Id))

for (f in FF){
  DD <- fread(f, colClasses="character")
  setkey(DD, Name)
  
  EE <- merge(X, DD, by.x="Entrez", by.y="Name")
  fwrite(file=f, EE, sep="\t")
}
