library(ngram)
library(data.table)
library(stringi)


txt <- multiread('/Users/kimon.froussios/Desktop/ngram' , extension ="*")

# Count characters
ng1 <- ngram(toupper(unlist(txt)), n=1, sep='')
F1 <- as.data.table(get.phrasetable(ng1))

# Count bigrams
ng2 <- ngram(toupper(unlist(txt)), n=2, sep='')
F2 <- as.data.table(get.phrasetable(ng2))

# Focus on letters
subF1 <- F1[grepl('[A-Z]', ngrams)]
subF2 <- F2[grepl('[A-Z] [A-Z] ', ngrams)]
subF2[, fst := substr(ngrams, 1, 1)]
subF2[, snd := substr(ngrams, 3, 3)]
setorder(subF2, fst, -freq)

# Sum the mirror-pair bigrams
subF2 <- subF2[fst != snd, ] # double letters do not create awkward movements 
subF2[, rev := paste0(snd, ' ', fst, ' ')]
subF2 <- merge(subF2, subF2[, .(rev, freq)], by.x='ngrams', by.y='rev', all.x=TRUE, suffixes=c('', '.y'))
subF2[, total := rowSums(subF2[, .(freq, freq.y)], na.rm=TRUE)]

# Proportion per letter
subF2[, frac := total/sum(.SD$total), by=fst]
setorder(subF2, fst, -total)





# View(subF1)
# View(subF2)
# View(topF2)

fwrite(subF2, file="/Users/kimon.froussios/Desktop/ngrams.txt", sep="\t", quote=FALSE)
