library(ngram)
library(data.table)
library(stringi)


txt <- multiread('C:/Users/jack_/Downloads/utility_scripts-master/utility_scripts-master' , extension ="*")

# Count characters
ng1 <- ngram(tolower(unlist(txt)), n=1, sep='')
F1 <- as.data.table(get.phrasetable(ng1))

# Count bigrams
ng2 <- ngram(tolower(unlist(txt)), n=2, sep='')
F2 <- as.data.table(get.phrasetable(ng2))

# Focus on letters
subF1 <- F1[grepl('[a-z]', ngrams)]
subF2 <- F2[grepl('[a-z] [a-z]', ngrams)]
subF2[, fst := substr(ngrams, 1, 2)]
subF2[, snd := substr(ngrams, 3, 4)]
setorder(subF2, fst, -freq)

# Sum the mirror-pair bigrams
revF2 <- copy(subF2)
revF2[, fst := substr(ngrams, 3, 4)]
revF2[, snd := substr(ngrams, 1, 2)]
revF2[, ngrams := paste0(fst, snd)]

subF2 <- merge(subF2, revF2, by='ngrams', all.x=TRUE)
subF2[, freq := freq.x]
subF2[fst.x != snd.x, freq := rowSums(subF2[fst.x != snd.x, .(freq.x, freq.y)], na.rm=TRUE)]
subF2 <- subF2[, .(ngrams, freq, fst.x, snd.x)]

topF2 <- subF2[, head(.SD, 6), by=fst.x]

# View(subF1)
# View(subF2)
# View(topF2)

fwrite(subF2, file="C:/Users/jack_/Downloads/utility_scripts-master/ngrams.txt")
