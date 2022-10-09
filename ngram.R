library(ngram)
library(data.table)
library(stringi)


txt1 <- multiread('C:/Users/jack_/Downloads/utility_scripts-master/utility_scripts-master' , extension ="Rmd")
ng1 <- ngram(tolower(unlist(txt1)), n=2, sep='')
F1 <- as.data.table(get.phrasetable(ng1))


txt2 <- multiread('C:/Users/jack_/Downloads/utility_scripts-master/utility_scripts-master' , extension ="py")
ng2 <- ngram(tolower(unlist(txt2)), n=2, sep='')
F2 <- as.data.table(get.phrasetable(ng2))


txt3 <- multiread('C:/Users/jack_/Downloads/utility_scripts-master/utility_scripts-master' , extension ="*")
ng3 <- ngram(tolower(unlist(txt3)), n=2, sep='')
F3 <- as.data.table(get.phrasetable(ng3))


subF1 <- F1[grepl('[a-z] [a-z]', ngrams)]
subF2 <- F2[grepl('[a-z] [a-z]', ngrams)]
subF3 <- F3[grepl('[a-z] [a-z]', ngrams)]

setkey(subF3, ngrams)

vapply(subF3$ngrams, function(x){
  # x <- 'b a'
  sum(subF3[paste0(c(x, stri_reverse(x)),' '), freq], na.rm=TRUE)
}, numeric(1))



expand.grid(letters, letters)
