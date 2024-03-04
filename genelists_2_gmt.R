library(data.table)

baseDir <- '/BIOINFORMATICS/REFERENCES/custom_genesets/'
genelists <- dir(baseDir)
genelists <- genelists[grepl('.genelist.txt$', genelists)]


rows <- vapply(genelists, function(f) {
  # f <- genelists[1]
  d <- fread(file.path(baseDir, f), header = FALSE)

  prefix <- sub('.genelist.txt$', '', f)
  phenotype = sub('\\..*', '', prefix)
  source = sub('.*\\.', '', prefix)

  paste(prefix, source, paste(d$V1, collapse = "\t"), sep = '\t')
}, character(1))


write(rows, file = '/BIOINFORMATICS/REFERENCES/custom_genesets/nlt_custom.gmt', ncolumns = 1)
