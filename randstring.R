#!/usr/bin/env Rscript

# LETTERS is a built-in
# letters is a built-in
numbers <- as.character(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
symbols <- c('~', '`', '!', '@', '#', '$', '%', '^', '&', '*', '(', ')', '_', '-', '+', '=', '{', '[', '}', ']', '|', '\\', ':', ';', '"', '<', ',', '>', '.', '?', '/')

library(getopt)

spec <- matrix(c(
  'length'      , 'l', 1, "numeric", "password length",
  'include'     , 'i', 1, "character", "something to include",
  'incsymb'     , 's', 0, "logical", "include symbols",
  'seed'        , 'r', 1, "numeric", "seed the randomizer"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

# opt <- list(length=18, include='foobar', incsymb=TRUE)

if (!is.null(opt$seed))
  set.seed(opt$seed) 


# Cover the character rules
if (opt$incsymb) {
  pool <- c(LETTERS, letters, numbers, symbols)
  oneofeach <- c(sample(LETTERS, 1), sample(letters, 1), sample(numbers, 1), sample(symbols, 1))
  length <- opt$length - 4
} else {
  pool <- c(LETTERS, letters, numbers)
  oneofeach <- c(sample(LETTERS, 1), sample(letters, 1), sample(numbers, 1))
  length <- opt$length - 3
}

# fill up rest of required length
fillup <- sample(pool, length, replace = TRUE)


# Combine and mix the two sets
passvec <- sample(c(oneofeach, fillup), opt$length, replace = FALSE)


# include thing at random location
loc <- sample(1:(opt$length + 1), 1)
if (loc == 1){
  passvec <- c(opt$include, passvec)
} else if (loc > opt$length) {
  passvec <- c(passvec, opt$include)
} else {
  passvec <- c(passvec[1:(loc-1)], opt$include, passvec[loc:opt$length])
}

# assemble
passw <- paste0(passvec, collapse='')

message(passw)
