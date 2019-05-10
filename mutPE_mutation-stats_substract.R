#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/Rscript

# Substracts one (or more) stats files from another stats file.
# Individual mutation types at each position are summed, and the sums are then substracted across files.

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]                            # Ouput stats file
fromfile <- args[2]                           # Input TSV stats to substract from (seq \t pos \t type \t count \t depth)
statsfiles <- args[3:length(args)]            # Input TSV stats to be substracted 

library(dplyr)
library(readr)

# Left-hand side table to substract from
lhs <- read_tsv(fromfile)
if (all(names(lhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
  names(lhs) <- c('seq', 'pos', 'type', 'count', 'depth')
stopifnot(all(names(lhs) == c('seq', 'pos', 'type', 'count', 'depth')))
# Aggregate frequencies per position.
lhs <- lhs %>% 
  mutate(freq = count / depth) %>%
  group_by(seq, pos) %>%
  mutate(aggrfreq = sum(freq), type="any", depth=1) %>%  # depth 1 neutralizes re-calculating frequencies in the plotting script
  ungroup() %>%
  select(seq, pos, type, aggrfreq, depth)  %>%
  unique() %>%
  arrange(seq, pos)
names(lhs) <- c('seq', 'pos', 'type', 'freq', 'depth')


# Substract cumulatively each right-hand side table
for (s in statsfiles) {
  # Right-hand side table to be substracted in this round
  rhs <- read_tsv(s)
  if (all(names(rhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    names(rhs) <- c('seq', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(rhs) == c('seq', 'pos', 'type', 'count', 'depth')))
  # Aggregate frequencies per position.
  rhs <- rhs %>% 
    mutate(freq = count / depth) %>%
    group_by(seq, pos) %>%
    mutate(aggrfreq = sum(freq), type="any") %>%
    ungroup() %>%
    select(seq, pos, aggrfreq)  %>%
    unique() %>%
    arrange(seq, pos)
  
  # Outer join, to deal with the possibility of missing positions or missing sequences in either table
  lhs <- lhs %>% 
    full_join(rhs, by=c('seq','pos'))
  # Fix NAs from the merge, if any
  lhs$freq[is.na(lhs$freq)] <- 0
  lhs$aggrfreq[is.na(lhs$aggrfreq)] <- 0
  lhs$type[is.na(lhs$type)] <- "any"
  
  # Substract
  lhs <- lhs %>%
    mutate(freq = freq - aggrfreq) %>%
    select(seq, pos, type, freq, depth)
}

names(lhs) <- c('seq', 'pos', 'type', 'count', 'depth')
write_tsv(lhs, outfile)
