#!/users/kimon.froussios/miniconda3/envs/bioinfo/bin/Rscript

################################################################################
# process library text file
# part of CRISPR / shRNA screen pre-processing pipeline
#
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2017/08/30
#
# Modified: Kimon Froussios
# IMP, Vienna, Austria
# 2019/04/02
################################################################################

library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(Biostrings)
library(purrr)

# library file should be tab-separated text file with three columns:
# 1) id      : unique id of sgRNA / shRNA
# 2) group   : group targeted by sgRNA / shRNA (e.g gene id or domain)
# 3) sequence: sequence of sgRNA / shRNA as it appears in the sequencing reads

### command line parameters
args         <- commandArgs(trailingOnly = TRUE)
padding_base <- toupper(args[1])
input_file   <- args[2]
output_file <- args[3]
#library_name <- stringr::str_replace(basename(input_file), ".txt", "")

### import
raw <- read_tsv(input_file)

### check for id and sequence duplication
stopifnot(!any(duplicated(raw$id)))
stopifnot(!any(duplicated(raw$sequence)))

### generate fasta file for bowtie2 index
seq_length <- max(nchar(raw$sequence))

raw$sequence %>%
  toupper %>%
  str_pad(pad = padding_base, width = seq_length, side = "left") %>%
  set_names(raw$id) %>%
  DNAStringSet() %>%
#  writeXStringSet(paste0(library_name, ".fasta"), format = "fasta")
  writeXStringSet(output_file, format = "fasta")

### generate SAF annotation file for featureCount (subread package)
# tibble::tibble(GeneID = raw$id,
#               Chr    = raw$id,
#               Start  = 1,
#               End    = seq_length,
#               Strand = "*") %>%
#  readr::write_tsv(paste0(library_name, ".saf"))
