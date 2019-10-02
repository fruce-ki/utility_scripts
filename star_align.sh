#!/usr/bin/env bash

mkdir -p "$3"
STAR --runThreadN "$4" --genomeDir ./aux/star_mm10 --readFilesIn "$1" "$2" --readFilesCommand zcat --outFileNamePrefix "$3" --outSAMtype BAM SortedByCoordinate
