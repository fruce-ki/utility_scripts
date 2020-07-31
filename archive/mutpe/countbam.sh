#!/usr/bin/env sh

printf '%s\t%d\n' $1 $(samtools view -c -F 4 $2) > $3
