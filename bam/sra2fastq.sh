#!/usr/bin/env sh

module load sra-toolkit/2.8.2-1-centos_linux64

fastq-dump --split-3 --split-spot --disable-multithreading --gzip -O $1 $2
