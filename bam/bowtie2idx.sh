#!/usr/bin/env sh

# module load bowtie2/2.2.9-foss-2017a

fa=$1 ; shift
pref=$1 ; shift
params=( "$@" )

bowtie2-build $fa $pref $params
