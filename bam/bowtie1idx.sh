#!/usr/bin/env sh

#module load bowtie/1.2.2_p1-foss-2017a

fa=$1 ; shift
pref=$1 ; shift
params=( "$@" )

bowtie-build $fa $pref $params
