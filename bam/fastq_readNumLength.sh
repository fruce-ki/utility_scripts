#!/usr/bin/env sh

# Number and length of reads in UNZIPPED FastQ.

awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}'
