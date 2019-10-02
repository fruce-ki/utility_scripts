#!/usr/bin/env bash


mkdir -p $(dirname "$2")
featureCounts --primary -a "$3" -o "$2" "$1"
