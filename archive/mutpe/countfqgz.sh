#!/usr/bin/env sh

printf '%s\t%d\n' $1 $(zcat $2 | grep -c '^@') > $3
