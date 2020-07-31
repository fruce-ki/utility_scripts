#!/usr/bin/env sh

printf '%s\t%d\n' $1 $(grep -c '^@' $2) > $3
