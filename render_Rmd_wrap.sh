#!/usr/bin/env bash

Rscript -e "rmarkdown::render('$1')"
