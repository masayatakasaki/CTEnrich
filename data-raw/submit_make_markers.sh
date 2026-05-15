#!/bin/bash
#$ -N make_markers
#$ -l mfree=4G
#$ -cwd
#$ -o data-raw/make_markers.log
#$ -e data-raw/make_markers.err

module load R/4.5.1
Rscript data-raw/make_markers.R
