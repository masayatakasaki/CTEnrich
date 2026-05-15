#!/bin/bash
#$ -N make_markers_human
#$ -l mfree=4G
#$ -o make_markers_human.log
#$ -e make_markers_human.err
#$ -cwd

module load R/4.5.1
Rscript /net/hamazaki/vol1/home/masaki03/R/CTEnrich/data-raw/make_markers_human.R
