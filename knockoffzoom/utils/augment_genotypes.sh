#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Package the knockoffs
#
# Authors: Matteo Sesia
# Date:    07/19/2018

GEN_BASENAME=$1
OUT_BASENAME=$2
RESOLUTION=$3
CHR_MIN=$4
CHR_MAX=$5

# List of chromosomes
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

# Make list of chromosomes to be merged with the first one
MERGE_LIST=$GEN_BASENAME"_mergelist.txt"
rm -f $MERGE_LIST
touch $MERGE_LIST
for CHR in ${CHR_LIST[@]}; do
  if [ $CHR != $CHR_MIN ]; then
    # Basename for the augmented genotype file for this chromosome
    CHR_BASENAME="../tmp/example_chr"$CHR"_res"$RESOLUTION
    echo $CHR_BASENAME >> $MERGE_LIST
  fi
done

# Merge the augmented data from all chromosomes
CHR_BASENAME_FIRST="../tmp/example_chr"$CHR_MIN"_res"$RESOLUTION
plink \
  --bfile $CHR_BASENAME_FIRST \
  --merge-list $MERGE_LIST \
  --make-bed \
  --memory 5000 \
  --out $OUT_BASENAME
