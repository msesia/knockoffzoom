#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Package the knockoffs
#
# Authors: Matteo Sesia
# Date:    07/19/2018

OUT_BASENAME=$1
SAMPLE_FILE=$2

# Convert PED to BED
plink \
  --file $OUT_BASENAME \
  --no-fid --no-parents --no-sex --no-pheno --keep-allele-order \
  --make-bed \
  --out $OUT_BASENAME

# Add subject information to FAM file
awk 'FNR > 2 { print $1,$4 }' $SAMPLE_FILE > $OUT_BASENAME".sex"
awk '{ print $1,$2,$3,$4,$6 }' $OUT_BASENAME".fam"> $OUT_BASENAME".fam.tmp"
join -1 1 -2 1 $OUT_BASENAME".fam.tmp" $OUT_BASENAME".sex" > $OUT_BASENAME".fam.tmp2"
awk '{ print $1,$2,$3,$4,$6,$5 }' $OUT_BASENAME".fam.tmp2"> $OUT_BASENAME".fam"
rm $OUT_BASENAME".fam.tmp"
rm $OUT_BASENAME".fam.tmp2"
rm $OUT_BASENAME".sex"

# Remove PED file
rm $OUT_BASENAME".ped"
