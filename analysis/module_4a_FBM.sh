#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Combine the original variants and the knockoffs into a single BED file, 
# then memory-map it to prepare for computation of test statistics.
#
# Authors: Matteo Sesia
# Date:    11/08/2018

################
# Parse input  #
################

RESOLUTION=$1

#########
# Setup #
#########

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

echo "----------------------------------------------------------------------------------------------------"
echo "Merging knockoff-augmented BED files the following parameters:"
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "  - Resolution        :" $RESOLUTION
echo "----------------------------------------------------------------------------------------------------"
echo ""

#############################
# Merge BED files for lasso #
#############################

# Input
GENO_FILE="/scratch/PI/candes/ukbiobank_tmp/knockoffs/"$RESOLUTION"_K50/ukb_gen_chr"

# Location for temporary data storage
AUGMENTED_DATA_DIR="/scratch/PI/candes/ukbiobank_tmp/augmented_data_big"
mkdir -p $AUGMENTED_DATA_DIR

# File name for temporary data storage
AUGMENTED_DATA=$AUGMENTED_DATA_DIR"/ukb_gen_"$RESOLUTION

if [ ! -s $AUGMENTED_DATA".bed" ]; then
  # Merge multiple BED files
  GENO_FILE="/scratch/PI/candes/ukbiobank_tmp/knockoffs/"$RESOLUTION"_K50/ukb_gen_chr"
  MERGE_LIST=$AUGMENTED_DATA_DIR"/"$RESOLUTION"_mergelist.txt"
  for CHR in $CHR_LIST; do
    if [ $CHR == 1 ]; then
      rm -rf $MERGE_LIST
      touch $MERGE_LIST
    else
      echo $GENO_FILE$CHR >> $MERGE_LIST
    fi
  done

  plink \
    --bfile $GENO_FILE"1" \
    --merge-list $MERGE_LIST \
    --make-bed \
    --memory 39000 \
    --out $AUGMENTED_DATA

  rm $AUGMENTED_DATA".log"
else
  echo "--------------------------------------------------"
  echo "Skipping merge because:"
  echo " "$AUGMENTED_DATA".bed exists."
  echo "--------------------------------------------------"
fi

##################################
# Make list of original variants #
##################################

Rscript --vanilla utils/list_original.R $RESOLUTION

######################
# Convert BED to FBM #
######################

# Call R script to convert BED to FBM
if [ ! -s $LASSO_DATA".bk" ] | [ ! -s $LASSO_DATA".rds" ]; then
  echo "--------------------------------------------------"
  echo "Converting BED to FBM"
  echo "--------------------------------------------------"
  source activate ukb
  Rscript --vanilla make_FBM.R $RESOLUTION
else
  echo "--------------------------------------------------"
  echo "Skipping conversion of BED to FBM because:"
  echo " "$LASSO_DATA".bk exists."
  echo " "$LASSO_DATA".rds exists."
  echo "--------------------------------------------------"
fi
