#!/bin/bash

################
# Parse input  #
################

PHENOTYPE=$1
RESOLUTION=$2

#########
# Setup #
#########

# Whether to run LASSO
RUN_LASSO=1

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

echo "----------------------------------------------------------------------------------------------------"
echo "Performing association analysis with knockoffs using the following parameters:"
echo "  - Phenotype         :" $PHENOTYPE
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "  - Resolution        :" $RESOLUTION
echo "  - Run-LASSO         :" $RUN_LASSO
echo "----------------------------------------------------------------------------------------------------"
echo ""

#################################
# Initialize input/output paths #
#################################

if [ $RUN_LASSO -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing importance measures with LASSO"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
  LASSO_DIR="/scratch/PI/candes/ukbiobank_tmp/analysis/knockoffs"
  mkdir -p $LASSO_DIR
  source activate ukb
  Rscript --vanilla utils/lasso.R $PHENOTYPE $RESOLUTION
else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping computation of importance measures with LASSO"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
fi
