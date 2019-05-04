#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Compute the KnockoffZoom test statistics
#
# Authors: Matteo Sesia
# Date:    07/19/2018

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR

# Storage of phenotype files
PHENO_DIR="../data/phenotypes"

# Storage of output files
OUT_DIR="../results"
mkdir -p $OUT_DIR

# List of chromosomes
CHR_MIN=21
CHR_MAX=21
CHR_LIST=($(seq $CHR_MIN $CHR_MAX))

# List of resolutions
RESOLUTION_LIST=("2" "5" "10" "20" "50" "100")

# Utility scripts
BED_TO_FBM="Rscript --vanilla utils/make_FBM.R"
COMPUTE_STATS="Rscript --vanilla utils/lasso.R"
AUGMENT_GENOTYPES="utils/augment_genotypes.sh"

# Which operations should we perform?
FLAG_MAKE_FBM=1
FLAG_COMPUTE_STATS=1

############
# Make FBM #
############

if [[ $FLAG_MAKE_FBM == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting augmented genotypes into FBM"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Basename for genotypes
    GEN_BASENAME="../tmp/example_res"$RESOLUTION
    # Basename for output FBM
    OUT_BASENAME=$TMP_DIR"/example_res"$RESOLUTION
    # Combine genotypes and knockoffs into bed
    $AUGMENT_GENOTYPES $GEN_BASENAME $OUT_BASENAME $RESOLUTION $CHR_MIN $CHR_MAX
    # Convert augmented BED to FBM
    $BED_TO_FBM $OUT_BASENAME $OUT_BASENAME

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of augmented genotypes into FBM"
  echo "----------------------------------------------------------------------------------------------------"
fi

###########################
# Compute test statistics #
###########################

if [[ $FLAG_COMPUTE_STATS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing test statistics"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Augmented genotypes in FBM format
    FBM_FILE=$TMP_DIR"/example_res"$RESOLUTION".rds"
    # Knockoff key basename (wildcard ? for chromosome number)
    KEY_BASENAME=$TMP_DIR"/example_chr?_res"$RESOLUTION".key"
    # Phenotype file
    PHENO_FILE=$PHENO_DIR"/phenotypes.tab"
    # Phenotype name
    PHENO_NAME="y"
    # Output file
    OUT_BASENAME=$TMP_DIR"/example_res"$RESOLUTION
    # Compute test statistics
    $COMPUTE_STATS $FBM_FILE $KEY_BASENAME $PHENO_FILE $PHENO_NAME $OUT_BASENAME

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping test statistics"
  echo "----------------------------------------------------------------------------------------------------"
fi
