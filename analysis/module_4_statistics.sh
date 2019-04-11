#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Compute the KnockoffZoom test statistics
#
# Authors: Matteo Sesia
# Date:    07/19/2018

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR

# Storage of phenotype files
PHENO_DIR="../data/phenotypes"

# Storage of output files
OUT_DIR="../results"
mkdir -p $OUT_DIR

# List of chromosomes
CHR_LIST=($(seq 21 22))

# List of resolutions
RESOLUTION_LIST=("2" "5" "10" "20" "50" "100")

# Utility scripts
BED_TO_FBM="Rscript --vanilla ../utils/make_FBM.R"
COMPUTE_STATS="Rscript --vanilla ../utils/lasso.R"

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

    # Basename for output FBM
    OUT_BASENAME=$TMP_DIR"/example_res"$RESOLUTION

    # Make list of chromosomes to be merged with the first one
    GEN_BASENAME="../tmp/example_res"$RESOLUTION
    MERGE_LIST=$GEN_BASENAME"_mergelist.txt"
    rm -f $MERGE_LIST
    touch $MERGE_LIST
    for CHR in ${CHR_LIST[@]}; do
      if [ $CHR != ${CHR_LIST[0]} ]; then
        # Basename for the augmented genotype file for this chromosome    
        CHR_BASENAME="../tmp/example_chr"$CHR"_res"$RESOLUTION
        echo $CHR_BASENAME >> $MERGE_LIST
      fi
    done

    # Merge the augmented data from all chromosomes
    CHR_FIRST=${CHR_LIST[0]}
    CHR_BASENAME_FIRST="../tmp/example_chr"$CHR_FIRST"_res"$RESOLUTION
    plink \
      --bfile $CHR_BASENAME_FIRST \
      --merge-list $MERGE_LIST \
      --make-bed \
      --memory 5000 \
      --out $OUT_BASENAME

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
