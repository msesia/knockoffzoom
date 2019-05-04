#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Generate knockoff negative-controls
#
# Authors: Matteo Sesia
# Date:    07/19/2018

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR

# Storage of output files
OUT_DIR="../results"
mkdir -p $OUT_DIR

# List of chromosomes
CHR_LIST=$(seq 21 22)

# List of resolutions
RESOLUTION_LIST=("2" "5" "10" "20" "50" "100")

# Utility scripts
GENERATE_KNOCKOFFS="Rscript --vanilla utils/knockoffs.R"
PACKAGE_KNOCKOFFS="utils/package_knockoffs.sh"

# Which operations should we perform?
FLAG_GENERATE_KNOCKOFFS=1

######################
# Generate knockoffs #
######################

if [[ $FLAG_GENERATE_KNOCKOFFS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Generating knockoffs"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do
    for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing chromosome "$CHR" at resolution "$RESOLUTION" ..."
    echo ""

    # Basename for the HMM files produced by fastPHASE
    HMM_BASENAME="../tmp/example_chr"$CHR
    # Basename for the true HMM used to generate this data
    #HMM_BASENAME="../data/hmm/example_chr"$CHR

    # Basename for the INP haplotype files
    INP_BASENAME="../tmp/example_chr"$CHR

    # Partition file
    GROUPS_FILE=$TMP_DIR"/example_chr"$CHR"_groups"$RESOLUTION".txt"

    # Basename for output knockoff PED files
    OUT_BASENAME=$TMP_DIR"/example_chr"$CHR"_res"$RESOLUTION

    # Generate knockoffs and save the augmented genetic data in PED format
    $GENERATE_KNOCKOFFS $HMM_BASENAME $INP_BASENAME $GROUPS_FILE $OUT_BASENAME

    # Sample information file
    SAMPLE_FILE=$TMP_DIR"/example_chr"$CHR".sample"

    # Convert knockoffs to BED format
    $PACKAGE_KNOCKOFFS $OUT_BASENAME $SAMPLE_FILE

    done
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping generation of knockoffs"
  echo "----------------------------------------------------------------------------------------------------"
fi
