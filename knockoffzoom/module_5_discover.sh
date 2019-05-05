#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Threshold the the KnockoffZoom test statistics with the knockoff filter
# and report discoveries
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
CHR_LIST=($(seq 21 22))

# List of resolutions
RESOLUTION_LIST=("2" "5" "10" "20" "50")
#RESOLUTION_LIST=("5")

# Utility scripts
FILTER_STATS="Rscript --vanilla utils/filter_stats.R"

# Which operations should we perform?
FLAG_FILTER=1

#################################
# Filter the test statistics    #
#################################

if [[ $FLAG_FILTER == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Filtering the test statistics"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Stats basename
    STATS_BASENAME=$TMP_DIR"/example_res"$RESOLUTION

    # Partition file
    GROUPS_BASENAME=$TMP_DIR"/example_chr?_groups"$RESOLUTION".txt"

    # Basename for output files
    OUT_BASENAME=$OUT_DIR"/example_res"$RESOLUTION

    # Threshold the test statistics and report discoveries
    $FILTER_STATS $STATS_BASENAME $GROUPS_BASENAME $OUT_BASENAME

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping filtering of the test statistics"
  echo "----------------------------------------------------------------------------------------------------"
fi
