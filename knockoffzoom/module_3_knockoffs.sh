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
RESOLUTION_LIST=("2" "5" "10" "20" "50")
#RESOLUTION_LIST=("5")

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

    # Convert PED to BED
    plink \
      --file $OUT_BASENAME \
      --no-fid --no-parents --no-sex --no-pheno --keep-allele-order \
      --make-bed \
      --out $OUT_BASENAME

    # Add subject information to FAM file
    SEX_FILE="../tmp/example_chr"$CHR".sample"
    awk 'FNR > 2 { print $1,$4 }' $SEX_FILE > $OUT_BASENAME".sex"
    awk '{ print $1,$2,$3,$4,$6 }' $OUT_BASENAME".fam"> $OUT_BASENAME".fam.tmp"
    join -1 1 -2 1 $OUT_BASENAME".fam.tmp" $OUT_BASENAME".sex" > $OUT_BASENAME".fam.tmp2"
    awk '{ print $1,$2,$3,$4,$6,$5 }' $OUT_BASENAME".fam.tmp2"> $OUT_BASENAME".fam"
    rm $OUT_BASENAME".fam.tmp"
    rm $OUT_BASENAME".fam.tmp2"
    rm $OUT_BASENAME".sex"

    # Remove PED file
    rm $OUT_BASENAME".ped"

    # Convert knockoffs to BED format
    #$PACKAGE_KNOCKOFFS $OUT_BASENAME $SAMPLE_FILE

    done
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping generation of knockoffs"
  echo "----------------------------------------------------------------------------------------------------"
fi
