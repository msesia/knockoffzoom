#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Prepare phased haplotype for use with fastPHASE
#
# Authors: Matteo Sesia
# Date:    07/19/2018

# Utility scripts
BGEN_TO_HAPST="../utils/bgen_to_hapst.sh"
HAPST_TO_INP="../utils/hapst_to_inp.sh"

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR

# List of chromosomes
CHR_LIST=$(seq 21 22)

# Which operations should we perform?
CONVERT_HAPLOTYPES_HAP=0
CONVERT_HAPLOTYPES_INP=1

##########################################
# Convert BGEN v1.2 into transposed HAPS #
##########################################

if [[ $CONVERT_HAPLOTYPES_HAP == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes into HAPS.T format"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Basename for the input haplotype files (BGEN format)
    BGEN_BASENAME="../data/haplotypes/example_chr"$CHR

    # List of individuals that passed QC
    BGEN_SAMPLE="../data/qc/samples_qc.txt"

    # List of variants that passed QC
    BGEN_VARIANT="../data/qc/variants_qc.txt"

    # Basename for output haplotype files (transposed HAPS format)
    HAP_BASENAME=$TMP_DIR"/example_chr"$CHR

    # Convert BGEN haplotypes to transposed HAPS file
    $BGEN_TO_HAPST -i $BGEN_BASENAME -c $CHR -n $BGEN_SAMPLE -v $BGEN_VARIANT -b 500 -o $HAP_BASENAME

  done
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping coversion of haplotypes into HAPS.T format"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Convert HAPS.T into INP
########################################
if [[ $CONVERT_HAPLOTYPES_INP == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes into INP format"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Basename for output haplotype files (transposed HAPS format)
    HAP_BASENAME=$TMP_DIR"/example_chr"$CHR

    # Convert haplotypes
    $HAPST_TO_INP -i $HAP_BASENAME -o $HAP_BASENAME

  done
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes into INP format"
  echo "----------------------------------------------------------------------------------------------------"
fi

exit

########################################
# Convert BED into RAW
########################################
if [[ ! -f $BASENAME_GEN".raw" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting genotypes into RAW format"
  echo "----------------------------------------------------------------------------------------------------"
  plink2 --bfile $GENOTYPES --keep $INDIVIDUALS --extract $VARIANTS --export A --out $BASENAME_GEN
  rm $BASENAME_GEN".log"
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes into RAW format because"
  echo $BASENAME_GEN".raw exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Cross-reference haplotypes and genotypes
########################################
if [[ ! -f $BASENAME_HAP".ref" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Cross-referencing haplotypes and genotypes"
  echo "----------------------------------------------------------------------------------------------------"
  ml gcc
  ml R
  Rscript --vanilla utils/crossref_alleles.R $CHR
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping Cross-reference of haplotypes and genotypes because"
  echo $BASENAME_HAP".ref exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
