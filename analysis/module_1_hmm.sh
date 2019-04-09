#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Convert phased haplotypes and estimate HMM with fastPHASE
#
# Authors: Matteo Sesia
# Date:    07/19/2018

# Utility scripts
BGEN_TO_HAPST="../utils/bgen_to_hapst.sh"
HAPST_TO_INP="../utils/hapst_to_inp.sh"
VERIFY_HAPS="Rscript --vanilla ../utils/verify_haps.R"

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR

# List of chromosomes
CHR_LIST=$(seq 21 22)

# Which operations should we perform?
FLAG_CONVERT_HAP=0
FLAG_CONVERT_INP=0
FLAG_CHECK_INP=1
FLAG_RUN_FASTPHASE=0

##########################################
# Convert BGEN v1.2 into transposed HAPS #
##########################################

if [[ $FLAG_CONVERT_HAP == 1 ]]; then
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
    QC_SAMPLES="../data/qc/samples_qc.txt"

    # List of variants that passed QC
    QC_VARIANTS="../data/qc/variants_qc.txt"

    # Basename for output haplotype files (transposed HAPS format)
    HAP_BASENAME=$TMP_DIR"/example_chr"$CHR

    # Convert BGEN haplotypes to transposed HAPS file
    $BGEN_TO_HAPST -i $BGEN_BASENAME -c $CHR -n $QC_SAMPLES -v $QC_VARIANTS -b 500 -o $HAP_BASENAME

  done
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes into HAPS.T format"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Convert HAPS.T into INP
########################################
if [[ $FLAG_CONVERT_INP == 1 ]]; then
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

############################################
# Cross-reference haplotypes and genotypes #
############################################
if [[ $FLAG_CHECK_INP == 1 ]]; then

  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Cross-referencing haplotypes and genotypes"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Basename for output haplotype files (transposed HAPS format)
    HAP_BASENAME=$TMP_DIR"/example_chr"$CHR

    # Basename for the input genotype files (PLINK format)
    GENO_BASENAME="../data/genotypes/example_chr"$CHR

    # List of individuals that passed QC
    QC_SAMPLES="../data/qc/samples_qc.txt"

    # List of variants that passed QC
    QC_VARIANTS="../data/qc/variants_qc.txt"

    # Check whether the reference alleles match
    $VERIFY_HAPS $HAP_BASENAME $GENO_BASENAME $HAP_BASENAME $QC_SAMPLES $QC_VARIANTS

  done
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping cross-reference of haplotypes and genotypes"
  echo "----------------------------------------------------------------------------------------------------"
fi


###############################
# Estimate HMM with fastPHASE #
###############################
if [[ $FLAG_RUN_FASTPHASE == 1 ]]; then

  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Estimating HMM with fastPHASE"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # fastPHASE parameters
    FP_K=20         # Number of haplotype motifs
    FP_IT=10        # Number of EM iterations
    FP_SEED=123     # Random seed (this doesn't actually do anything because of bug in fastPhase)

    # Haplotypes in INP format
    HAP_INP=$TMP_DIR"/example_chr"$CHR".inp"

    # Basename for output files
    OUT_BASENAME=$TMP_DIR"/example_chr"$CHR

    # Run fastPHASE on this chromosome
    CMD="fastphase -Pp -T1 -K"$FP_K" -g -H-4 -B -C"$FP_IT" -S"$FP_SEED" -o"$OUT_BASENAME" "$HAP_INP
    $CMD

  done
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping estimation of HMM with fastPHASE"
  echo "----------------------------------------------------------------------------------------------------"
fi
