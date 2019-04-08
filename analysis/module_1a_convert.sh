#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Prepare phased haplotype for use with fastPHASE
#
# Authors: Matteo Sesia
# Date:    07/19/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -n holdout -b batch_size [-r]"
    echo "  -c chromosome   specify chromosome"
    echo "  -n holdout      size of test set"
    echo "  -b batch size   specify size of batches for conversion to haps"
    echo "  -r              delete all intermediate results and recompute"
    exit 1
}

# Default values for optional input arguments
RESET=0

# If less than 3 arguments supplied, display usage
if [  $# -le 4 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":c:b:n:r" opt; do
  case $opt in
    c)
      echo "  - chromosome            : $OPTARG" >&2
      CHR=$OPTARG
      ;;
    n)
      echo "  - size of test set      : $OPTARG" >&2
      N_TEST=$OPTARG
      ;;
    b)
      echo "  - batch size            : $OPTARG" >&2
      BATCH_SIZE=$OPTARG
      ;;
    r)
      echo "  - reset                 : TRUE" >&2
      RESET=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

############################################################
# Initialize variables
############################################################

# Load plink and datamash
export PATH=$PI_HOME/bin/:$PATH

# Location of input
DAT_DIR="/scratch/PI/candes/ukbiobank"               # Location of original data
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"           # Location of intermediate files
GENOTYPES=$DAT_DIR"/genotypes/ukb_gen_chr"$CHR
HAPLOTYPES=$DAT_DIR"/haplotypes/ukb_hap_chr"$CHR

# List of individuals and variants
INDIVIDUALS=$TMP_DIR"/QC_output/individuals_QC.txt"
VARIANTS=$TMP_DIR"/QC_output/QC_chr"$CHR".snplist"

# Location of output
FP_DIR=$TMP_DIR"/fastphase"
mkdir -p $FP_DIR
mkdir -p $FP_DIR"/data"
BASENAME_HAP=$FP_DIR"/data/ukb_hap_chr"$CHR
BASENAME_GEN=$FP_DIR"/data/ukb_gen_chr"$CHR
mkdir -p $FP_DIR"/train"
mkdir -p $FP_DIR"/test"
BASENAME_TRAIN=$FP_DIR"/train/ukb_hap_chr"$CHR
BASENAME_TEST=$FP_DIR"/test/ukb_hap_chr"$CHR

########################################
# Convert BGEN v1.2 into HAPS.T
########################################

#CONVERT_BGEN_HAPS=0

CONVERT_BGEN_HAPS=1
# Check whether the HAPS.T file has the correct number of rows
if [[ -s $BASENAME_HAP".haps.t" ]] && [[ $RESET -eq 0 ]]; then
  echo "Counting number of lines in "$BASENAME_HAP".haps.t ..."
  N_LINES=($(wc -l $BASENAME_HAP".haps.t"))
  N_LINES=${N_LINES[0]}
  echo "Number of lines = "$N_LINES
  if [[ $N_LINES == 700238 ]]; then
    CONVERT_BGEN_HAPS=0
  fi
fi

if [[ $CONVERT_BGEN_HAPS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes into HAPS.T format"
  echo "----------------------------------------------------------------------------------------------------"
  utils/bgen_to_hapst.sh -i $HAPLOTYPES -c $CHR -n $INDIVIDUALS -v $VARIANTS -b $BATCH_SIZE -o $BASENAME_HAP
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes into HAPS.T format because"
  echo $BASENAME_HAP".haps.t"" exists and has $N_LINES lines"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Convert HAPS.T into INP
########################################
if [[ ! -f $BASENAME_HAP".inp" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes into INP format"
  echo "----------------------------------------------------------------------------------------------------"
  utils/hapst_to_inp.sh -i $BASENAME_HAP -o $BASENAME_HAP
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes into INP format because"
  echo $BASENAME_HAP".inp"" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

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
