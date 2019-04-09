#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Run fastPHASE on each chromosome separately.
#
# Author: Matteo Sesia
# Date:   06/20/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -k K [-n samples] [-u] "
    echo "  -c chromosome   specify chromosome"
    echo "  -k K            number of motifs in fastPHASE model"
    echo "  -u              whether to interpret the genotypes as unphased"
    exit 1
}

# Default values for optional input arguments
PHASE="phased"
PHASE_FLAG="-B"

# If less than 2 arguments supplied, display usage
if [  $# -le 3 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":c:k:ur" opt; do
  case $opt in
    c)
      echo "  - chromosome            : $OPTARG" >&2
      CHR=$OPTARG
      ;;
    k)
      echo "  - fastPhase motifs      : $OPTARG" >&2
      K=$OPTARG
      ;;
    u)
      echo "  - unphased              : TRUE" >&2
      PHASE="unphased"
      PHASE_FLAG=""
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
# Run fastPHASE
############################################################

# Load plink
export PATH=$PI_HOME/bin/:$PATH

# Location of input
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"
INP_DIR=$TMP_DIR"/fastphase/train"
INP_FILE=$INP_DIR/"ukb_hap_chr"$CHR".inp"

# Location of output
OUT_DIR=$TMP_DIR"/fastphase/"$PHASE"_K"$K
mkdir -p $OUT_DIR
chmod ugo-rwx $OUT_DIR
chmod u+rwx $OUT_DIR

OUT_FILE=$OUT_DIR"/ukb_hap_chr"$CHR

# Other fastPHASE parameters
NIT=10           # Number of iterations (DO NOT USE MORE THAN 20 ish ITERATIONS or bad things will happen)
SED=123          # Random seed (this doesn't actually do anything because of bug in fastPhase)

# Run fastPHASE on this chromosome
CMD="fastphase -Pp -T1 -K"$K" -g -H-4 "$PHASE_FLAG" -C"$NIT" -S"$SED" -o"$OUT_FILE" "$INP_FILE
$CMD

# Evaluate imputation error
./test_impute.sh -c $CHR -k $K
