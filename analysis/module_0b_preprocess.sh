#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Pre-process data downloaded from UK Biobank server
#
# Authors: Matteo Sesia
# Date:    06/24/2018

###################
# Main parameters #
###################

# Data storage location
DATA_DIR="/scratch/PI/candes/ukbiobank"

# Application number
APP=xxx

##############################
# Parse input
##############################
PROGRAM_NAME=$0
function display_usage {
    echo "Usage: $PROGRAM_NAME [--help --resample]"
    echo "  --help          print this usage message and quit"
    echo "  --reset         force re-download of all data"
    exit 1
}
# Set default options
RESET=0
HELP=0

# Call getopt to validate theinput.
options=$(getopt -o h --long help,reset -- "$@")
[ $? -eq 0 ] || {
    display_usage
}
eval set -- "$options"
while true; do
    case "$1" in
    --reset)
        RESET=1
        ;;
    -h | --help)
        HELP=1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done

if [[ $HELP -eq 1 ]]; then
  display_usage
fi

########################################
# Haplotypes
########################################
# Set storage directory for genotype data
OUT_DIR=$DATA_DIR"haplotypes/"
mkdir -p $OUT_DIR
for CHR in {1..21}
do
  # Compute BIM files for haplotypes
  BASENAME="ukb_hap_chr"$CHR
  OUT_FILE=$BASENAME".bim"
  if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Creating haplotype map (BIM) for chromosome "$CHR
    echo "----------------------------------------------------------------------------------------------------"
    plink2 --bgen $OUT_DIR$BASENAME".bgen" --sample $OUT_DIR$BASENAME".sample" --oxford-single-chr $CHR \
           --make-just-bim --out $OUT_DIR$BASENAME
    rm $OUT_DIR$BASENAME".log"
  else
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Skipping haplotype map (BIM) for chromosome "$CHR" because"
    echo $OUT_DIR$OUT_FILE" exists"
    echo "----------------------------------------------------------------------------------------------------"
  fi
  
done

