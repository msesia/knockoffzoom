#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Partition the variants.
#
# Authors: Matteo Sesia
# Date:    08/09/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -m method -r factor [--reset]"
    echo "  -c chromosome         specify chromosome"
    echo "  -m partition method   LD measure and clustering method (Radj, Rsingle, Dadj, Dsingle)"
    echo "  -r resolution factor  resolution percentage (0,100)"
    echo "  --reset               delete all intermediate results and recompute"
    exit 1
}

# Set default options
RESET=0
METHOD="Rsingle"

# Call getopt to validate theinput.
options=$(getopt -o c:m:r: --long reset -- "$@")
[ $? -eq 0 ] || {
    echo "Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    -c)
      shift;
      CHR=$1
      [[ ! $CHR =~ 1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22 ]] && {
        echo "Incorrect value for option -c provided"
        exit 1
      }
      ;;
    -m)
      shift;
      METHOD=$1
      [[ ! $METHOD =~ "Rsingle"|"Radj"|"Dsingle"|"Dadj" ]] && {
        echo "Incorrect value for option -m provided"
        exit 1
      }
      ;;
    -r)
      shift;
      RESOLUTION=$1
      ;;
    --reset)
        RESET=1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done

# Verify that all required parameters have been supplied
if [[ -z "$CHR" ]] || [[ -z "$METHOD" ]] || [[ -z "$RESOLUTION" ]] || [[ -z "$K" ]]; then
  display_usage
fi

# Initialize environment
DAT_DIR="/scratch/PI/candes/ukbiobank"             # Location of original data
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"         # Location of intermediate files

##############################
# Resolution
##############################

ml gcc
ml R/3.5.1
export PATH=$PI_HOME/bin/:$PATH                    # Load PLINK

# Location of data for this chromosome
GENOTYPES=$DAT_DIR"/genotypes/ukb_gen_chr"$CHR
HAPLOTYPES=$DAT_DIR"/haplotypes/ukb_hap_chr"$CHR

# List of individuals and variants that were fed in to fastPHASE
INDIVIDUALS=$TMP_DIR"/fastphase/data/individuals.txt"
VARIANTS=$TMP_DIR"/fastphase/data/ukb_hap_chr"$CHR".legend"

# Create directory for resolution files
GRP_DIR=$TMP_DIR"/resolution"
mkdir -p $GRP_DIR

# Create directory for dendrograms
mkdir -p $GRP_DIR"/"$METHOD

# Create directory for group files
OUT_DIR=$GRP_DIR"/"$METHOD$RESOLUTION
mkdir -p $OUT_DIR

# Location of group files
OUT_FILE_1=$OUT_DIR"/grp_chr"$CHR".txt"
OUT_FILE_2=$OUT_DIR"/rep_chr"$CHR".txt"

# Compute variant frequencies if not already present
if [[ ! -f $OUT_FILE_1 ]] || [[ ! -f $OUT_FILE_2 ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Partitioning variants for chromosome "$CHR" with "$RESOLUTION"% resolution and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  Rscript --vanilla utils/partition.R $CHR $RESOLUTION $METHOD
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping partitioning for chromosome "$CHR" at "$RESOLUTION"% resolution and method "$METHOD" because"
  echo $OUT_FILE_1" exists"
  echo $OUT_FILE_2" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
