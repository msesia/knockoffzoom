#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Generate the knockoffs.
#
# Authors: Matteo Sesia
# Date:    08/09/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -m method -r factor -k K [-n subjects] [--reset --resample]"
    echo "  -c chromosome         specify chromosome"
    echo "  -m partition method   LD measure and clustering method (Radj, Rsingle, Dadj, Dsingle)"
    echo "  -r resolution factor  resolution percentage (0,100)"
    echo "  -k motifs             specify number of haoplotype motifs in the fastPHASE model"
    echo "  -n subjects           number of subjects to construct knockoffs for (default: all)"
    echo "  --reset               delete all intermediate results and recompute"
    echo "  --resample            delete knockoffs and resample"
    exit 1
}

# Set default options
RESET=0
RESAMPLE=0
METHOD="Rsingle"
NSUBJECTS=-1

# Call getopt to validate theinput.
options=$(getopt -o c:m:r:k:n: --long reset,resample -- "$@")
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
    -k)
      shift;
      K=$1
      ;;
    -n)
      shift;
      NSUBJECTS=$1
      ;;
    --reset)
        RESET=1
        ;;
    --resample)
        RESAMPLE=1
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
# HMM knockoffs
##############################
# Create directory for knockoff files
KNO_DIR=$TMP_DIR"/knockoffs"
mkdir -p $KNO_DIR
OUT_DIR=$KNO_DIR"/"$METHOD$RESOLUTION"_K"$K
mkdir -p $OUT_DIR

# Location of this group's files
OUT_FILE=$OUT_DIR"/ukb_gen_chr"$CHR

# Generate knockoffs if not already present
if [[ ! -f $OUT_FILE".bed" ]] || [[ ! -f $OUT_FILE".bim" ]] || [[ ! -f $OUT_FILE".fam" ]] || [[ ! -f $OUT_FILE".key" ]] || [[ $RESET -eq 1 ]] || [[ $RESAMPLE -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Making knockoffs for chromosome "$CHR" with "$RESOLUTION"% resolution and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  # Generate knockoffs and save them in PED format
  Rscript --vanilla utils/knockoffs.R $CHR $RESOLUTION $METHOD $K $NSUBJECTS
  # Convert PED to BED
  plink --file $OUT_FILE --no-fid --no-parents --no-sex --no-pheno --keep-allele-order --make-bed --out $OUT_FILE
  # Add subject information to FAM file  
  SEX_FILE=$TMP_DIR"/fastphase/data/ukb_hap_chr"$CHR".sample"  
  awk 'FNR > 2 { print $1,$4 }' $SEX_FILE > $OUT_FILE".sex"
  awk '{ print $1,$2,$3,$4,$6 }' $OUT_FILE".fam"> $OUT_FILE".fam.tmp"
  join -1 1 -2 1 $OUT_FILE".fam.tmp" $OUT_FILE".sex" > $OUT_FILE".fam.tmp2"
  awk '{ print $1,$2,$3,$4,$6,$5 }' $OUT_FILE".fam.tmp2"> $OUT_FILE".fam"
  rm $OUT_FILE".fam.tmp"
  rm $OUT_FILE".fam.tmp2"
  rm $OUT_FILE".sex"
  # Remove PED file
  rm $OUT_FILE".ped"
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping knockoffs for chromosome "$CHR" with "$RESOLUTION"% resolution and method "$METHOD
  echo $OUT_FILE".bed exists"
  echo $OUT_FILE".bim exists"
  echo $OUT_FILE".fam exists"
  echo $OUT_FILE".key exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
