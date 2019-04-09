#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Convert phased haplotype files from BGEN 1.2 to transposed HAPS
#
# Authors: Matteo Sesia
# Date:    07/19/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -i input -c chromosome -n individuals -v variants -b batch -o output"
    echo "  -i input           basename for BGEN/sample files"
    echo "  -c chromosome      chromosome number"
    echo "  -n individiduals   list of individuals"
    echo "  -v variants        list of variants"
    echo "  -b batch           batch size"
    echo "  -o output          basename for output HAPS.T file"
    exit 1
}

# Default values for optional input arguments
RESET=0

# If less than 6 arguments supplied, display usage
if [  $# -le 7 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":i:c:n:v:b:o:" opt; do
  case $opt in
    i)
      echo "  - input                 : $OPTARG" >&2
      INPUT=$OPTARG
      ;;
    c)
      echo "  - chromosome            : $OPTARG" >&2
      CHR=$OPTARG
      ;;
    n)
      echo "  - individuals           : $OPTARG" >&2
      INDIVIDUALS=$OPTARG
      ;;
    v)
      echo "  - variants              : $OPTARG" >&2
      VARIANTS=$OPTARG
      ;;
    b)
      echo "  - batch                 : $OPTARG" >&2
      BATCH_SIZE=$OPTARG
      ;;
    o)
      echo "  - output                : $OPTARG" >&2
      OUTPUT=$OPTARG
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
echo ""

# Divide individuals into batches
INDIVIDUALS_BATCHES=$OUTPUT".individuals.batch"
echo "----------------------------------------------------------------------------------------------------"
split -d -l $BATCH_SIZE --verbose $INDIVIDUALS $INDIVIDUALS_BATCHES
echo "----------------------------------------------------------------------------------------------------"

# Remove old HAPS.T file
HAPS_T=$OUTPUT".haps.t"
if [ -f $HAPS_T ]; then
  rm $HAPS_T
fi

# Remove old sample file
SAM_FILE=$OUTPUT".sample"
if [ -f $SAM_FILE ]; then
  rm $SAM_FILE
fi

# Create transposed HAPS files in batches
for INDIVIDUALS_BATCH in $INDIVIDUALS_BATCHES*; do
  BATCH_IDX="${INDIVIDUALS_BATCH##*.}"

  # Convert haplotypes into HAPS format for all individuals in this batch
  echo ""
  echo "Converting BGEN haplotypes into transposed HAPS format for $BATCH_IDX ..."
  HAPS_BATCH=$OUTPUT".batch"
  plink2 --bgen $INPUT".bgen" --sample $INPUT".sample" --oxford-single-chr $CHR \
         --keep $INDIVIDUALS_BATCH --extract $VARIANTS --export hapslegend --out $OUTPUT".batch"

  # Transpose HAPS for this batch and append them to the full HAPS.T file
  echo "Transposing haplotypes for "$BATCH_IDX" ..."
  datamash transpose -W --output-delimiter=" " < $OUTPUT".batch.haps" > $HAPS_T"."$BATCH_IDX

  # Append sample file for this batch to the full sample file
  SAM_FILE=$OUTPUT".sample"
  if [ ! -f $SAM_FILE ]; then
    head -n 2 $OUTPUT".batch.sample" > $SAM_FILE
  fi
  tail -n +3 $OUTPUT".batch.sample" >> $SAM_FILE

  # Append HAPS.T file for this batch to the full HAPS.T file
  echo "Appending transposed HAPS file for "$BATCH_IDX" to the full file for chromosome "$CHR" ..."
  cat $HAPS_T"."$BATCH_IDX >> $HAPS_T

  # Remove temporary files
  rm $HAPS_T"."$BATCH_IDX
  rm $INDIVIDUALS_BATCH
  rm $OUTPUT".batch.haps"
  rm $OUTPUT".batch.sample"
  rm $OUTPUT".batch.log"

done

echo ""
echo "Done. Results written on: "$HAPS_T

# Rename temporary legend file
mv $OUTPUT".batch.legend" $OUTPUT".legend"
