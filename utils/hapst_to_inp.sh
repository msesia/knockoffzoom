#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Convert phased haplotype files transposed HAPS to INP
#
# Authors: Matteo Sesia
# Date:    07/19/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -i input -o output"
    echo "  -i input           basename for HAPS.T/sample/legend files"
    echo "  -o output          basename for output INP file"
    exit 1
}

# Default values for optional input arguments
RESET=0

# If less than 2 arguments supplied, display usage
if [  $# -le 3 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":i:o:" opt; do
  case $opt in
    i)
      echo "  - input                 : $OPTARG" >&2
      INPUT=$OPTARG
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

# Convert transposed HAPS file into INP file
HAPS_T=$INPUT".haps.t"
SAMPLE=$INPUT".sample"
LEGEND=$INPUT".legend"
INP_FILE=$OUTPUT".inp"

# Print header of inp file
TMP1_FILE=$INP_FILE".tmp1"
TMP2_FILE=$INP_FILE".tmp2"
awk 'END{print (NR-2)}' $SAMPLE > $TMP1_FILE     # Number of individuals
awk 'END{print (NR-1)}' $LEGEND >> $TMP1_FILE    # Number of variants

# Print list of variant positions
awk 'BEGIN{print "P"} NR>1{print $2}' $LEGEND | datamash transpose -W --output-delimiter=" " >> $TMP1_FILE

# Add haplotypes with spacing between individuals to inp file
awk 'BEGIN { OFS=""} {$1=$1; if(NR%2) {printf "# ID \n"$0"\n"} else {print $0}}' $HAPS_T >> $TMP1_FILE

# Add individual IDs to the inp file
awk 'BEGIN {printf "\n\n\n"} NR>2{printf $1"\n\n\n"}' $SAMPLE > $TMP2_FILE
paste -d "" $TMP1_FILE $TMP2_FILE > $INP_FILE

# Clean up remporary files
rm $TMP1_FILE $TMP2_FILE

echo ""
echo "Done. Results written on: "$INP_FILE
