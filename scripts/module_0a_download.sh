#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Download data from UK Biobank server
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
    echo "Usage: $PROGRAM_NAME [-c chromosome --help --resample]"
    echo "  -c chromosome   specify chromosome"
    echo "  -r              force re-download of all data"
    exit 1
}
# Set default options
RESET=0
CHR=0

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":c:hr" opt; do
  case $opt in
    c)
      echo "  - chromosome            : $OPTARG" >&2
      CHR=$OPTARG
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

# Create target storage location
mkdir -p $DATA_DIR

# Set working directory
WORK_DIR=$DATA_DIR"/tmp"
mkdir -p $WORK_DIR
CURR_DIR="$PWD"

# Copy useful scripts
cp .ukbkey $WORK_DIR"/"
cp ukbfetch $WORK_DIR"/"
cp ukbgene $WORK_DIR"/"

# Change working directory
echo "Changing working directory from:"
echo " "$CURR_DIR
echo " to:"
cd $WORK_DIR
echo $PWD

if [ $CHR == 0 ]; then
  CHROMOSOMES=({1..22} "X" "Y" "XY" "MT")
  CHROMOSOMES_HAP=({1..22})
else 
  CHROMOSOMES=($CHR)
  CHROMOSOMES_HAP=($CHR)
fi

RELEASE=s488349

########################################
# Genotypes
########################################
# Set storage directory for genotype data
OUT_DIR=$DATA_DIR"/genotypes/"
mkdir -p $OUT_DIR
# Download BIM files for all chromosomes
function check_bim_exist {
  # Check whether any BIM files are missing
  for CHR in ${CHROMOSOMES[@]}
  do
    if [[ ! -f $OUT_DIR"ukb_gen_chr"$CHR".bim" ]] ; then
      BIM_EXIST=0
      return
    fi
  done  
  BIM_EXIST=1
}

check_bim_exist
if [[ ! $BIM_EXIST -eq 1 ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Downloading genotype data (BIM) for all chromosomes"
  echo "----------------------------------------------------------------------------------------------------"
  OUT_FILE="ukb_snp_bim.tar"
  wget -nd biobank.ndph.ox.ac.uk/showcase/auxdata/$OUT_FILE
  tar -xvf $OUT_FILE --directory $OUT_DIR
  for CHR in ${CHROMOSOMES[@]}
  do
    mv $OUT_DIR"ukb_snp_chr"$CHR"_v2.bim" $OUT_DIR"ukb_gen_chr"$CHR".bim"
  done
  rm $OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping download of genotype data (BIM) for all chromosomes because"
  for CHR in ${CHROMOSOMES[@]}
  do
    echo $OUT_DIR"ukb_gen_chr"$CHR".bim exists"
  done  
  echo "----------------------------------------------------------------------------------------------------"
fi

for CHR in ${CHROMOSOMES[@]}
do
  # Download BED file for this chromosome
  OUT_FILE="ukb_gen_chr"$CHR".bed"
  if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Downloading genotype data (BED) for chromosome "$CHR
    echo "----------------------------------------------------------------------------------------------------"
    ./ukbgene cal -c$CHR -v
    mv "ukb_cal_chr"$CHR"_v2.bed" $OUT_DIR$OUT_FILE
  else
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Skipping download of genotype data (BED) for chromosome "$CHR" because"
    echo $OUT_DIR$OUT_FILE" exists"
    echo "----------------------------------------------------------------------------------------------------"
  fi

  # Download FAM file for this chromosome
  OUT_FILE="ukb_gen_chr"$CHR".fam"
  if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Downloading genotype data (FAM) for chromosome "$CHR
    echo "----------------------------------------------------------------------------------------------------"
    ./ukbgene cal -c$CHR -m -v
    mv "ukb"$APP"_cal_chr"$CHR"_v2_"$RELEASE".fam" $OUT_DIR$OUT_FILE
  else
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Skipping download of genotype data (FAM) for chromosome "$CHR" because"
    echo $OUT_DIR$OUT_FILE" exists"
    echo "----------------------------------------------------------------------------------------------------"
  fi
done

########################################
# Haplotypes
########################################
# Set storage directory for genotype data
OUT_DIR=$DATA_DIR"/haplotypes/"
mkdir -p $OUT_DIR
for CHR in ${CHROMOSOMES_HAP[@]}
do
  # Download BGEN file for this chromosome
  OUT_FILE="ukb_hap_chr"$CHR".bgen"
  if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Downloading haplotype data (BGEN) for chromosome "$CHR
    echo "----------------------------------------------------------------------------------------------------"
    ./ukbgene hap -c$CHR -v
    mv "ukb_hap_chr"$CHR"_v2.bgen" $OUT_DIR$OUT_FILE
  else
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Skipping download of haplotype data (BGEN) for chromosome "$CHR" because"
    echo $OUT_DIR$OUT_FILE" exists"
    echo "----------------------------------------------------------------------------------------------------"
  fi

  # Download SAMPLE file for this chromosome
  OUT_FILE="ukb_hap_chr"$CHR".sample"
  if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Downloading haplotype data (SAMPLE) for chromosome "$CHR
    echo "----------------------------------------------------------------------------------------------------"
    ./ukbgene hap -c$CHR -m -v
    mv "ukb"$APP"_hap_chr"$CHR"_v2_s487381.sample" $OUT_DIR$OUT_FILE
  else
    echo ""
    echo "----------------------------------------------------------------------------------------------------"
    echo "Skipping download of haplotype data (SAMPLE) for chromosome "$CHR" because"
    echo $OUT_DIR$OUT_FILE" exists"
    echo "----------------------------------------------------------------------------------------------------"
  fi  
done

# Download BGI files
function check_bgi_exist {
  # Check whether any BIM files are missing
  for CHR in ${CHROMOSOMES_HAP[@]}
  do
    if [[ ! -f $OUT_DIR"ukb_hap_chr"$CHR".bgen.bgi" ]] ; then
      BGI_EXIST=0
      return
    fi
  done  
  BGI_EXIST=1
}
check_bgi_exist
if [[ ! $BGI_EXIST -eq 1 ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Downloading haplotype data (BGI) for all chromosomes"
  echo "----------------------------------------------------------------------------------------------------"
  wget -nd biobank.ndph.ox.ac.uk/showcase/auxdata/ukb_hap_bgi.tgz
  tar -xvf ukb_hap_bgi.tgz --directory $OUT_DIR
  for CHR in ${CHROMOSOMES_HAP[@]}
  do
    mv $OUT_DIR"ukb_hap_chr"$CHR"_v2.bgen.bgi" $OUT_DIR"ukb_hap_chr"$CHR".bgen.bgi"
  done
  rm ukb_hap_bgi.tgz
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping download of haplotype data (BGI) for all chromosomes because"
  for CHR in ${CHROMOSOMES_HAP[@]}
  do
    echo $OUT_DIR"ukb_hap_chr"$CHR".bgen.bgi"
  done  
  echo "----------------------------------------------------------------------------------------------------"
fi


########################################
# Quality control
########################################

# Download relatedness data
OUT_DIR=$DATA_DIR"/quality-control/"
OUT_FILE="ukb_rel.dat"
if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Downloading relatedness data"
  echo "----------------------------------------------------------------------------------------------------"
  ./ukbgene rel -v
  mkdir -p $OUT_DIR
  mv "ukb"$APP"_rel_"$RELEASE".dat" $OUT_DIR$OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping download of relatedness data because"
  echo $OUT_DIR$OUT_FILE" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

# Download SNP quality control data
OUT_DIR=$DATA_DIR"/quality-control/"
OUT_FILE="ukb_snp_qc.txt"
if [[ ! -f $OUT_DIR$OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Downloading SNP quality control data"
  echo "----------------------------------------------------------------------------------------------------"
  mkdir -p $OUT_DIR
  wget -nd biobank.ndph.ox.ac.uk/showcase/auxdata/ukb_snp_qc.txt -O $OUT_DIR$OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping download of SNP quality control data"
  echo $OUT_DIR$OUT_FILE" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Phenotypes
########################################
# These must be downloaded from the UK Biobank data management website

