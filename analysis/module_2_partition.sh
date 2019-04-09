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
CLUSTER_VARIANTS="Rscript --vanilla ../utils/cluster.R"
PARTITION_VARIANTS="Rscript --vanilla ../utils/partition.R"

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR

# List of chromosomes
CHR_LIST=$(seq 21 22)

# List of resolutions
RESOLUTION_LIST=("2" "10")

# Which operations should we perform?
FLAG_COMPUTE_LD=0
FLAG_CLUSTER=0
FLAG_PARTITION=1

#####################
# Compute LD matrix #
#####################

if [[ $FLAG_COMPUTE_LD == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing variant statistics"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Basename for the input genotype files (PLINK format)
    GENO_BASENAME="../data/genotypes/example_chr"$CHR

    # List of individuals that passed QC
    QC_SAMPLES="../data/qc/samples_qc.txt"

    # List of variants that passed QC
    QC_VARIANTS="../data/qc/variants_qc.txt"

    # List of variants that passed QC
    BGEN_VARIANT="../data/qc/variants_qc.txt"

    # Basename for output haplotype files (transposed HAPS format)
    OUT_BASENAME=$TMP_DIR"/example_chr"$CHR

    # Compute variant frequencies and LD table for this chromosome
    plink \
      --bfile $GENO_BASENAME \
      --keep $QC_SAMPLES --extract $QC_VARIANTS \
      --r2 dprime --ld-window 1000 --ld-window-kb 1000 --ld-window-r2 0.01 \
      --freq \
      --out $OUT_BASENAME

  done
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping variant statistics"
  echo "----------------------------------------------------------------------------------------------------"
fi

##############
# Clustering #
##############

if [[ $FLAG_CLUSTER == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Clustering variants"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Basename for the input genotype files (PLINK format)
    GENO_BIM="../data/genotypes/example_chr"$CHR".bim"

    # List of variants that passed QC
    QC_VARIANTS="../data/qc/variants_qc.txt"

    # Basename for output haplotype files (transposed HAPS format)
    LD_TABLE=$TMP_DIR"/example_chr"$CHR".ld"

    # Basename for output dendrogram file
    OUT_FILE=$TMP_DIR"/example_chr"$CHR".RData"

    # Compute adjacent clustering dendrogram
    $CLUSTER_VARIANTS $LD_TABLE $GENO_BIM $QC_VARIANTS $OUT_FILE
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping variant clustering"
  echo "----------------------------------------------------------------------------------------------------"
fi

################
# Partitioning #
################
if [[ $FLAG_PARTITION == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Partitioning variants"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do
    for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing chromosome "$CHR" at resolution "$RESOLUTION" ..."
    echo ""

    # Basename for the input genotype files (PLINK format)
    GENO_BIM="../data/genotypes/example_chr"$CHR".bim"

    # List of variants that passed QC
    QC_VARIANTS="../data/qc/variants_qc.txt"

    # Basename for output haplotype files (transposed HAPS format)
    LD_TABLE=$TMP_DIR"/example_chr"$CHR".ld"

    # Basename for dendrogram file
    DENDRO_FILE=$TMP_DIR"/example_chr"$CHR".RData"

    # Basename for output partition file
    GROUPS_FILE=$TMP_DIR"/example_chr"$CHR

    # Compute adjacent clustering dendrogram
    $PARTITION_VARIANTS $DENDRO_FILE $GENO_BIM $QC_VARIANTS $RESOLUTION $GROUPS_FILE

    done
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping variant clustering"
  echo "----------------------------------------------------------------------------------------------------"
fi
