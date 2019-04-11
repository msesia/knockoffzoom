#!/bin/bash
#
#

# Range of chromosomes to include in the analysis
CHR_MIN=21
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

GENO_FILE="../data/genotypes/example_chr"
FAM_FILE="../data/genotypes/example_chr"$CHR_MIN".fam"

PHENO_FILE="../data/phenotypes/phenotypes.tab"
PHENO_NAME="y"

STATS_FILE="../results/example_lmm.txt"

# Stuff for bolt
LD_TABLE="/home/groups/candes/Software/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz"
MAP_TABLE="/home/groups/candes/Software/BOLT-LMM_v2.3.2/tables/genetic_map_hg17_withX.txt.gz"

bolt \
    --bed=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bed" \
    --bim=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bim" \
    --fam=$FAM_FILE \
    --maxMissingPerSnp=1 \
    --phenoFile=$PHENO_FILE \
    --phenoCol=$PHENO_NAME \
    --covarFile=$PHENO_FILE \
    --covarCol="sex" \
    --LDscoresFile=$LD_TABLE \
    --geneticMapFile=$MAP_TABLE \
    --statsFile=$STATS_FILE \
    --lmm \
    --numThreads=1

echo "Output file:"
echo $STATS_FILE
