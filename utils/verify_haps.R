#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Input arguments
args = commandArgs(trailingOnly=TRUE)

hap.basename  <- as.character(args[1])
geno.basename <- as.character(args[2])
ref.basename  <- as.character(args[3])

cat(sprintf("Checking references alleles...\n"))

# Load libraries
suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
source("../utils/utils.R")

#############################
## Check reference alleles ##
#############################

# Load haplotypes from INP file
cat(sprintf("Loading haplotypes... "))
Haplotypes <- read.inp(hap.basename, progress=FALSE)
subjects <- rownames(Haplotypes)
Haplotypes <- Haplotypes %>% as_tibble() %>% mutate(Subject=subjects) %>% dplyr::select(Subject, everything())
cat(sprintf("done.\n"))

# Load genotypes
bed.file <- sprintf("%s.bed", geno.basename)
bim.file <- sprintf("%s.bim", geno.basename)
fam.file <- sprintf("%s.fam", geno.basename)
geno.dat <- read.plink(bed.file, bim.file, fam.file)

# Convert to numeric matrices
H <- as.matrix(dplyr::select(Haplotypes, -Subject))
X <- as(geno.dat$genotypes, "numeric")
H.X <- H[seq(1,nrow(H),by=2),] + H[seq(2,nrow(H),by=2),]

# Match alleles
cat(sprintf("Matching alleles... "))
flipped <- sapply(1:ncol(X), function(j) {
    not.missing <- which(!is.na(X[,j]))
    error.original <- sum(X[not.missing,j] != H.X[not.missing,j])
    error.flipped  <- sum(X[not.missing,j] != 2-H.X[not.missing,j])
    if(error.original == 0) return(0)
    else if(error.flipped == 0) return(1)
    else return(-1)
})
cat(sprintf("done.\n"))
# Verify that there are no errors, then flip the alleles
stopifnot(sum(flipped==-1)==0)
cat(sprintf("%d alleles out of %d will be flipped. ", sum(flipped==1), length(flipped)))
cat("No errors were found in the converted data.\n")
Map <- tibble(SNP=colnames(geno.dat$genotypes), Flip=as.integer(flipped))

# Check whether any alleles are ambiguous
is.ambiguous <- sapply(1:ncol(X), function(j) {
    not.missing <- which(!is.na(X[,j]))
    ambiguous <- sum(X[not.missing,j]==1) == length(not.missing)
    return(ambiguous)
})
n.ambiguous <- sum(is.ambiguous)
cat(sprintf("The reference allele for %d out of %d SNPs is ambiguous.\n", n.ambiguous, ncol(X)))
if(n.ambiguous > 0) {
    cat("WARNING! You should cross-reference haplotypes and genotypes again, using more data. \n")
}

# Write map of flipped alleles to file
ref.file <- sprintf("%s.ref", ref.basename)
write_delim(Map, ref.file, delim="\t")
cat(sprintf("List of alleles to be flipped written on: %s\n", ref.file))
