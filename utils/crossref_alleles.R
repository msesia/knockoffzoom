#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Input arguments
chr = args[1]
n.load = 10000

# Load libraries
suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
source("utils.R")

# Location of data
fp.dir <- "/scratch/PI/candes/ukbiobank_tmp/fastphase/data"

# Load haplotypes from INP file
cat(sprintf("Loading haplotypes... "))
basename <- paste(fp.dir, "/ukb_hap_chr", chr, sep="")
Haplotypes <- read.inp(basename, progress=FALSE, n.load=n.load)
subjects <- rownames(Haplotypes)
Haplotypes <- Haplotypes %>% as_tibble() %>% mutate(Subject=subjects) %>% dplyr::select(Subject, everything())
cat(sprintf("done.\n"))

# Load genotypes
cat(sprintf("Loading genotypes... "))
basename <- paste(fp.dir, "/ukb_gen_chr", chr, sep="")
Genotypes <- read.raw(basename, progress=FALSE, n.load=n.load)
Genotypes <- dplyr::select(Genotypes, colnames(Haplotypes))
Genotypes <- filter(Genotypes, Subject %in% Haplotypes$Subject)
cat(sprintf("done.\n"))

# Convert to numeric matrices
H <- as.matrix(dplyr::select(Haplotypes, -Subject))
X <- as.matrix(dplyr::select(Genotypes, -Subject))
X.H <- H[seq(1,nrow(H),by=2),] + H[seq(2,nrow(H),by=2),]

# Match alleles
cat(sprintf("Matching alleles... "))
flipped <- sapply(1:ncol(X), function(j) {
    not.missing <- which(!is.na(X[,j]))
    error.original <- sum(X[not.missing,j] != X.H[not.missing,j])
    error.flipped  <- sum(X[not.missing,j] != 2-X.H[not.missing,j])
    if(error.original == 0) return(0)
    else if(error.flipped == 0) return(1)
    else return(-1)
})
cat(sprintf("done.\n"))
# Verify that there are no errors, then flip the alleles
stopifnot(sum(flipped==-1)==0)
cat(sprintf("%d alleles out of %d will be flipped. ", sum(flipped==1), length(flipped)))
cat("No errors were found in the converted data.\n")
Map <- tibble(Variant=colnames(Genotypes)[-1], Flip=as.integer(flipped))

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
basename <- paste(fp.dir, "/ukb_hap_chr", chr, sep="")
write_delim(Map, paste(basename, ".ref", sep=""), delim="\t")
