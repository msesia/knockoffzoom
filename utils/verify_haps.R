#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

LocationOfThisScript = function() {
    # Solution from:
    # https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
    
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
local.source <- function(script.name) {
    current.dir <- LocationOfThisScript()
    source(sprintf("%s/%s", current.dir, script.name))
}

# Default arguments (for debugging)
hap.basename     <- "../tmp/example_chr21"
geno.basename    <- "../data/genotypes/example_chr21"
ref.basename     <- "../tmp/example_chr21"
qc.samples.file  <- "../data/qc/samples_qc.txt"
qc.variants.file <- "../data/qc/variants_qc.txt"

# Input arguments
args = commandArgs(trailingOnly=TRUE)
hap.basename     <- as.character(args[1])
geno.basename    <- as.character(args[2])
ref.basename     <- as.character(args[3])
qc.samples.file  <- as.character(args[4])
qc.variants.file <- as.character(args[5])

cat(sprintf("Checking references alleles...\n"))

# Load libraries
suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
local.source("utils.R")

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

# Load results of qc
samples.qc <- read_delim(qc.samples.file, delim=" ", col_types=cols(.default=col_character()),
                         col_names=c("FID","IID"))
variants.qc <- read_delim(qc.variants.file, delim=" ", col_types=cols(.default=col_character()), col_names="SNP")

# Extract genotypes that passed QC
keep.rows <- which(geno.dat$fam$pedigree %in% parse_number(samples.qc$FID))
keep.cols <- which(geno.dat$map$snp.name %in% variants.qc$SNP)

# Convert to numeric matrices
H <- as.matrix(dplyr::select(Haplotypes, -Subject))
X <- as(geno.dat$genotypes[keep.rows, keep.cols], "numeric")
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
Map <- tibble(SNP=colnames(geno.dat$genotypes[keep.cols]), Flip=as.integer(flipped))

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
