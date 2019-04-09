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
ld.file    <- "../tmp/example_chr21.ld"
bim.file   <- "../data/genotypes/example_chr21.bim"
qc.file    <- "../data/qc/variants_qc.txt"
out.file   <- "../tmp/example_chr21.RData"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
ld.file  <- as.character(args[1])
bim.file <- as.character(args[2])
qc.file  <- as.character(args[3])
out.file <- as.character(args[4])

# Load libraries
suppressMessages(library(tidyverse))
local.source("utils.R")

# Standard parameters
ld.measure <- "R"  # R^2

# Load list of variants
Variants <- read_delim(bim.file, delim="\t", col_names=c("CHR", "SNP", "X", "BP", "A1", "A2"), col_types=cols())

# Remove variants that did not pass QC
variants.qc <- read_delim(qc.file, delim=" ", col_names="SNP", col_types=cols())
Variants <- Variants %>% inner_join(variants.qc, by="SNP")

# Load sparse covariance matrix
cat("Loading covariance matrix... ")   
LD <- read_table2(ld.file, col_types=cols())
LD <- filter(LD, BP_A %in% Variants$BP, BP_B %in% Variants$BP)
cat("done.\n")

# Convert LD to a symmetric square matrix
cat("Converting LD table to LD matrix... ")   
LD.variants <- unique(c(LD$SNP_A, LD$SNP_B))
LD.positions <- unique(c(LD$BP_A, LD$BP_B))
Sigma <- ld.to.mat(LD, ld.measure, Variants$BP)
cat("done.\n")

# Convert covariance matrix to distance matrix
cat("Converting LD matrix to distance matrix... ")   
#Sigma.dist <- as.dist(1-abs(Sigma))                 # Faster, but requires lots of memory
Sigma.dist <- mat.to.dist(tril(Sigma,-1))            # Slower, but requires less memory
cat("done.\n")

cat("Computing clustering dendrogram... ")   
if(nrow(Sigma) > 10000) {
    Sigma.clust <- adjClust(Sigma.dist,h=10000)       # Makes contiguous groups
} else {
    Sigma.clust <- adjClust(Sigma.dist)               # Makes contiguous groups
}
cat("done.\n")

# Save dendrogram
save(list=c("Sigma.clust"), file=out.file)
cat(sprintf("Dendrogram saved in: %s\n", out.file))
