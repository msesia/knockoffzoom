#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Install packages bigsnpr and bigstatsr
# devtools::install_github("privefl/bigstatsr")
# devtools::install_github("privefl/bigsnpr")

# Documentation here:
# https://privefl.github.io/bigsnpr/reference/index.html
# https://privefl.github.io/bigstatsr/reference/big_spLinReg.html

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(devtools))
suppressMessages(library(bigsnpr))

# Default arguments (for debugging)
stats.basename  <- "../tmp/example_res2"
groups.basename <- "../tmp/example_chr?_groups2.txt"
out.basename    <- "../results/example_res2"
                                        
# Input arguments
args <- commandArgs(trailingOnly=TRUE)
stats.basename  <- as.character(args[1])
groups.basename <- as.character(args[2])
out.basename    <- as.character(args[3])

##############################
## Load the test statistics ##
##############################

# Loading test statistics
stats.file <- sprintf("%s_stats.txt", stats.basename)
Stats <- read_delim(stats.file, delim=" ", col_types=cols())

# Load variant partitions
chr.list <- unique(Stats$CHR)
Variants <- lapply(chr.list, function(chr) {
    group.file <- gsub(fixed("[?]"), chr, groups.basename)
    Variants.chr <- read_delim(group.file, delim=" ", col_types=cols()) %>%
        select(CHR, SNP, BP, Group)
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants) %>% arrange(CHR,BP)

# Cross reference stats and list of variants
Stats <- Stats %>% left_join(Variants, by = c("CHR", "Group")) %>%
    group_by(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
    summarise(BP.min=min(BP), BP.max=max(BP)) %>%
    ungroup() %>%
    arrange(desc(abs(W)))

###############################
## Apply the knockoff filter ##
###############################

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
knockoff.filter <- function(Stats, fdr=0.1, offset=1) {
    W.thres <- knockoff.threshold(Stats$W, fdr=fdr, offset=offset)
    Selected <- Stats %>% filter(W >= W.thres)
    return(Selected)
}
Selections <- Stats %>% knockoff.filter(fdr=0.1, offset=1)

# Give preview
cat(sprintf("Selections:\n"))
Selections %>% print()

# Save selections
out.file <- sprintf("%s.txt", out.basename)
Selections %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written on: %s\n", out.file))
