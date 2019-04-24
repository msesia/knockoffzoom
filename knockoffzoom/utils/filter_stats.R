#!/usr/bin/env Rscript

# Load packages
suppressMessages(library(tidyverse))

# Default arguments (for debugging)
stats.basename  <- "../../tmp/example_res2"
groups.basename <- "../../tmp/example_chr?_groups2.txt"
out.basename    <- "../../results/example_res2"

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
    summarise(BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min) %>%
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

# Save list of discoveries
out.file <- sprintf("%s_discoveries.txt", out.basename)
Selections %>% write_delim(out.file, delim=" ")
cat(sprintf("Discoveries written on: %s\n", out.file))

# Save list of test statistics
out.file <- sprintf("%s_stats.txt", out.basename)
Stats %>% write_delim(out.file, delim=" ")
cat(sprintf("Test statistics written on: %s\n", out.file))
