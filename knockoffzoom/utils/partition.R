#!/usr/bin/env Rscript

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
dendro.file  <- "../../tmp/example_chr21.RData"
bim.file     <- "../../data/genotypes/example_chr21.bim"
qc.file      <- "../../data/qc/variants_qc.txt"
resolution   <- 2
out.basename <- "../../tmp/example_chr21"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
dendro.file  <- as.character(args[1])
bim.file     <- as.character(args[2])
qc.file      <- as.character(args[3])
resolution   <- round(as.numeric(args[4]))
out.basename <- as.character(args[5])

# Load libraries
suppressMessages(library(tidyverse))
local.source("utils.R")

# Load list of variants
Variants <- read_delim(bim.file, delim="\t", col_names=c("CHR", "SNP", "X", "BP", "A1", "A2"), col_types=cols())

# Remove variants that did not pass QC
variants.qc <- read_delim(qc.file, delim=" ", col_names="SNP", col_types=cols())
Variants <- Variants %>% inner_join(variants.qc, by="SNP")

# Load clustering dendrogram
cat("Loading clustering dendrogram... ")
load(dendro.file)
cat("done.\n")

# Choose groups by cutting the dendrogram
groups <- cutree(Sigma.clust, k = round(0.01*resolution*nrow(Variants)))
n.groups <- max(groups)
cat(sprintf("Divided %d variables into %d groups\n", length(groups), n.groups))

# Verify that the groups are nested
cat(sprintf("Checking whether the groups are nested... "))
are.nested <- is.dendrogram.nested(Sigma.clust, resolution.list=seq(1,100,length.out=10))
cat(sprintf("%s\n",are.nested))

# Show distribution of group sizes
cat(sprintf("Distribution of group sizes at resolution %d%%:\n", round(resolution)))
group.sizes <- sapply(1:n.groups, function(g) {sum(groups==g)} )
table(group.sizes)

# Write list of grouped variants
Variants.grouped <- mutate(Variants, Group = groups)
out.file <- sprintf("%s_groups%d.txt", out.basename, round(resolution))
write_delim(Variants.grouped, out.file)
cat(sprintf("Variant partition written on: %s\n", out.file))
