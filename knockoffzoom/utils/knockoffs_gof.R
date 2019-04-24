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
stats.basename <- "../../tmp/example_chr21_res2"
key.file       <- "../../tmp/example_chr21_res2.key"
groups.file    <- "../../tmp/example_chr21_groups2.txt"
out.basename   <- "../../results/example_chr21_res2"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
stats.basename <- as.character(args[1])
key.file       <- as.character(args[2])
groups.file    <- as.character(args[3])
out.basename   <- as.character(args[4])

# Load libraries
cat("Loading R libraries... ")
suppressMessages(library(tidyverse))
cat("done.\n")

# Load variable grouping
cat("Loading list of variants with clumping information... ")
Variants <- read_delim(groups.file, delim=" ", progress = FALSE, col_types = cols()) %>%
    mutate(Group = as.integer(Group))
cat("done.\n")

# Load knockoff keys
Key <- read_delim(key.file, delim=" ", col_types=cols())

# Load variant frequency table
frq.file <- sprintf("%s.frq", stats.basename)
Frq <- read_table(frq.file, col_types=cols())
# Compute diagnostics
maf <- Frq$MAF
maf.x  <- maf[seq(1,length(maf), by=2)]
maf.xk <- maf[seq(2,length(maf), by=2)]
Diagnostics <- tibble(index=seq(length(maf.x)), x=maf.x, xk=maf.xk) %>%
    mutate(error=(x-xk)^2)

# Plot frequency diagnostics
p.frq <- Diagnostics %>%
    ggplot(aes(x=x, y=xk)) +
    geom_point(alpha=0.2) +
    geom_abline(intercept = 0, slope = 1, color="red", linetype=2) +
    theme_bw()
# Save frequency plot
frq.pfile <- sprintf("%s_frq.png", out.basename)
ggsave(frq.pfile, plot=p.frq, width = 4, height = 4, dpi = 300, units = c("in", "cm", "mm"))
cat(sprintf("Plot saved on: %s\n", frq.pfile))

# Load LD table
ld.file <- sprintf("%s.ld", stats.basename)
LD <- read_table2(ld.file, col_types=cols()) %>%
    mutate(CHR=CHR_A) %>% select(-CHR_A, -CHR_B)

# Add grouping information
LD <- LD %>%
    mutate(SNP=SNP_A) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP)) %>%
    inner_join(Variants %>% select(CHR, SNP, Group), by = c("CHR", "SNP")) %>%
    mutate(Group_A=Group) %>% select(-SNP, -Group) %>%
    mutate(SNP=SNP_B) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP)) %>%
    inner_join(Variants %>% select(CHR, SNP, Group), by = c("CHR", "SNP")) %>%
    mutate(Group_B=Group) %>% select(-SNP, -Group)

# Add knockoff key information
LD <- LD %>%
    mutate(SNP=SNP_A) %>%
    inner_join(Key %>% select(CHR, SNP, Knockoff), by = c("CHR", "SNP")) %>%
    mutate(Knockoff_A=Knockoff) %>% select(-SNP, -Knockoff) %>%
    mutate(SNP=SNP_B) %>%
    inner_join(Key %>% select(CHR, SNP, Knockoff), by = c("CHR", "SNP")) %>%
    mutate(Knockoff_B=Knockoff) %>% select(-SNP, -Knockoff)

# Create correlation tables between different groups
group.range <- seq(0,10)
LD.XX <- LD %>%
    filter(abs(Group_B-Group_A) %in% group.range, Knockoff_A==FALSE, Knockoff_B==FALSE) %>%
    mutate(R.XX=R2) %>%
    mutate(SNP_A=str_replace(SNP_A,".A",""), SNP_A=str_replace(SNP_A,".B","")) %>%
    mutate(SNP_B=str_replace(SNP_B,".A",""), SNP_B=str_replace(SNP_B,".B","")) %>%
    select(Group_A, Group_B, SNP_A, SNP_B, R.XX) %>%
    distinct(Group_A, Group_B, SNP_A, SNP_B, R.XX)
LD.XkXk <- LD %>%
    filter(abs(Group_B-Group_A) %in% group.range, Knockoff_A==TRUE, Knockoff_B==TRUE) %>%
    mutate(R.XkXk=R2) %>%
    mutate(SNP_A=str_replace(SNP_A,".A",""), SNP_A=str_replace(SNP_A,".B","")) %>%
    mutate(SNP_B=str_replace(SNP_B,".A",""), SNP_B=str_replace(SNP_B,".B","")) %>%
    select(Group_A, Group_B, SNP_A, SNP_B, R.XkXk) %>%
    distinct(Group_A, Group_B, SNP_A, SNP_B, R.XkXk)
LD.XXk <- LD %>%
    filter((Group_B-Group_A) %in% seq(1,10), Knockoff_A*Knockoff_B==FALSE) %>%
    mutate(R.XXk=R2) %>%
    mutate(SNP_A=str_replace(SNP_A,".A",""), SNP_A=str_replace(SNP_A,".B","")) %>%
    mutate(SNP_B=str_replace(SNP_B,".A",""), SNP_B=str_replace(SNP_B,".B","")) %>%
    select(Group_A, Group_B, SNP_A, SNP_B, R.XXk) %>%
    distinct(Group_A, Group_B, SNP_A, SNP_B, R.XXk)

# Plot originality
LD.cross <- inner_join(LD.XX, LD.XkXk, by = c("Group_A", "Group_B", "SNP_A", "SNP_B"))
p.orig <- LD.cross %>%
    mutate(Distance = as.factor(abs(Group_A-Group_B))) %>%
    ggplot(aes(x=abs(R.XX), y=abs(R.XkXk))) +
    geom_abline(color="red") +
    geom_point(alpha=0.1) +
    xlim(0,1) + ylim(0,1) +
    theme_bw()
p.orig

# Save plot
orig.pfile <- sprintf("%s_orig.png", out.basename)
ggsave(orig.pfile, plot=p.orig, width = 4, height = 4, dpi = 300, units = c("in", "cm", "mm"))
cat(sprintf("Plot saved on: %s\n", orig.pfile))

# Plot exchangeability
options(repr.plot.width=4, repr.plot.height=3)
LD.cross <- inner_join(LD.XX, LD.XXk, by = c("Group_A", "Group_B", "SNP_A", "SNP_B")) %>%
    filter(Group_A!=Group_B)
p.exch <- LD.cross %>%
    mutate(Distance = as.factor(abs(Group_A-Group_B))) %>%
    ggplot(aes(x=abs(R.XX), y=abs(R.XXk))) +
    geom_abline(color="red") +
    geom_point(alpha=0.1) +
    xlim(0,1) + ylim(0,1) +
    theme_bw()
p.exch

#Save plot
exch.pfile <- sprintf("%s_exch.png", out.basename)
ggsave(exch.pfile, plot=p.exch, width = 4, height = 4, dpi = 300, units = c("in", "cm", "mm"))
cat(sprintf("Plot saved on: %s\n", exch.pfile))

# Plot histogram of self-correlations
p.self <- LD %>% filter(BP_A==BP_B) %>%
    ggplot(aes(x=R2)) +
    geom_histogram(bins=30) +
    theme_bw()
p.self

# Save plot
self.pfile <- sprintf("%s_self.png", out.basename)
ggsave(self.pfile, plot=p.self, width = 4, height = 4, dpi = 300, units = c("in", "cm", "mm"))
cat(sprintf("Plot saved on: %s\n", self.pfile))
