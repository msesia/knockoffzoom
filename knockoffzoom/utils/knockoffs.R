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
hmm.basename <- "../../tmp/example_chr21"
inp.basename <- "../../tmp/example_chr21"
groups.file  <- "../../tmp/example_chr21_groups2.txt"
out.basename <- "../../tmp/example_chr21"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
hmm.basename <- as.character(args[1])
inp.basename <- as.character(args[2])
groups.file  <- as.character(args[3])
out.basename <- as.character(args[4])

# Global parameters
block.size <- 1000  # Number of individuals to process simultaneously
numCores <- 1       # Use only 1 core because parallel computation is giving some problems (FIXME)

# Load libraries
cat("Loading R libraries... ")
local.source("utils.R")
suppressMessages(library(tidyverse))
suppressMessages(library(snpStats))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
cat("done.\n")

# Make cluster for parallel computing
cl <- makeCluster(numCores)

# Inp haplotype file
inp.file <- sprintf("%s.inp", inp.basename)

# Load HMM fitted by fastPhase
cat("Loading HMM fitted by fastPHASE... ")
r.file <- sprintf("%s_rhat.txt", hmm.basename)
alpha.file <- sprintf("%s_alphahat.txt", hmm.basename)
theta.file <- sprintf("%s_thetahat.txt", hmm.basename)
origchars.file <- sprintf("%s_origchars", hmm.basename)
hmm <- SNPknock::loadHMM(r.file, alpha.file, theta.file, origchars.file)
hmm$positions <- scan(inp.file, what=integer(), skip=2, nlines=1, sep=" ", na.strings=c("P"), quiet=T)
hmm$positions <- hmm$positions[-1]
cat("done.\n")

# Load variable grouping
cat("Loading list of variants with clumping information... ")
Variants <- read_delim(groups.file, delim=" ", progress = FALSE, col_types = cols()) %>%
    mutate(Group = as.integer(Group))
cat("done.\n")

# Transform HMM for subset of variants
cat("Transforming HMM for subset of variants... ")
hmm <- transform.hmm(hmm, Variants$BP)
cat("done.\n")

# Specify location of output
knock.file.gen <- sprintf("%s.ped", out.basename)
knock.file.map <- sprintf("%s.map", out.basename)
knock.file.key <- sprintf("%s.key", out.basename)

# Count number of subjects
n <- as.numeric(read_lines(inp.file, n_max=1))

# Set random seed
seed.value <- as.integer(Variants$CHR[1])
set.seed(seed.value)
cat(sprintf("Random seed: %d\n", seed.value))

# Randomly swap genotypes and knockoffs
cat("Generating random swaps... ")
n.groups <- max(Variants$Group)
group.swap <- sort(sample(n.groups, round(n.groups/2)))
Variants.X <- Variants %>% mutate(Swap = (Group %in% group.swap)) %>%
    mutate(Suffix = ifelse(Swap, "B", "A")) %>%
    mutate(SNP = paste(SNP, Suffix, sep=".")) %>%
    select(-Suffix)
Variants.Xk <- Variants %>% mutate(Swap = (Group %in% group.swap)) %>%
    mutate(Suffix = ifelse(Swap, "A", "B")) %>%
    mutate(SNP = paste(SNP, Suffix, sep=".")) %>%
    select(-Suffix)
# Save MAP file
Variants.A <- rbind(Variants.X, Variants.Xk) %>% select(CHR, SNP, X, BP, Group, Swap) %>%
    mutate(Index=row_number()) %>% arrange(BP, SNP)
Variants.A %>% select(CHR, SNP, X, BP) %>%
    write_delim(knock.file.map, delim=" ", col_names=FALSE, append=FALSE)
# Save swap key
Swap.key <- Variants.A %>%
    separate(SNP, c("SNP", "Suffix"), sep="\\.") %>%
    mutate(SNP = paste(SNP, Suffix, sep=".")) %>%
    mutate(Knockoff = as.logical((Suffix=="A")*Swap + (Suffix=="B")*(1-Swap))) %>%
    select(CHR,SNP,BP,Group,Knockoff,Swap)
write_delim(Swap.key, knock.file.key, delim=" ")
cat("done.\n")

# Generate knockoffs in blocks
cat(paste("Generating knockoff copies for ", n, " subjects and ", nrow(Variants), " variants.\n",sep=""))
blocks <- split(1:n, ceiling(seq_along(1:n)/block.size))

for(block in blocks) {
    block.start <- block[1]
    block.length <- length(block)

    # Load haplotypes for this block
    cat(paste("\nLoading haplotypes for individuals ",block.start,"through",block.start+block.length-1,"... "))
    H <- read.inp(inp.basename, n.skip=block.start-1, n.load=block.length, progress=FALSE)
    H.subjects <- rownames(H)
    keep.cols <- which(colnames(H) %in% Variants$SNP)
    H <- as.matrix(H[,keep.cols])
    cat("done.\n")

    # Generate knockoff haplotypes
    cat(sprintf("Generating knockoff copies on %d parallel threads... ", numCores))
    Hk <- SNPknock::knockoffHaplotypes(H, hmm$r, hmm$alpha, hmm$theta, groups=Variants$Group, cluster=cl,
                                       seed=block.start, display_progress=FALSE)
    cat("done.\n")

    # Converting haplotypes to genotypes
    X  <- H[seq(1,nrow(H),2),] + H[seq(2,nrow(H),2),]
    Xk <- Hk[seq(1,nrow(Hk),2),] + Hk[seq(2,nrow(Hk),2),]

    # Write genotypes in PED format
    cat("Writing augmented genotypes... ")
    colnames(X)  <- Variants.X$SNP
    colnames(Xk) <- Variants.Xk$SNP
    rownames(X)  <- H.subjects[seq(1,length(H.subjects),2)]
    rownames(Xk) <- rownames(X)
    X.all <- cbind(X,Xk)
    X.all <- X.all[,Variants.A$Index]
    A.SnpMatrix <- as(X.all, "SnpMatrix")
    write.SnpMatrix(A.SnpMatrix,knock.file.gen, as.alleles=TRUE, quote=FALSE, col.names=FALSE,
                    append=(block.start!=1))
    cat("done.\n")
}
cat("\n")

cat(sprintf("Done. Results written on: %s\n", knock.file.gen))
stopCluster(cl)


if(FALSE) {

    # Compare SNP column means
    plot(colMeans(X),colMeans(Xk),col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

    # Compare correlations between consecutive SNPs
    corrX = unlist(lapply(seq(1,10), function(dj) { sapply((1+dj):ncol(X), function(j) cor(X[,j-dj],X[,j])) }))
    corrXk = unlist(lapply(seq(1,10), function(dj) { sapply((1+dj):ncol(X), function(j) cor(Xk[,j-dj],Xk[,j])) }))
    plot(corrX,corrXk,col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)
    plot(corrX^2,corrXk^2,col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

}
