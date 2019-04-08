#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
args = commandArgs(trailingOnly=TRUE)

chr <- 6
clumping.factor <- 100
clumping.method <- "Radj"
K <- 50
n.load <- 1000

# Read input arguments
chr <- as.character(args[1])
clumping.factor <- as.numeric(args[2])
clumping.method <- as.character(args[3])
K <- as.numeric(args[4])
n.load <- as.numeric(args[5])

if(n.load == -1) n.load <- Inf

# Global parameters
block.size <- 1000  # Number of individuals to process simultaneously
numCores <- 1       # Use only 1 core because parallel computation is giving some problems (FIXME)

cat(paste("\nThis script is generating knockoffs for chromosome ",chr,".\n",
          "Running with K = ",K, ", clumping method ", clumping.method,
          " and clumping percentage = ", clumping.factor, ".\n", sep=""))

cat("Loading R libraries... ")
source("utils.R") # Load util functions
suppressMessages(library(tidyverse))
library(SNPknockG)
suppressMessages(library(snpStats))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
cat("done.\n")

# Make cluster for parallel computing
cl <- makeCluster(numCores)

# Specify location of data
pi.dir <- "/scratch/PI/candes"
rep.dir <- paste(pi.dir, "/ukbiobank_tmp/clumping/", clumping.method, clumping.factor, sep="")
fp.dir.in  <- paste(pi.dir, "/ukbiobank_tmp/fastphase/data", sep="")
fp.dir.out <- paste(pi.dir, "/ukbiobank_tmp/fastphase/phased", "_K", K, sep="")

# Specify location of output
knock.dir <- paste(pi.dir, "/ukbiobank_tmp/knockoffs/", clumping.method, clumping.factor, "_K", K, "/", sep="")
dir.create(knock.dir, showWarnings = FALSE)
knock.file.gen <- paste(knock.dir, "ukb_gen_chr", chr, ".ped", sep="")
knock.file.map <- paste(knock.dir, "ukb_gen_chr", chr, ".map", sep="")
knock.file.key <- paste(knock.dir, "ukb_gen_chr", chr, ".key", sep="")

# Load variable grouping
cat("Loading list of variants with clumping information... ")
grp.file <- paste(rep.dir, "/grp_chr", chr, ".txt", sep="")
Variants <- read_delim(grp.file, delim=" ", progress = FALSE, col_types = cols())
Variants$Group <- as.integer(Variants$Group)
cat("done.\n")

# Load HMM fitted by fastPhase
cat("Loading HMM fitted by fastPHASE... ")
alpha.file <- paste(fp.dir.out, "/ukb_hap_chr", chr, "_alphahat.txt", sep = "")
r.file <- paste(fp.dir.out, "/ukb_hap_chr", chr, "_rhat.txt", sep = "")
theta.file <- paste(fp.dir.out, "/ukb_hap_chr", chr, "_thetahat.txt", sep = "")
origchars.file <- paste(fp.dir.out, "/ukb_hap_chr", chr, "_origchars", sep = "")
hmm <- SNPknock.fp.loadFit(r.file, alpha.file, theta.file, origchars.file)
inp.file <- paste(fp.dir.in, "/ukb_hap_chr", chr, ".inp", sep = "")
hmm$positions <- scan(inp.file, what=integer(), skip=2, nlines=1, sep=" ", na.strings=c("P"), quiet=T)
hmm$positions <- hmm$positions[-1]
cat("done.\n")

# Transform HMM for subset of variants
cat("Transforming HMM for subset of variants... ")
hmm <- transform.hmm(hmm, Variants$Position)
cat("done.\n")

# Count number of subjects
n <- as.numeric(read_lines(inp.file, n_max=1))
if(n.load == Inf) {
  n.load <- n
}

# Randomly swap genotypes and knockoffs
cat("Generating random swaps... ")
n.groups <- max(Variants$Group)
group.swap <- sort(sample(n.groups, round(n.groups/2)))
Variants.X <- Variants %>% mutate(Swap = (Group %in% group.swap)) %>%
    mutate(Suffix = ifelse(Swap, "B", "A")) %>%
    mutate(Variant = paste(Variant, Suffix, sep=".")) %>%
    select(-Suffix)
Variants.Xk <- Variants %>% mutate(Swap = (Group %in% group.swap)) %>%
    mutate(Suffix = ifelse(Swap, "A", "B")) %>%
    mutate(Variant = paste(Variant, Suffix, sep=".")) %>%
    select(-Suffix)
# Save MAP file
Variants.A <- rbind(Variants.X, Variants.Xk) %>% select(Chr, Variant, X1, Position, Group, Swap) %>%
    mutate(Index=row_number()) %>% arrange(Position, Variant)
Variants.A %>% select(Chr, Variant, X1, Position) %>%
    write_delim(knock.file.map, delim=" ", col_names=FALSE, append=FALSE)
# Save swap key
Swap.key <- Variants.A %>%
    separate(Variant, c("Variant", "Suffix"), sep="\\.") %>%
    mutate(Variant = paste(Variant, Suffix, sep=".")) %>%
    mutate(Knockoff = as.logical((Suffix=="A")*Swap + (Suffix=="B")*(1-Swap))) %>%
    select(Chr,Variant,Position,Group,Knockoff,Swap)
write_delim(Swap.key, knock.file.key, delim=" ")
cat("done.\n")

# Generate knockoffs in blocks
cat(paste("Generating knockoff copies for ", n.load, " subjects and ", nrow(Variants), " variants.\n",sep=""))
blocks <- split(1:n.load, ceiling(seq_along(1:n.load)/block.size))

for(block in blocks) {

    block.start <- block[1]
    block.length <- length(block)

    # Load haplotypes for this block
    cat(paste("\nLoading haplotypes for individuals ",block.start,"through",block.start+block.length-1,"... "))
    basename <- paste(fp.dir.in, "/ukb_hap_chr", chr, sep="")
    H <- read.inp(basename, n.skip=block.start-1, n.load=block.length, progress=FALSE)
    H.subjects <- rownames(H)
    keep.cols <- which(colnames(H) %in% Variants$Variant)
    H <- as.matrix(H[,keep.cols])
    cat("done.\n")

    # Generate knockoff haplotypes
    cat(sprintf("Generating knockoff copies on %d parallel threads... ", numCores))
    Hk <- SNPknock.knockoffHaplotypes_group(H, hmm$r, hmm$alpha, hmm$theta, Variants$Group, cluster=cl,
                                            seed=block.start, display_progress=FALSE)
    cat("done.\n")

    # Converting haplotypes to genotypes
    X  <- H[seq(1,nrow(H),2),] + H[seq(2,nrow(H),2),]
    Xk <- Hk[seq(1,nrow(Hk),2),] + Hk[seq(2,nrow(Hk),2),]

    # Write genotypes in PED format
    cat("Writing augmented genotypes... ")
    colnames(X)  <- Variants.X$Variant
    colnames(Xk) <- Variants.Xk$Variant
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

stopCluster(cl)
