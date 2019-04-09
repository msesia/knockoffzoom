#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
chr = as.integer(args[1])
partitions.factor = as.numeric(args[2])/100
partitions.method <- as.character(args[3])
    
suppressMessages(library(adjclust))
suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
source("utils.R")

# Parse input
ld.measure <- substring(partitions.method, 1, 1)
clustering <- substring(partitions.method, 2)

# Paths to data
dat.dir <- "/scratch/PI/candes/ukbiobank"
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

# Load list of variants for phased haplotypes that have passed QC
positions.filename <- paste(scratch, "/fastphase/data/ukb_hap_chr", chr, ".legend", sep="")
variant.names <- read_delim(positions.filename, delim=" ", skip=1, col_types=cols(),
                           col_names=c("Variant", "Position", "A1", "A2"))
bim.filename <- paste(dat.dir, "/genotypes/ukb_gen_chr", chr, ".bim", sep="")
Variants <- as_tibble(read.table(bim.filename, sep="", header=FALSE, stringsAsFactors=FALSE,
                                col.names=c("Chr", "Variant", "X1", "Position", "A1", "A2")))
Variants <- semi_join(Variants, variant.names, by=c("Variant"))

# Load MAF
frq.filename <- sprintf("%s/stats/ukb_gen_chr%s.frq", scratch, chr)
Frequencies <- as_tibble(read.table(frq.filename, sep="", header=TRUE, stringsAsFactors=FALSE,
                                   col.names=c("Chr","Variant","A1","A2","MAF","NCHROBS")))
Variants <- inner_join(Variants, Frequencies %>% select(c("Variant", "MAF", "NCHROBS")), by=c("Variant"))

# Compute dendrogram if not already present
dend.filename <- sprintf("%s/partitions/%s/grp_chr%s.RData", scratch, partitions.method, chr)
if(!file.exists(dend.filename)) {
    cat(sprintf("File %s not found.\n", dend.filename))
    cat("Loading covariance matrix from... ")   
    # Load sparse covariance matrix
    corr.filename <- sprintf("%s/stats/ukb_gen_chr%s.ld", scratch, chr)
    LD <- read_table2(corr.filename, col_types=cols())
    LD <- filter(LD, BP_A %in% Variants$Position, BP_B %in% Variants$Position)
    cat("done.\n")

    cat("Converting covariance matrix... ")   
    # Convert LD to a symmetric square matrix
    LD.variants <- unique(c(LD$SNP_A, LD$SNP_B))
    LD.positions <- unique(c(LD$BP_A, LD$BP_B))
    Sigma <- ld.to.mat(LD, ld.measure, Variants$Position)

    # Convert covariance to distance matrix
    #Sigma.dist <- as.dist(1-abs(Sigma))                   # Faster, but requires lots of memory
    Sigma.dist <- mat.to.dist(tril(Sigma,-1))            # Slower, but requires less memory
    cat("done.\n")
    
    cat("Computing clustering dendrogram... ")   
    if(clustering == "single") { 
        # Clustering
        Sigma.clust <- hclust(Sigma.dist, method = "single")  # Does not make contiguous groups
    } else if(clustering == "adj") {
        if(nrow(Sigma) > 10000) {
            Sigma.clust <- adjClust(Sigma.dist,h=10000)       # Makes contiguous groups
        } else {
            Sigma.clust <- adjClust(Sigma.dist)               # Makes contiguous groups
        }
    } else {
        stop(sprintf("Error: unknown clustering method %s\n", clustering))
    }

    # Save dendrogram
    save(list=c("Sigma.clust"), file=dend.filename)
    cat("done.\n")
} else {
    cat(sprintf("File %s found.\n", dend.filename))
    cat("Loading clustering dendrogram... ")
    load(dend.filename)
    cat("done.\n")
}

# Choose groups by cutting the dendrogram
groups <- cutree(Sigma.clust, k = round(partitions.factor*nrow(Variants)))
n.groups <- max(groups)
cat(sprintf("Divided %d variables into %d groups\n", length(groups), n.groups))

# Verify that the groups are nested
cat(sprintf("Checking whether the groups are nested... "))
are.nested <- is.dendrogram.nested(Sigma.clust, partitions.list=seq(1,100,length.out=10))
cat(sprintf("%s\n",are.nested))

# Show distribution of group sizes
cat(sprintf("Distribution of group sizes at partitions %.4f:\n",partitions.factor))
group.sizes <- sapply(1:n.groups, function(g) {sum(groups==g)} )
table(group.sizes)

# Pick representatives based on MAF
group.repres <- sapply(1:n.groups, function(g) {
    maf <- rep(0, nrow(Variants))
    maf[which(groups==g)] <- Variants$MAF[which(groups==g)]
    which.max(maf)
})

# Look at remaining correlations
#Sigma.small <- Sigma[group.repres,group.repres]
#cat(sprintf("Largest remaining correlation: %s\n", max(abs(Sigma.small))))

# Compute outer correlations
## outcorr <- sapply(1:max(groups), function(g) {
##     g.idx <- which(groups==g)
##     max(abs(Sigma[,-g.idx][g.idx,]))
## })
## cat(sprintf("Typical remaining correlation: %s\n", round(mean(outcorr),3)))

# Create table of grouped variants
#Variants.grouped <- mutate(Variants, Group = groups, Outcorr = outcorr[groups])
Variants.grouped <- mutate(Variants, Group = groups)
grp.filename <- sprintf("%s/partitions/%s%s/grp_chr%s.txt",
                        scratch, partitions.method, partitions.factor*100, chr)
write_delim(Variants.grouped, grp.filename)

# Create table of group representatives
Variants.rep <- Variants[group.repres,]
Variants.rep <- mutate(Variants.rep, Size=group.sizes)
Variants.rep <- arrange(Variants.rep, Position)
rep.filename <- sprintf("%s/partitions/%s%s/rep_chr%s.txt",
                        scratch, partitions.method, partitions.factor*100, chr)
write_delim(Variants.rep, rep.filename)
