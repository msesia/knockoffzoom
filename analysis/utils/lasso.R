#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Install packages bigsnpr and bigstatsr
# devtools::install_github("privefl/bigstatsr")
# devtools::install_github("privefl/bigsnpr")

# Documentation here:
# https://privefl.github.io/bigsnpr/reference/index.html
# https://privefl.github.io/bigstatsr/reference/big_spLinReg.html

# Load packages
library(tidyverse)
library(devtools)
library(bigsnpr)

# Default arguments
phenotype <- "respiratory"
resolution <- "Radj2"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
phenotype <- as.character(args[1])
resolution <- as.character(args[2])

# Other parameters
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1,22)
R2.max <- 0.99
fit.lasso <- TRUE
use.PCs <- TRUE

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    # Load list of variants
    key.file <- sprintf("%s/knockoffs/%s_K50/ukb_gen_chr%d.key", scratch, resolution, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    # Load LD table
    ld.file <- sprintf("%s/knockoff_diagnostics/%s_K50/ukb_gen_chr%d.ld", scratch, resolution, chr)
    LD.chr <- read_table(ld.file, col_types=cols(), guess_max=Inf) %>%
        filter(BP_A==BP_B) %>%
        mutate(CHR=CHR_A, BP=BP_A) %>%
        select(CHR, BP, R2)
    # Combine list of variants with LD and MAF tables
    Variants.chr <- Variants.chr %>% left_join(LD.chr, by = c("CHR", "BP"))
    cat("done.\n")
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

####################
## Load genotypes ##
####################

rds.file <- sprintf("%s/augmented_data_big/ukb_gen_%s.rds", scratch, resolution)
if(file.exists(rds.file)){
    cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
    cat(sprintf("Could not find FBM in %s.\n", rds.file))
    quit()
}

# Attach the "bigSNP" object in R session
ptm <- proc.time() # Start the clock!
cat("Attaching bigSNP object... ")
obj.bigSNP <- snp_attach(rds.file)
cat("done.\n")
proc.time() - ptm # Stop the clock

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- obj.bigSNP$map$marker.ID

# Compute scaling factor for the genotypes
scale.file <- sprintf("%s/augmented_data_big/ukb_gen_%s_scale.txt", scratch, resolution)
if(!file.exists(scale.file)){
    ptm <- proc.time() # Start the clock!
    cat("Computing scaling factors for all variants... ")
    scaler <- big_scale()
    G.scale <- scaler(G)
    scaling.factors <- G.scale$scale
    cat("done.\n")
    proc.time() - ptm # Stop the clock
    # Save the scaling factors to file
    write.table(scaling.factors, scale.file, sep = "\t", row.names=F, col.names=F)
    cat(sprintf("Saved scaling factors to %s\n", scale.file))
} else {
    cat(sprintf("Loading scaling factors from %s... ", scale.file))
    scaling.factors <- read.table(scale.file)$V1
    cat("done.\n")
}

if(fit.lasso==TRUE) {

    #####################
    ## Load phenotypes ##
    #####################
    cat("Reading phenotype file... ")

    # Load list of subjects
    fam.file <- sprintf("%s/QC_output/individuals_QC.txt", scratch)
    Subjects <- read_tsv(fam.file, col_types=cols(), col_names = c("FID", "IID"))

    # Load response values
    pheno.file <- sprintf("%s/phenotypes/phenotypes_qc.tab", scratch)
    Phenotypes <- read_delim(pheno.file, delim="\t", col_types=cols())
    Phenotypes <- Phenotypes %>% right_join(Subjects, by=c("FID", "IID"))

    # Make sure that the rows of the genotypes match the rows of the phenotypes
    Phenotypes <- Phenotypes %>%
        right_join(transmute(obj.bigSNP$fam, FID=family.ID, IID=sample.ID), by = c("FID", "IID"))

    cat("done.\n")

    ###################
    ## Fit the lasso ##
    ###################

    # Extract response variable
    ind.train <- which(!is.na(Phenotypes[[phenotype]]))
    cat(sprintf("%d individuals out of %d have missing phenotype.\n",
                nrow(Subjects)-length(ind.train), nrow(Subjects)))
    y <- Phenotypes[[phenotype]][ind.train]

    # Extract covariates
    covariate.file <- sprintf("%s/phenotypes/analysis.tab", scratch)
    Analysis <- read_tsv(covariate.file, col_types=cols())
    covariate.names <- Analysis %>% filter(Name==phenotype) %>% select(Covariates) %>%
        as.character() %>% strsplit(",")
    if(use.PCs) {
        covariate.names <- c(covariate.names[[1]], paste("PC", seq(5), sep="."))
    } else {
        covariate.names <- covariate.names[[1]]
    }
    Covariates <- Phenotypes %>% select(covariate.names)
    covar.train <- as.matrix(Covariates)[ind.train,]

    # Find the class of response (numeric or factor)
    phenotype.class <- Analysis %>% filter(Name==phenotype) %>% select(Class) %>% as.character()
    if(phenotype.class=="factor") {
        y <- factor(y, levels=c(1,2), labels=c(0,1))
        y <- as.numeric(levels(y))[y]
    }

    # Fit the lasso
    ptm <- proc.time() # Start the clock!
    # Use only linear regression (debug)
    #phenotype.class <- "numeric"
    # Set dfmax (may want to increase to 20k)
    dfmax <- 10000
    if(phenotype.class=="factor") {
        cat(sprintf("Fitting sparse logistic regression with %d observations, %d variants and %d covariates... ",
                    length(y), ncol(G), ncol(covar.train)))
        lasso.fit <- big_spLogReg(G, y01.train=y, ind.train=ind.train, covar.train=covar.train,
                                  dfmax=dfmax, ncores=10)
    } else {
        cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and %d covariates... ",
                    length(y), ncol(G), ncol(covar.train)))
        lasso.fit <- big_spLinReg(G, y.train=y, ind.train=ind.train, covar.train=covar.train,
                                  dfmax=dfmax, ncores=10)
    }
    cat("done.\n")
    proc.time() - ptm # Stop the clock

    # Extract beta from each fold and combine them
    cat("Extracting importance measures... ")
    beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)

    # Saving beta matrix
    if(use.PCs) {
        beta.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso_beta.txt", scratch, phenotype, resolution)
        beta %>% as_tibble() %>% write_delim(beta.file, delim=" ")
        cat(sprintf("Lasso coefficient matrix (%d x %d) saved to:\n%s\n", nrow(beta), ncol(beta), beta.file))
    }
} else {
    if(use.PCs) {
        beta.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso_beta.txt", scratch, phenotype, resolution)
        beta <- read_delim(beta.file, delim=" ") %>% as.matrix()
    } else {
        error("The lasso coefficient matrix without PCs was not stored.\n")
    }
}

# Separate coefficients for variants from coefficients for covariates
beta.variants <- beta[1:ncol(G),]
beta.covariates <- beta[(ncol(G)+1):nrow(beta),]

# Undo scaling of lasso coefficients
beta.variants <- beta.variants * scaling.factors
Beta <- cbind(tibble(CHR=Variants$CHR,
                     SNP=Variants$SNP, BP=Variants$BP),
                     as_tibble(beta.variants)) %>% as_tibble()
colnames(Beta) <- c("CHR", "SNP", "BP", paste("K", seq(ncol(beta.variants)),sep=""))
Beta <- Beta %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) #%>%
#    filter(Nonzero==10)
# Extract the estimated coefficients
Lasso.res <- Beta %>% mutate(Resolution=resolution) %>%
    inner_join(Variants, by = c("CHR", "SNP", "BP")) %>%
    filter(Z!=0) %>%
    arrange(desc(Z))
cat("done.\n")

# Print summary of fit
Lasso.res %>% summarise(Nonzero=sum(Z!=0))

##########################
## Save lasso estimates ##
##########################
if(use.PCs) {
    out.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso.txt", scratch, phenotype, resolution)
} else {
    out.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso_nopc.txt", scratch, phenotype, resolution)
}        
Lasso.res %>% mutate(Importance=abs(Z)) %>%
    select(Resolution, CHR, SNP, BP, Group, Importance) %>%
    write_delim(out.file, delim=" ")
cat(sprintf("Results saved to %s.\n", out.file))

##################################
## Assemble knockoff statistics ##
##################################
cat("Computing knockoff statistics... ")

# Compute the knockoff statistics
W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
}

Stats <- Lasso.res %>% select("Resolution", "CHR", "Group", "SNP", "BP", "Z") %>%
    filter(Z!=0) %>% 
    left_join(Variants, by = c("CHR", "Group", "SNP", "BP")) %>%
    group_by(Resolution, CHR, Group) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n(), R2=max(R2)) %>%
    mutate(SNP.lead = gsub(".A", "", SNP.lead), SNP.lead = gsub(".B", "", SNP.lead)) %>%
    ungroup() %>%
    arrange(desc(abs(W))) %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W, R2, Resolution) %>%
    filter(W!=0)
cat("done.\n")
# Show a preview of the knockoff stats
cat("Knockoff statistics:\n")
Stats %>% print(n=10)

if(use.PCs) {
    stats.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso_stats.txt", scratch, phenotype, resolution)
} else {
    stats.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso_stats_nopc.txt", scratch, phenotype, resolution)
}        
Stats %>% write_delim(stats.file, delim=" ")
cat(sprintf("Knockoff statistics saved to:\n%s\n", stats.file))

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

Selections <- Stats %>% filter(is.na(R2)|R2<=R2.max) %>% knockoff.filter(fdr=0.1, offset=1)

Selections <- Selections %>%
    select(Resolution, CHR, Group, SNP.lead, BP.lead, W) %>% inner_join(Variants, by = c("CHR", "Group")) %>%
    group_by(Resolution, CHR, Group, SNP.lead, BP.lead) %>%
    summarise(W=mean(W), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n()/2,) %>%
    ungroup() %>%
    mutate(Method="Knockoffs") %>%
    arrange(desc(W))

# Give a preview of the selections
cat("Knockoff selections:\n")
Selections %>% print()

cat("Number of discoveries made by knockoffs:\n")
Selections %>% summarise(Discoveries=n())

# Save list of discoveries
if(use.PCs) {
    discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s.txt", scratch, phenotype, resolution)
} else {
    discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s_nopc.txt", scratch, phenotype, resolution)
}        
Selections %>% select(CHR, SNP.lead, BP.lead, W, BP.min, BP.max, BP.width, Size, Group) %>%
    write_delim(discoveries.file, delim=" ")
cat(sprintf("Saved list of %d discoveries to %s\n", nrow(Selections), discoveries.file))
