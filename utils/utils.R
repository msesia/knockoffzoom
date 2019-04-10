suppressMessages(library(adjclust))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))

read.inp <- function(basename, n.skip=0, n.load=Inf, progress=FALSE) {
    # Read haplotype sequences from INP file
    H.file <- paste(basename, ".inp", sep="")
    # Read header of INP file    
    n <- as.numeric(read_lines(H.file, n_max=1, skip=0))
    p <- as.numeric(read_lines(H.file, n_max=1, skip=1))
    # Check whether file contains variant positions
    positions <- scan(H.file, what=character(), skip=2, nlines=1, sep=" ", quiet=TRUE)
    if(startsWith(positions[1], "P")) {
        positions <- scan(H.file, what=integer(), skip=2, nlines=1, sep=" ", na.strings=c("P"), quiet=TRUE)[-1]
        skip.pos <- 1
    }
    else {
        positions <- rep(NA, p)
        skip.pos <- 0
    }
    
    if(n.load == Inf) {
        n.load <- n
    }
    # Read body of INP file (fast)
    H <- fread(H.file, skip=2+skip.pos+3*n.skip+1, nrows=3*n.load-1,
               header=F, na.strings="?", showProgress=progress) %>%
        as_tibble() %>% filter(substr(V1,1,1)!="#")
    H <- matrix(unlist(strsplit(H$V1,'')) %>% as.integer(), p) %>% t()
   
    # Read variant legend
    legend.file <- paste(basename, ".legend", sep="")
    if (file.exists(legend.file)) {
        Legend <- read_delim(legend.file, delim=" ",
                             col_types=cols(position = col_integer(), .default = col_character()))
        colnames(Legend) <- c("Variant", "Position", "A1", "A2")
        stopifnot(all.equal(Legend$Position,positions))
        colnames(H) <- Legend$Variant
    }
    # Read subject IDs
    sample.file <- paste(basename, ".sample", sep="")
    if (file.exists(sample.file)) {
        Subjects <- read_delim(sample.file, delim=" ", skip=2+n.skip, n_max=n.load, col_names=FALSE,
                               col_types=cols(.default = col_character()))
        colnames(Subjects) <- c("Pedigree","Family","Missing","Sex")
        rownames(H) <- rep(Subjects$Pedigree, each=2)
    }
    return(H)
}

read.haps <- function(basename, progress=FALSE) {
    # Read haplotype sequences from HAPS file
    H.file <- paste(basename, ".haps", sep="")
    H <- fread(H.file, sep=" ",showProgress=progress)
    # Read variant legend
    legend.file <- paste(basename, ".legend", sep="")
    Legend <- read_delim(legend.file, delim=" ",
                         col_types=cols(position = col_integer(), .default = col_character()))
    colnames(Legend) <- c("Variant", "Position", "A1", "A2")
    # Read subject IDs
    sample.file <- paste(basename, ".sample", sep="")
    Subjects <- read_delim(sample.file, delim=" ", skip=2, col_names=FALSE,
                           col_types=cols(X4 = col_character(), .default = col_integer()))
    colnames(Subjects) <- c("Pedigree","Family","Missing","Sex")
    # Transpose H and use variant IDs as column names
    H <- as_tibble(t(H))
    colnames(H) <- Legend$Variant
    # Add column with individual IDs
    H$Subject <- rep(Subjects$Pedigree, each=2)
    H <- dplyr::select(H, Subject, everything())
    return(H)
}

read.raw <- function(basename, n.skip=0, n.load=Inf, progress=FALSE) {
    # Read genotype sequences
    X.file <- paste(basename, ".raw", sep="")
    X <- read_delim(X.file, delim="\t", skip=n.skip, n_max=n.load,
                    col_types = cols(.default = col_integer()), progress=progress)
    variants <- colnames(X)[7:ncol(X)]
    variants <- gsub("\\_.*", "", variants)
    # Transform X into same format as output of read.haps
    X <- dplyr::select(X, -c(1,3,4,5,6))
    colnames(X) <- c("Subject", variants)
    return(X)
}

transform.hmm <- function(hmm, positions) {
    select.snps <- which(hmm$positions %in% positions)
    select.snps <- select.snps[order(hmm$positions[select.snps])]
    hmm.repr <- c()
    hmm.repr$alpha <- hmm$alpha[select.snps,]
    hmm.repr$theta <- hmm$theta[select.snps,]    
    hmm.repr$r <- sapply(1:length(select.snps), function(j) {
        if(j==1) {
            return(hmm$r[1])
        } else {
            return(sum(hmm$r[(select.snps[j-1]+1):select.snps[j]]))
        }
    })
    hmm.repr$positions <- positions
    return(hmm.repr)
}

ld.to.mat <- function (LD, ld.measure, positions) {    
    stopifnot(ld.measure %in% c("R", "D"))

    stopifnot(length(unique(LD$CHR_A))<=1)
    stopifnot(length(unique(LD$CHR_A))<=1)
    chr <- LD$CHR_A[1]
    
    if(ld.measure == "R") {
        LD <- add_row(LD, R2 = 1,
                      CHR_A=chr, SNP_A=NA, BP_A = tail(positions,1),
                      CHR_B=chr, SNP_B=NA, BP_B = tail(positions,1))
        Sigma = sparseMatrix(i = match(LD$BP_A,positions),
                             j = match(LD$BP_B,positions),
                             x = LD$R2, symmetric=TRUE)
    } else {
        LD <- add_row(LD, DP = 1,
                      CHR_A=chr, SNP_A=NA, BP_A = tail(positions,1),
                      CHR_B=chr, SNP_B=NA, BP_B = tail(positions,1))
        Sigma = sparseMatrix(i = match(LD$BP_A,positions),
                             j = match(LD$BP_B,positions),
                             x = LD$DP, symmetric=TRUE)

    }
    diag(Sigma) = rep(1,ncol(Sigma))
    colnames(Sigma) = positions
    rownames(Sigma) = positions
    return(Sigma)
}

mat.to.dist <- function(S, progress=FALSE) {
  S <- tril(S,-1)
  p.small <- 10
  if(nrow(S)<p.small) {
    D <- as.dist(1-abs(S))
  }
  else
  {
    p <- nrow(S)
    D <- as.dist(1)
    D[p*(p-1)/2]=1
    attr(D, "Size") <- p
    attr(D, "Labels") <- colnames(S)
    nonzero <- which(S!=0)
    nonzero.row <- (nonzero-1) %% p + 1
    nonzero.col <- (nonzero - 1) %/% p + 1
    if(progress) pb <- txtProgressBar(min = 0, max = length(nonzero), style = 3)
    for(k in 1:length(nonzero)) {
      position <- nonzero[k]
      i <- nonzero.row[k]
      j <- nonzero.col[k]
      position.new <- position - j*(j+1)/2
      D[position.new] <- 1-abs(S[i,j])
      if(progress) {
          if(k %% 1000 == 0 | k==length(nonzero)) setTxtProgressBar(pb, k)
      }
    }
    D[is.na(D)] <- 1
  }
  if(progress) cat("\n")
  return(D)
}

cut.dendrogram <- function(dendrogram, h=NULL, k=NULL) {
    # Verify that either h or r parameter was provided
    stopifnot(!(is.null(h) & is.null(r)))
    
    if(is.null(h)) {
        k.min <- -0.1
        k.max <- 1.1
        while( abs(h.max-h.min) > 1e-12) {
            k.center <- (k.min+k.max)/2
            groups <- cutree(dendrogram, k=k.center)
            n.groups <- max(groups)
            if( n.groups <= k ) {
                k.max <- k.center
            }
            if( n.groups >= k ) {
                k.min <- k.center
            }
        }        
    }
    else {
        groups <- cutree(dendrogram, h=h)
    }
    return(groups)
}

is.dendrogram.nested <- function(dendrogram, resolution.list=c(1,2,5,10,20,50)) {
    # Choose groups by cutting the dendrogram
    resolution.list <- resolution.list/100
    groups.list <- sapply(resolution.list, function(resolution) {
        cutree(Sigma.clust, k = round(resolution*nrow(Variants)))
    })
    colnames(groups.list) <- resolution.list

    all.contained <- sapply(seq(ncol(groups.list)-1), function(j) {
        is.contained <- sapply(unique(groups.list[,j+1]), function(val) {
            val.idx <- which(groups.list[,j+1]==val)
            contained <- length(unique(groups.list[val.idx,j]!=val))==1
        })
        sum(!is.contained)==0
    })

    return(sum(!all.contained)==0)
}
