.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
source("utils_clumping.R")
source("utils_shiny.R")
source("utils_manhattan.R")
source("utils_annotations.R")

data_dir = "../results"
chromosomes = 1:22
phenotype = "example"

load_association_results = function(data_dir, phenotype){
  # Knockoffs
  resolution.list <- c("2", "5", "10", "20", "50", "100")
  phenotype.list <- c(phenotype)

  Params <- expand.grid(Resolution=resolution.list, Phenotype=phenotype.list) %>% as_tibble()

  Discoveries <- lapply(1:nrow(Params), function(idx) {
    resolution <- Params$Resolution[idx]
    phenotype <- Params$Phenotype[idx]
    #cat(sprintf("Loading list of discoveries for %s at resolution %s...\n", phenotype, resolution))

    # Load knockoffs discoveries (if available)
    knockoffs.file <- sprintf("%s/%s_res%s.txt", data_dir, phenotype, resolution)
    if(file.exists(knockoffs.file)) {
      Discoveries.knockoffs <- read_delim(knockoffs.file, delim=" ", col_types=cols()) %>%
        mutate(Phenotype=phenotype, Method="Knockoffs", Importance=W, Resolution=resolution) %>%
        select(-c("W", "Group"))
      #cat(sprintf("Found %d discoveries for %s made with knockoffs at resolution %s.\n",
      #            nrow(Discoveries.knockoffs), phenotype, resolution))
    } else {
      #cat(sprintf("Discoveries for %s made with knockoffs at resolution %s are not available.\n",
      #    phenotype, resolution))
      Discoveries.knockoffs <- tibble()
    }

    # Combine results
    return(Discoveries.knockoffs)
  })
  Discoveries <- do.call("rbind", Discoveries) %>%
    mutate(Resolution=factor(Resolution, labels=c("GWAS", resolution.list), levels=c("GWAS", resolution.list)))

  # Load LMM p-values
  lmm.file <- sprintf("%s/%s_lmm.txt", data_dir, phenotype)
  LMM <- read_delim(lmm.file, delim="\t", col_types=cols()) %>% as_tibble()
  if("P_BOLT_LMM" %in% colnames(LMM)) {
    LMM <- LMM %>% mutate(P=P_BOLT_LMM)
  } else {
    LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
  }

  association_results = c()
  association_results$Discoveries = Discoveries
  association_results$LMM = LMM

  return(association_results)
}

association_results <- load_association_results(data_dir, phenotype)

#Make annotation plot
source("utils_annotations.R")
window.chr = 22
window.center = mean(filter(association_results$LMM, CHR==window.chr)$BP/1e6)
window.width = 1
window.left <- 1e6*max(0, window.center - window.width)
window.right <- 1e6*min(window.center + window.width,
                    1e-6*max(filter(association_results$LMM, CHR==window.chr)$BP))
plot_sears_tower(window.chr, window.left, window.right,
                 association_results$Discoveries,
                 association_results$LMM)
