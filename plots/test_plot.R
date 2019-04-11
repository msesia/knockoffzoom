.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
source("utils_shiny.R")
source("utils_clumping.R")
source("utils_manhattan.R")
source("utils_plotting.R")

data_dir = "../results"
chromosomes = 1:22
phenotype = "example"

annotations.dir <- "../data"
annotations <- load_annotations(annotations.dir)

load_association_results = function(data_dir, phenotype){
    # Knockoffs
    resolution.list <- c("2", "5", "10", "20", "50", "100")
    phenotype.list <- c(phenotype)

    # Load KnockoffZoom discoveries (if available)
    Discoveries <- lapply(resolution.list, function(resolution) {
        knockoffs.file <- sprintf("%s/%s_res%s.txt", data_dir, phenotype, resolution)
        if(file.exists(knockoffs.file)) {
            Discoveries.knockoffs <- read_delim(knockoffs.file, delim=" ", col_types=cols()) %>%
                mutate(Phenotype=phenotype, Method="Knockoffs", Importance=W, Resolution=resolution) %>%
                select(-c("W", "Group"))
        } else {
            Discoveries.knockoffs <- tibble()
        }
        return(Discoveries.knockoffs)
    })
    Discoveries <- do.call("rbind", Discoveries)
    
    # Load LMM p-values
    lmm.file <- sprintf("%s/%s_lmm.txt", data_dir, phenotype)
    LMM <- read_delim(lmm.file, delim="\t", col_types=cols()) %>% as_tibble()
    if("P_BOLT_LMM" %in% colnames(LMM)) {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM)
    } else {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
    }

    # Load clumped LMM results
    lmm.file <- sprintf("%s/%s_lmm_clumped.txt", data_dir, phenotype)
    LMM.clumped <- read_delim(lmm.file, delim=" ", col_types=cols()) %>%
        mutate(Phenotype=phenotype, Method="LMM", Importance=-log10(P), Resolution="GWAS") %>%
        select(-c("P"))
    
    association_results = c()
    association_results$Discoveries = Discoveries
    association_results$LMM = LMM
    association_results$LMM.clumped = LMM.clumped

    return(association_results)
}

association_results <- load_association_results(data_dir, phenotype)

#Make annotation plot
source("utils_plotting.R")
window.chr = 22
window.center = mean(filter(association_results$LMM, CHR==window.chr)$BP/1e6)
window.width = 1
window.left <- 1e6*max(0, window.center - window.width)
window.right <- 1e6*min(window.center + window.width,
                    1e-6*max(filter(association_results$LMM, CHR==window.chr)$BP))
plot_combined(window.chr, window.left, window.right,
              association_results$Discoveries,
              association_results$LMM, association_results$LMM.clumped,
              Annotations.func=annotations$Annotations.func,
              Exons.canonical=annotations$Exons.canonical,
              highlight.gene="KIF1B", max.gene.rows=3)
