load_annotations <- function(data_dir) {
    
    # Functional annotations
    annotations.func.file <- sprintf("%s/annotations/GM12878_ChromHMM_hg19.txt", data_dir)
    if(file.exists(annotations.func.file)) {
        Annotations.func.raw <- read_tsv(annotations.func.file,
                                         col_types=cols(`#bin`=col_integer(), 
                                                        chromStart=col_integer(), 
                                                        chromEnd=col_integer(),
                                                        name=col_character(),
                                                        score=col_integer(), 
                                                        thickStart=col_integer(), 
                                                        thickEnd=col_integer()))
        Annotations.func.raw <- Annotations.func.raw %>% mutate(name=as.factor(name))
        
        # Remove weird chromosomes and convert colors
        valid.chrom <- paste("chr", seq(1,22), sep="")
        Annotations.func <- Annotations.func.raw %>%
            filter(chrom %in% valid.chrom) %>% mutate(chrom=parse_number(as.character(chrom))) %>%
            separate(itemRgb, into=c("itemR", "itemG", "itemB"), sep=",", convert=T) %>%
            mutate(itemColor = rgb(red=itemR, blue=itemB, green=itemG, maxColorValue=255)) %>%
            arrange(chrom, chromStart)
        
        # Extract color map
        annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>% 
            ungroup() %>%
            mutate(name.num=parse_number(as.character(name))) %>% 
            mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
            arrange(name.num)
        
        # Convert names to factors according to color maps
        Annotations.func <- Annotations.func %>%
            mutate(name=factor(name, levels=annotation.color.map$name, labels=annotation.color.map$name))
    } else {
        Annotations.func <- NULL
    }
        
    # Gene annotations
    annotations.genes.file <- sprintf("%s/annotations/NCBI_genes_hg19.txt", data_dir)
    if(file.exists(annotations.genes.file)) {
        Annotations.genes.raw <- read_tsv(annotations.genes.file,
                                          col_types=cols(`#bin`=col_integer(), 
                                                         txStart=col_integer(), 
                                                         txEnd=col_integer(), 
                                                         cdsStart=col_integer(), 
                                                         cdsEnd=col_integer(), 
                                                         exonCount=col_integer(), 
                                                         score=col_integer(), 
                                                         .default=col_character()))
        
        # Remove weird chromosomes
        valid.chrom <- paste("chr", seq(1,22), sep="")
        Annotations.genes <- Annotations.genes.raw %>%
            filter(chrom %in% valid.chrom) %>% mutate(chrom=parse_number(as.character(chrom))) %>%
            arrange(chrom)
        
        # Split rows corresponding to same gene but different exons
        Exons <- Annotations.genes %>% 
            separate_rows(exonStarts, exonEnds, exonFrames, sep=",", convert=TRUE) %>%
            drop_na()
        
        # Pick the canonical transcripts
        # That is, for each unique "name2", keep only the rows corresponding to the "name" 
        # with the largest sum of exon lengths
        Exons.canonical <- Exons %>%
            mutate(exonLength=exonEnds-exonStarts) %>%
            group_by(name, name2) %>% summarise(Length=sum(exonLength)) %>%
            ungroup() %>% group_by(name2) %>% top_n(1, Length) %>%
            inner_join(Exons, by=c("name", "name2")) 
    } else {
        Exons.canonical <- NULL
    }

    # Return annotations
    annotations <- c()
    annotations$Annotations.func <- Annotations.func
    annotations$Exons.canonical <- Exons.canonical    
    return(annotations)
}

load_association_results <- function(data.dir, phenotype){
    # Knockoffs
    resolution.list <- c("2", "5", "10", "20", "50", "100")
    phenotype.list <- c(phenotype)

    # Load KnockoffZoom discoveries (if available)
    Discoveries <- lapply(resolution.list, function(resolution) {
        knockoffs.file <- sprintf("%s/%s_res%s_discoveries.txt", data.dir, phenotype, resolution)
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
    lmm.file <- sprintf("%s/%s_lmm.txt", data.dir, phenotype)
    if(file.exists(lmm.file)) {
        LMM <- read_delim(lmm.file, delim="\t", col_types=cols()) %>% as_tibble()
        if("P_BOLT_LMM" %in% colnames(LMM)) {
            LMM <- LMM %>% mutate(P=P_BOLT_LMM)
        } else {
            LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
        }
    } else {
        LMM <- NULL
    }

    # Load clumped LMM results
    lmm.file <- sprintf("%s/%s_lmm_clumped.txt", data.dir, phenotype)
    if(file.exists(lmm.file)) {
        LMM.clumped <- read_delim(lmm.file, delim=" ", col_types=cols()) %>%
            mutate(Phenotype=phenotype, Method="LMM", Importance=-log10(P), Resolution="GWAS") %>%
            select(-c("P"))
    } else {
        LMM.clumped <- NULL
    }

    # Load knockoff stats
    stats.file <- sprintf("%s/%s_res%s_stats.txt", data.dir, phenotype, "2")
    if(file.exists(stats.file)) {
        Stats <- read_delim(stats.file, delim=" ", col_types=cols())
    } else {
        Stats <- NULL
    }

    association_results = c()
    association_results$Discoveries = Discoveries
    association_results$Stats = Stats
    association_results$LMM = LMM
    association_results$LMM.clumped = LMM.clumped

    return(association_results)
}

plot_combined_state <- function(state, annotations){
    
    if(is.null(annotations$Annotations.func)) {
        showNotification("Missing annotations file.")
    }
    if(is.null(annotations$Exons.canonical)) {
        showNotification("Missing genes file.")
    }
        
    plot_combined(state$chr, state$window.left, state$window.right,
                  state$association_results$Discoveries,
                  state$association_results$LMM,
                  state$association_results$LMM.clumped,
                  Annotations.func=annotations$Annotations.func,
                  Exons.canonical=annotations$Exons.canonical,
                  highlight.gene=state$highlight.gene)
}

adjust_window <- function(state){
    # Extract knockoff discoveries within this window
    Knockoffs.window <- state$association_results$Discoveries %>% filter(Method=="Knockoffs") %>%
        filter(CHR==state$chr, BP.min<=state$window.right, BP.max>=state$window.left)

    # Update window limits
    if(nrow(Knockoffs.window)>0) {
        state$window.left <- min(state$window.left, min(Knockoffs.window$BP.min))
        state$window.right <- max(state$window.right, max(Knockoffs.window$BP.max))
    }
  
  return(state)
}
