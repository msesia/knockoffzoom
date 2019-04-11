load_annotations = function(data_dir){
  # Functional annotations
  annotations.func.file <- sprintf("%s/annotations/GM12878_ChromHMM_hg19.txt", data_dir)
  Annotations.func.raw <- read_tsv(annotations.func.file, col_types=cols(`#bin`=col_integer(), 
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
  
  # Gene annotations
  annotations.genes.file <- "../data/annotations/NCBI_genes_hg19.txt"
  annotations.genes.file <- sprintf("%s/annotations/NCBI_genes_hg19.txt", data_dir)
  Annotations.genes.raw <- read_tsv(annotations.genes.file, col_types=cols(`#bin`=col_integer(), 
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
  
  #colnames(Annotations.genes)
  #Annotations.genes %>% head
  
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
  
  annotations = c()
  annotations$Annotations.func = Annotations.func
  annotations$Exons.canonical = Exons.canonical
  
  return(annotations)
}

load_association_results = function(data_dir, phenotype){
  # Knockoffs
  resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
  phenotype.list <- c(phenotype)
  
  Params <- expand.grid(Resolution=resolution.list, Phenotype=phenotype.list) %>% as_tibble()
  
  Discoveries <- lapply(1:nrow(Params), function(idx) {
    resolution <- Params$Resolution[idx]
    phenotype <- Params$Phenotype[idx]
    
    # Load LMM discoveries (if available)
    lmm.file <- sprintf("%s/discoveries/%s_lmm_regions.txt", data_dir, phenotype)
    if(resolution=="Radj2" && file.exists(lmm.file)) {
      Discoveries.LMM <- read_delim(lmm.file, delim=" ", col_types=cols()) %>%
        mutate(Phenotype=phenotype, Method="LMM", Importance=-log10(P), Resolution="GWAS") %>%
        select(-c("P"))
    } else {
      Discoveries.LMM <- tibble()
    }
    
    # Load knockoffs discoveries (if available)
    knockoffs.file <- sprintf("%s/discoveries/%s_knockoffs_%s.txt", data_dir, phenotype, resolution)
    if(file.exists(knockoffs.file)) {
      Discoveries.knockoffs <- read_delim(knockoffs.file, delim=" ", col_types=cols()) %>%
        mutate(Phenotype=phenotype, Method="Knockoffs", Importance=W, Resolution=resolution) %>%
        select(-c("W", "Group"))
    } else {
      Discoveries.knockoffs <- tibble()
    }
    
    # Combine results
    return(rbind(Discoveries.LMM,Discoveries.knockoffs))
  })
  Discoveries <- do.call("rbind", Discoveries) %>%
    mutate(Resolution=factor(Resolution, labels=c("GWAS", resolution.list), levels=c("GWAS", resolution.list)))
  
  Discoveries %>% mutate(BP.width=BP.max-BP.min) %>%
    group_by(Phenotype, Resolution) %>%
    summarise(Discoveries=n(), Size=median(Size), BP.width=median(BP.width)/1e6)
  
  # Load LMM p-values
  lmm.file <- sprintf("%s/discoveries/%s_stats.txt", data_dir, phenotype)
  LMM <- read_delim(lmm.file, delim="\t", col_types=cols()) %>% as_tibble()
  if("P_BOLT_LMM" %in% colnames(LMM)) {
    LMM <- LMM %>% mutate(P=P_BOLT_LMM)
  } else {
    LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
  }
  
  # Load clumped LMM results
  lmm.file <- sprintf("%s/discoveries/%s_lmm_variants.txt", data_dir, phenotype)
  LMM.clumped <- read_delim(lmm.file, delim=" ", col_types=cols()) %>%
    mutate(Phenotype=phenotype, Method="LMM", Importance=-log10(P), Resolution="GWAS") %>%
    select(-c("P"))
  
  # Load knockoff stats
  stats.file <- sprintf("../data/discoveries/%s_Radj2_lasso_stats.txt", phenotype)
  Stats <- read_delim(stats.file, delim=" ", col_types=cols())
  Stats <- Stats %>% filter(is.na(R2)|R2<=0.99)
  
  association_results = c()
  association_results$Discoveries = Discoveries
  association_results$Stats = Stats
  association_results$LMM = LMM
  association_results$LMM.clumped = LMM.clumped
  
  return(association_results)
}

plot_sears_tower_state = function(state, annotations){
  plot_sears_tower(state$chr, state$window.left, state$window.right,
                   state$association_results$Discoveries,
                   state$association_results$LMM,
                   state$association_results$LMM.clumped,
                   annotations$Annotations.func,
                   annotations$Exons.canonical,
                   highlight.gene=state$highlight.gene)
}

adjust_window = function(state){
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