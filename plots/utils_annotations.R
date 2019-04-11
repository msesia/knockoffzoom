library(ggrepel)
library(latex2exp)
library(egg)
library(grid)

font.size <- 20

place_segments <- function(Segments, gap=5e4, verbose=FALSE) {

  # Define function that checks whether a new segment would fit in an existing row
  segment_fits <- function(start, end, row, gap) {
    if(length(row)==0) {
      return(TRUE)
    }
    for(segment in row) {
      if(start<=segment[2]+gap && end>=segment[1]-gap){
        return(FALSE)
      }
    }
    return(TRUE)
  }

  # Count number of segments
  n.segments <- nrow(Segments)
  segment.heights <- rep(NA, n.segments)
  segment.rows <- list()

  n.rows <- 0

  for(j in 1:n.segments) {
    start <- Segments$start[j]
    end <- Segments$end[j]

    if(length(segment.rows)==0) {
      # Add the first segment to the first row
      segment.rows[[1]] <- list()
      segment.rows[[1]][[1]] <- c(start,end)
      segment.heights[j] <- 1
      n.rows <- n.rows + 1
    } else {
      would.fit <- sapply(1:length(segment.rows), function(g) {
        segment_fits(start, end, segment.rows[[g]], gap=gap)
      })
      if(any(would.fit)) {
        # Add segment to existing row
        k <- min(which(would.fit))
        #print(length(segment.rows[[k]]))
        row.length <- length(segment.rows[[k]])
        if(verbose) cat(sprintf("Appending (%d,%d) to row %d \n",start,end,k))
        segment.rows[[k]][[row.length+1]] <- c(start,end)
        segment.heights[j] <- k
      } else {
        # Add segment to new row
        n.rows <- n.rows + 1
        if(verbose) cat(sprintf("No room to append (%d,%d). Creating new row %d. \n",start,end,n.rows))
        segment.rows[[n.rows]] <- list()
        segment.rows[[n.rows]][[1]] <- c(start,end)
        segment.heights[j] <- n.rows
      }
    }
  }

  # Assign heights to each segment
  Segments$Height <- segment.heights

  # Return modified data
  return(Segments)
}

plot_sears_tower <- function(window.chr, window.left, window.right, Discoveries, LMM, LMM.clumped,
                             Annotations.func, Exons.canonical, highlight.gene=NULL, max.gene.rows=4) {
    # Significance threshold
    p.significant <- 5e-8

    # Extract color map
    annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>%
        ungroup() %>%
        mutate(name.num=parse_number(as.character(name))) %>%
        mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
        mutate(label=as.factor(label)) %>%
        arrange(name.num)

    # Convert names to factors according to color maps
    Annotations.func <- Annotations.func %>%
        mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
        mutate(label=factor(label, levels=annotation.color.map$label, labels=annotation.color.map$label))

    # Extract knockoff discoveries within this window
    Knockoffs.window <- Discoveries %>% filter(Method=="Knockoffs") %>%
        filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))

    # # Update window limits
    # if(nrow(Knockoffs.window)>0) {
    #   window.left <- min(window.left, min(Knockoffs.window$BP.min))
    #   window.right <- max(window.right, max(Knockoffs.window$BP.max))
    # }

    # Extract LMM pvalues within this window
    LMM.clumped.window <- LMM.clumped %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left)
    LMM.window <- LMM %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left) %>%
        left_join(LMM.clumped.window, by = c("SNP", "CHR", "BP")) %>%
        mutate(BP.lead=factor(BP.lead))
    
    cat(sprintf("There are %d LMM pvalues within this window, %d of which are significant.\n",
                nrow(LMM.window), sum(LMM.window$P<p.significant)))

    # Select exons within this windows
    Exons.window <- Exons.canonical %>%
        filter(chrom==window.chr, txStart<=window.right, txEnd>=window.left)
    cat(sprintf("There are %d exons within this window, divided into %d genes.\n",
                nrow(Exons.window), length(unique(Exons.window$name2))))

    # Select functional annotations within this window
    Functional.window <- Annotations.func %>%
        filter(chrom==window.chr, chromStart<=window.right, chromEnd>=window.left)
    cat(sprintf("There are %d functional annotations within this window.\n",
                nrow(Functional.window)))
    
    # Significance level for pvalues
    Window.nominal <- LMM.window %>% mutate(Importance = -log10(p.significant))

    # Minimal theme
    theme_minimal <- theme_bw() +
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.border=element_blank(),
              panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray")
             )

    # Manhattan plot
    LMM.window$BP.lead <- round((LMM.window$BP.lead %>% as.character %>% parse_number)/1e6,3) %>% factor
    clump.lead.snps <- unique(LMM.clumped.window$SNP.lead)

    if(all(is.na(LMM.window$BP.lead))) {
        p.manhattan <- LMM.window %>%
            #filter(P<1e-2) %>%
            mutate(P=pmax(1e-300,P)) %>%
            ggplot(aes(x=BP, y=-log10(P))) +
            geom_point(color="black", alpha=0.25)
    } else {
        p.manhattan <- LMM.window %>%
            #filter(P<1e-2) %>%
            mutate(SNP.lead=factor(SNP.lead, levels=clump.lead.snps, labels=clump.lead.snps)) %>%
            mutate(P=pmax(1e-300,P)) %>%
            ggplot(aes(x=BP, y=-log10(P), color=SNP.lead, alpha=is.na(SNP.lead))) +
            geom_point() +
            scale_colour_discrete(na.value = "black", name="Clump lead SNP", guide=FALSE)
    }
    if(min(LMM.window$P)<1e-100) {
        manhattan.y.max <- 300
        manhattan.y.breaks <- c(0,7.3,100,300)
    } else {
        manhattan.y.max <- 100
        manhattan.y.breaks <- c(0,7.3,50,100)
    }

    p.manhattan <- p.manhattan +
        geom_hline(data=Window.nominal, aes(yintercept=Importance), linetype="dashed", color = "black") +
        scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 1), guide=FALSE) +
        scale_x_continuous(expand=c(0.01,0.01), limits=c(window.left,window.right),
                           labels=function(x) {x*1e-6}) +
        scale_y_continuous(limits=c(1.8,manhattan.y.max),
                           breaks=manhattan.y.breaks, trans="sqrt") +
        theme_bw() +
        xlab(sprintf("Chromosome %d (Mb)", window.chr)) + ylab(TeX("$-\\log_{10}(p)$")) +
        theme(text = element_text(size=font.size),
              legend.text = element_text(size=font.size),
              panel.border = element_blank(),
              panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray")
              ) +
        ggtitle(sprintf("Manhattan plot (LMM)"))

    # Plot clumped LMM discoveries
    if(nrow(LMM.clumped.window)>0) {
        LMM.clumps.window <- LMM.clumped.window %>% group_by(CHR, SNP.lead, BP.lead) %>%
            summarise(BP.min=min(BP), BP.max=max(BP)) %>%
            ungroup() %>%
            mutate(BP.lead=factor(BP.lead)) %>%
            mutate(SNP.lead=factor(SNP.lead, levels=clump.lead.snps, labels=clump.lead.snps)) %>%
            mutate(start=BP.min, end=BP.max) %>% mutate(width=end-start) %>% arrange(desc(width)) %>%
            place_segments(gap=1e4) %>%
            mutate(Height=-Height)

        p.clumped <- LMM.clumps.window %>%
            ggplot() +
            geom_segment(aes(x=pmax(BP.min,window.left), y=Height,
                             xend=pmin(BP.max,window.right), yend=Height, color=SNP.lead)) +
            geom_segment(aes(x=BP.min, y=Height-0.4, xend=BP.min, yend=Height+0.4, color=SNP.lead)) +
            geom_segment(aes(x=BP.max, y=Height-0.4, xend=BP.max, yend=Height+0.4, color=SNP.lead)) +
            scale_colour_discrete(na.value = "black", name="Clump lead SNP") +
            scale_x_continuous(expand=c(0.01,0.01), limits=c(window.left,window.right),
                               labels=function(x) {x*1e-6}) +
            scale_y_continuous(breaks = NULL) +
            ylab("") + xlab("") +
            theme_void() +
            theme(text = element_text(size=font.size),
                  legend.text = element_text(size=font.size),
                  legend.key.size = unit(0.5,"line"))

      # Determine whether we would list all lead SNPs in the legend
      if(length(unique(LMM.clumps.window$SNP.lead))>20) {
        p.clumped <- p.clumped +
          theme(legend.position = "none")
      }
    } else {
      LMM.clumps.window <- tibble()
      p.clumped <- ggplot(tibble()) + geom_blank()
    }

    # Plot knockoff discoveries
    resolution.list <- c("2","5","10","20","50","100") %>% rev
    resolution.labels <- resolution.list
    resolution.heights <- seq(length(resolution.list))
    names(resolution.heights) <- resolution.list

    if(nrow(Knockoffs.window)>0) {
      p.knockoffs <- Knockoffs.window %>%
          mutate(Resolution=as.character(Resolution)) %>%
          mutate(Height=resolution.heights[Resolution]) %>%
          ggplot() +
          geom_rect(aes(xmin=pmax(BP.min,window.left), xmax=pmin(BP.max,window.right), ymin=Height-0.5, ymax=Height+0.5),
                    alpha=0.5, fill="black", color="black") +
          ylab("Resolution (%)") + xlab("") +
          scale_x_continuous(expand=c(0.01,0.01), limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
          scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                             labels=resolution.labels, breaks=resolution.heights) +
          ggtitle("Chicago plot (KnockoffZoom)") +
          theme_bw() +
          theme(panel.grid.minor.y = element_blank(),
                axis.line=element_blank(), axis.text.x=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                panel.border=element_blank(),
                panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
                panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
                text = element_text(size=font.size)
                )
    } else {
      p.knockoffs <- ggplot(tibble()) + geom_blank()
    }
    
    # Plot functional annotations
    myColors <- annotation.color.map$itemColor
    names(myColors) <- annotation.color.map$label

    p.functional <- Functional.window %>%
        mutate(chromStart=pmax(chromStart, window.left), chromEnd=pmin(chromEnd, window.right)) %>%
        ggplot() +
        geom_rect(aes(xmin=chromStart, xmax=chromEnd, ymin=0.5, ymax=1.5, fill=label)) +
        ylab("") + xlab("") +
        scale_x_continuous(expand=c(0.01,0.01), limits=c(window.left,window.right),
                           labels=function(x) {x*1e-6}) +
        scale_color_manual(values=myColors, guide=FALSE) +
        scale_fill_manual(values=myColors, name="Annotation") +
        ggtitle("Functional annotations") +
        theme_void() +
        theme(legend.key.size = unit(0.5,"line"),
              text = element_text(size=font.size),
              legend.text = element_text(size=font.size))

    # Find out how many genes there are and determine whether we would plot all of them
    Genes.window <- Exons.window %>% group_by(name, name2, strand) %>%
        summarise(txStart=min(txStart), txEnd=max(txEnd)) %>%
        mutate(txStart=max(txStart, window.left), txEnd=min(txEnd, window.right)) %>%
        mutate(start=txStart, end=txEnd)

    # Highlight special gene, if available
    if(is.null(highlight.gene)) {
        Genes.window$highlight <- FALSE
    } else {
        Genes.window <- Genes.window %>% mutate(highlight = (name2==highlight.gene))
    }

    # Do not attempt to place more than a max number of genes
    max.genes <- 50
    if(length(unique(Genes.window$name))<=max.genes) {
      plot.all.genes <- TRUE
    } else {
      plot.all.genes <- FALSE
      genes.toshow <- unique(Genes.window$name)[1:max.genes]
      Genes.window <- Genes.window %>% filter((name %in% genes.toshow) || (highlight))
    }

    # Sort the genes by length, making sure that the highlighted gene is first;
    # then, place them on different rows so they don't overlap
    Genes.window <- Genes.window %>% mutate(width = end-start) %>% arrange(desc(highlight), desc(width)) %>%
        place_segments()
    
    # Remove genes that do not fit in 3 rows
    if(length(unique(Genes.window$Height))>max.gene.rows) {
        plot.all.genes <- FALSE
        Genes.window <- Genes.window %>% filter(Height<=max.gene.rows)
    }

    # Rescale the gene heights
    Genes.window <- Genes.window %>% mutate(Height=-(Height-0.5)) %>%
        mutate(txCenter=(txStart+txEnd)/2, Height=Height/2)
            
    # Extract genes within this window
    strand.levels <- c("+", "-")
    strand.labels <- c("->", "<-")
    Genes.window <- Genes.window %>%
        mutate(strand.latex=factor(strand, labels=strand.labels, levels=strand.levels)) %>%
        mutate(strand.latex=as.character(strand.latex))

    # Plot exons and genes
    p.genes <- Exons.window %>%
        filter(exonStarts>=window.left, exonEnds<=window.right) %>%
        inner_join(Genes.window %>% select(name, name2, Height), by = c("name", "name2")) %>%
        ggplot() +
        geom_rect(aes(xmin=exonStarts, xmax=exonEnds, ymin=Height-0.1, ymax=Height+0.1),
                  alpha=1, color="black", fill="black") +
        geom_segment(data=Genes.window, aes(x=txStart, y=Height, xend=txEnd, yend=Height, group=name2),
                     color="black") +
        geom_label_repel(data=Genes.window, aes(x=txCenter, y=Height,
                                                label=paste(name2,strand.latex,sep=" "),
                                                fill=highlight),
                         size=5, direction="both", force=1, max.iter=10000,
                         box.padding=0.1,
                         point.padding=0.1,
                         label.padding=0.15, 
                         segment.color = 'grey50', segment.alpha=0.5, seed=2019) +
        ylab("") + xlab("") +
        scale_x_continuous(expand=c(0.01,0.01), limits=c(window.left,window.right),
                           labels=function(x) {x*1e-6}) +
        scale_y_continuous(expand=c(0.01,0.01), breaks = NULL,
                           limits=c(min(Genes.window$Height)-0.25,max(Genes.window$Height)+0.25)) +
        scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "yellow"), guide=FALSE) +
        theme_void() +
        ggtitle("Genes") +
        theme(text = element_text(size=font.size),
              panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray")
              )
    p.genes

    if(!plot.all.genes) {
      p.genes <- p.genes + ggtitle("Genes (not all genes are shown)")
    }

    # Determine relative heights of each subplot
    height.manhattan <- 0.75
    height.clumps <- 0.3
    height.knockoffs <- 1.25
    height.functional <- 0.3
    height.genes <- 0.75*length(unique(Genes.window$Height))
    heights <- c(height.manhattan, height.clumps, height.knockoffs, height.functional, height.genes)

    debug.lines <- FALSE
    g1 <-  ggplotGrob(p.manhattan)
    g2 <-  ggplotGrob(p.clumped + theme(legend.position="none"))
    g3 <-  ggplotGrob(p.knockoffs)
    g4 <-  ggplotGrob(p.functional + theme(legend.position="none"))
    g5 <-  ggplotGrob(p.genes)
    fg1 <- gtable_frame(g1, width = unit(1, "null"), height = unit(heights[1], "null"), debug = debug.lines)
    fg2 <- gtable_frame(g2, width = unit(1, "null"), height = unit(heights[2], "null"), debug = debug.lines)
    fg3 <- gtable_frame(g3, width = unit(1, "null"), height = unit(heights[3], "null"), debug = debug.lines)
    fg4 <- gtable_frame(g4, width = unit(1, "null"), height = unit(heights[4], "null"), debug = debug.lines)
    fg5 <- gtable_frame(g5, width = unit(1, "null"), height = unit(heights[5], "null"), debug = debug.lines)

    fg.l <- gtable_frame(gtable_rbind(fg1, fg2, fg3, fg4, fg5),
                         width = unit(4, "null"),
                         height = unit(1, "null"))

    # Extract legends
    g6 <- ggplotGrob(ggplot())
    try(g6 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.clumped))+
                         theme(text = element_text(size=font.size))
                         ), silent=TRUE)
    g7 <- ggplotGrob(ggplot())
    try(g7 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.functional))+
                         theme(text = element_text(size=font.size))
                         ), silent=TRUE)
    
    g0 <- ggplotGrob(ggplot())
    fg00 <- gtable_frame(g0, width = unit(1, "null"),
                        height = unit(0.1, "null"), debug = debug.lines)
    fg6 <- gtable_frame(g6, width = unit(1, "null"),
                        height = unit(1, "null"), debug = debug.lines)
    fg01 <- gtable_frame(g0, width = unit(1, "null"),
                        height = unit(0.15, "null"), debug = debug.lines)
    fg7 <- gtable_frame(g7, width = unit(1, "null"),
                        height = unit(1, "null"), debug = debug.lines)
    fg02 <- gtable_frame(g0, width = unit(1, "null"),
                        height = unit(1, "null"), debug = debug.lines)

    fg.r <- gtable_frame(gtable_rbind(fg00,fg6,fg01,fg7,fg02),
                         width = unit(1, "null"),
                         height = unit(1, "null"))
   
    grid.newpage()
    combined <- gtable_frame(gtable_cbind(fg.l, fg.r),
                         width = unit(1, "null"),
                         height = unit(1, "null"))
    p.final <- ggplotify::as.ggplot(combined)
    
    return(p.final)
}
