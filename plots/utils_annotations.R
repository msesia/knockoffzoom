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

plot_sears_tower <- function(window.chr, window.left, window.right, Discoveries, LMM) {
    # Significance threshold
    p.significant <- 5e-8

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
    LMM.window <- LMM %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left)

    cat(sprintf("There are %d LMM pvalues within this window, %d of which are significant.\n",
                nrow(LMM.window), sum(LMM.window$P<p.significant)))

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

    p.manhattan <- LMM.window %>%
        #filter(P<1e-2) %>%
        mutate(P=pmax(1e-300,P)) %>%
        ggplot(aes(x=BP, y=-log10(P))) +
        geom_point(color="black", alpha=0.25)

    if(min(LMM.window$P)<1e-100) {
        manhattan.y.max <- 300
        manhattan.y.breaks <- c(0,7.3,100,300)
    } else if(min(LMM.window$P)<1e-50) {
        manhattan.y.max <- 100
        manhattan.y.breaks <- c(0,7.3,50,100)
    } else if(min(LMM.window$P)<1e-20) {
        manhattan.y.max <- 50
        manhattan.y.breaks <- c(0,7.3,20,50)
    } else {
        manhattan.y.max <- 20
        manhattan.y.breaks <- c(0,7.3,10,20)
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

    # Determine relative heights of each subplot
    height.manhattan <- 0.75
    height.knockoffs <- 1.25

    debug.lines <- FALSE
    g1 <-  ggplotGrob(p.manhattan)
    g3 <-  ggplotGrob(p.knockoffs)
    fg1 <- gtable_frame(g1, width = unit(1, "null"), height = unit(height.manhattan, "null"),
                        debug = debug.lines)
    fg3 <- gtable_frame(g3, width = unit(1, "null"), height = unit(height.knockoffs, "null"),
                        debug = debug.lines)

    fg.l <- gtable_frame(gtable_rbind(fg1, fg3),
                         width = unit(4, "null"),
                         height = unit(1, "null"))
     
    grid.newpage()
    combined <- gtable_frame(fg.l,
                         width = unit(1, "null"),
                         height = unit(1, "null"))
    p.final <- ggplotify::as.ggplot(combined)
    
    return(p.final)
}
