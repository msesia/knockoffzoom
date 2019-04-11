plot_sears_tower <- function(window.chr, window.left, window.right, Selections, LMM, Clumped, plot.file=NULL) {

    p.significant <- 5e-8

    # Color palette
    causal.cols <- c("TRUE"="red", "FALSE"="black", "NA"="black")
    causal.alphas <- c("TRUE"=1, "FALSE"=0.25, "NA"=0.25)

    # Extract knockoff discoveries within this window
    Knockoffs.window <- Selections %>% filter(Method=="Knockoffs") %>%
        filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))

    # Update window limits
    window.left <- min(window.left, min(Knockoffs.window$BP.min))
    window.right <- max(window.right, max(Knockoffs.window$BP.max))

    # Extract LMM pvalues within this window
    LMM.clumped.window <- Clumped %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left) %>%
        select(-Causal)
    LMM.window <- LMM %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left) %>%
        left_join(LMM.clumped.window, by = c("SNP", "CHR", "BP")) %>%
        mutate(BP.lead=factor(BP.lead))

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
              panel.grid.minor.x = element_line(size = 0.1, colour = "lightgray"),
             )

    # Manhattan plot
    LMM.window$BP.lead <- round((LMM.window$BP.lead %>% as.character %>% parse_number)/1e6,3) %>% factor

    if(all(is.na(LMM.window$BP.lead))) {
        p.manhattan <- LMM.window %>%
            mutate(Causal=as.character(Causal)) %>%
            mutate(P=pmax(1e-300,P)) %>%
            ggplot(aes(x=BP, y=-log10(P), color=Causal, alpha=Causal)) +
            geom_point()
    } else {
        p.manhattan <- LMM.window %>%
            mutate(Causal=as.character(Causal)) %>%
            mutate(P=pmax(1e-300,P)) %>%
            ggplot(aes(x=BP, y=-log10(P), color=Causal, alpha=Causal)) +
            geom_point()
    }
    p.manhattan <- p.manhattan +
        geom_hline(data=Window.nominal, aes(yintercept=Importance), linetype="dashed", color = "black") +
        scale_alpha_manual(values = causal.alphas, guide=FALSE) +
        scale_color_manual(values = causal.cols) +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(limits=c(0,300), breaks=c(0,7.3,50,100,200,300), trans="sqrt") +        
        theme_bw() +
        xlab(sprintf("Chromosome %d (Mb)", window.chr)) + ylab("$-\\log_{10}(p)$") +
        theme(panel.grid.minor.y = element_blank(),
              panel.border = element_blank()) +
        theme(legend.justification = c(-0.0,1.0), legend.position = c(0,1),
          legend.title=element_text(size=9), legend.text=element_text(size=9),
          legend.key.width = unit(0.85, "cm"),
          legend.background = element_rect(fill=alpha('white', 0.0)),
          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor.x = element_line(size = 0.1, colour = "lightgray"),
          )
    p.manhattan

    # Plot clumped LMM discoveries
    LMM.clumps.window <- LMM.clumped.window %>% group_by(CHR, SNP.lead, BP.lead) %>%
        summarise(BP.min=min(BP), BP.max=max(BP)) %>%
        ungroup() %>%
        mutate(BP.lead=factor(BP.lead))

    clump.list <- levels(LMM.clumps.window$BP.lead)
    clump.heights <- seq(length(clump.list)) %>% rev
    names(clump.heights) <- clump.list

    p.clumped <- LMM.clumps.window %>%
        mutate(Height=clump.heights[BP.lead]) %>%
        ggplot() +
        geom_segment(aes(x=pmax(BP.min,window.left), y=Height, xend=pmin(BP.max,window.right),
                         yend=Height, color=BP.lead)) +
        geom_segment(aes(x=BP.min, y=Height-0.4, xend=BP.min, yend=Height+0.4, color=BP.lead)) +
        geom_segment(aes(x=BP.max, y=Height-0.4, xend=BP.max, yend=Height+0.4, color=BP.lead)) +
        scale_colour_discrete(na.value = "black", guide=FALSE) +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(breaks = NULL) +
        ylab("") + xlab("") +
        theme_void()

    # Plot knockoff discoveries
    resolution.list <- levels(Knockoffs.window$Resolution)
    resolution.heights <- seq(length(resolution.list))
    names(resolution.heights) <- resolution.list
    #resolution.labels <- paste(parse_number(resolution.list), "\\%", sep="")
    resolution.labels <- c("0.226", "0.088", "0.042", "0.018", "0.004", "0.001", "0.000")

    p.knockoffs <- Knockoffs.window %>%
        mutate(Causal=as.character(Causal)) %>%
        mutate(Height=resolution.heights[Resolution]) %>%
        ggplot() +
        geom_rect(aes(xmin=pmax(BP.min,window.left), xmax=pmin(BP.max,window.right),
                      ymin=Height-0.5, ymax=Height+0.5, fill=Causal, color=Causal),
                  alpha=0.25) +
        ylab("Resolution") + xlab("") +
        scale_fill_manual(values = causal.cols, guide=FALSE) +
        scale_color_manual(values = causal.cols, guide=FALSE) +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                           labels=resolution.labels, breaks=resolution.heights) +
        theme_bw() +
        theme(panel.grid.minor.y = element_blank(),
              axis.line=element_blank(), axis.text.x=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              panel.border=element_blank(),
              panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
              )
    p.knockoffs

    height.manhattan <- 1.25
    height.clumps <- 0.2
    height.knockoffs <- 1.25

    heights <- c(height.manhattan, height.knockoffs)

    p.towers <- egg::ggarrange(p.manhattan, p.knockoffs, ncol=1, draw=F, heights=heights)

    if(!is.null(plot.file)) {
        ggsave(plot.file, p.towers, width=7, height=sum(heights))
    }

    return(p.towers)
}
