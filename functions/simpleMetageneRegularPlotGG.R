simpleMetageneRegularPlotGG <- function(
                                    ## INPUT
                                    METAGENE.DT=NULL,
                                    SAMPLE.NAME=NULL,
                                    PLOT.BY.VALUE="frequency",
                                    
                                    ## SETTINGS
                                    ID.VAR="NAME",
                                    
                                    ## PLOT SETTINGS
                                    PIRNA.SIZE=26,
                                    Y.LIMITS=NULL,
                                    Y.BREAKS=seq(0,100,10),
                                    NUC.COLORS = rev(c("firebrick1","dodgerblue1","goldenrod1","forestgreen")),
                                    ASPECT.RATIO = 1,
                                    
                                    ## PLOT SETTINGS
                                    FAM="Helvetica",
                                    XYT=12,
                                    TCOL="black",
                                    LEGEND.POSITION="bottom"
                                    ){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Plot a stadard metagene plot
  ## 06.30.21; Version 6; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("ggrepel")
    library("ggplot2")})
  
  ## INPUT CHECKING
  if(is.null(METAGENE.DT)) stop("Provide a metagene frequency table !!!")
  if(sum(METAGENE.DT[[2]]) != 100) stop("This is not a frequency table !!!")
  if(is.null(SAMPLE.NAME)) stop("Provide metagene SAMPLE.NAME")
  
  ## PROCESS DATA
  message("Processing ...")
  MGm <- data.table::melt(METAGENE.DT,
                          id.vars = ID.VAR, 
                          value.name = PLOT.BY.VALUE,
                          variable.name = "POSITION")
  MGm[["POSITION"]] <- as.numeric(as.character(MGm[["POSITION"]]))
  MGm[["SAMPLE"]] <- dplyr::nth(data.table::tstrsplit(MGm[[ID.VAR]],split="__"),2)
  MGm[["NUCLEOTIDE"]] <- dplyr::nth(data.table::tstrsplit(MGm[[ID.VAR]],split="__"),-1)
  
  ## METRICS & SETTINGS
  MAX.FREQ.U <- MGm[,lapply(.SD, max), by = NUCLEOTIDE, .SDcols = "frequency"][NUCLEOTIDE %in% "U"][["frequency"]]
  FIRST.NUC <- MGm[NUCLEOTIDE %in% "U" & POSITION %in% 1]
  X.BREAKS <- c(-PIRNA.SIZE,1,PIRNA.SIZE)
  
  ## PLOT
  message(" plotting")
  SMG <- ggplot() + theme_bw() +
    
    ## PIRNA LINES
    geom_vline(xintercept = c(-PIRNA.SIZE,1,PIRNA.SIZE), colour = "black", linetype = "dashed") +
    
    ## DATA
    geom_line(data = MGm, aes_string(x = "POSITION", y = PLOT.BY.VALUE, colour = "NUCLEOTIDE", group = "NUCLEOTIDE")) +
    
    ## U FREQ MAX
    geom_text_repel(data = FIRST.NUC, aes(x = POSITION, y = eval(parse(text = PLOT.BY.VALUE)), 
                                          label = round(eval(parse(text = PLOT.BY.VALUE)), digits = 1)), 
                    box.padding = 2, nudge_x = 2) +
    
    ## TITLES
    ggtitle(SAMPLE.NAME) + xlab("distance from 1N (nt)") + ylab("nucleotide frequency (%)") +
    
    ## SCALES
    scale_colour_manual(values = NUC.COLORS) +
    scale_x_continuous(breaks = X.BREAKS) +
    scale_y_continuous(limits = Y.LIMITS, breaks = Y.BREAKS) +
    
    ## THEME
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = LEGEND.POSITION,
          aspect.ratio = ASPECT.RATIO,
          plot.title = element_text(family = FAM, size = XYT, colour = TCOL),
          axis.text = element_text(family = FAM, size = XYT, colour = TCOL),
          axis.title = element_text(family = FAM, size = XYT, colour = TCOL))
  
  message("Done.")
  return(SMG)
  
  
}
