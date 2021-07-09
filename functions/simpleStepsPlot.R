simpleStepsPlot <- function(
                            ## INPUT
                            GRL = NULL,
                            
                            ## OPTIONS
                            Y.PPM = FALSE,
                            
                            ## SETTINGS
                            COLORS = scale_color_brewer(palette = "Dark2"),
                            LWD.STEP = 1.75, 
                            Y.TRANS = "log10",
                            X.LIMS = NULL, 
                            Y.LIMS = NULL,
                            FAM = "Helvetica", 
                            TCOL = "black", 
                            XYT = 12,
                            ASPECT.RATIO = 1,
                            LEGEND.POSITION = "bottom",
                            
                            ## RETURN
                            RETURN.ALL = FALSE){
  
  ## AUTHOR: Pavol Genzor & Daniel Stoyko
  ## Use: Calculate and plot steps plots
  ## 06.28.21; Version 3; refined

  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr");
    library("ggplot2");library("ggpubr")})
  
  ## INPUT CHECKING
  if(is.null(GRL)) stop("Please provide a GRL !")
  if(!is.list(GRL)) stop("GRL needs to be a list !")
  if(is.null(names(GRL))) stop("GRL has to be named !")
  
  ## LOOP THROUGH SAMPLES
  SAMPLE.L <- lapply(names(GRL), function(s){
    
    ## Subset
    GR <- GRL[[s]]
    if(!"MULT" %in% colnames(mcols(GR))) stop("MULT column is missing !")
    
    ## Process
    DT <- as.data.table(mcols(GR)[["MULT"]])[,.N, by = "V1"]
    colnames(DT) <- c("MULT","SEQ")
    DT[,READ := MULT*SEQ]
    DT[,fSEQ := (SEQ/sum(SEQ))]
    DT[,fREAD := (READ/sum(READ))]
    setorderv(x = DT, cols = "MULT", order = -1)
    DT[,afSEQ := Reduce('+', fSEQ, accumulate = T)]
    DT[,afREAD := Reduce('+', fREAD, accumulate = T)]
    DT[["ORG.MULT"]] <- DT[["MULT"]]
    DT[["SAMPLE"]] <- s
    
    ## CONVERT ABUNDANCE TO PPM
    if(isTRUE(Y.PPM)){ 
      message("... MULT to PPM conversion (adjust Y.LIMS) ")
      DT[["MULT"]] <- (DT[["MULT"]] / sum(DT[["READ"]])) * 1000000}
    
    ## Return
    return(DT)})
  
  ## Y LABELS: MULT vs PPM
  if(isTRUE(Y.PPM)){ Y.LAB <- ylab("MULT.PPM") } else { Y.LAB<- ylab("MULT") }
  
  ## Combine samples
  LDT <- rbindlist(SAMPLE.L)
  
  ## Melt & Prepare
  DTM <- melt.data.table(data = LDT, 
                         id.vars = c("MULT","SAMPLE"), 
                         measure.vars = c("afSEQ","afREAD"),
                         value.name = "FRACTION", 
                         variable.name = "TYPE")
  DTM[["SAMPLE"]] <- factor(DTM[["SAMPLE"]], levels = unique(DTM[["SAMPLE"]]))
  
  ## PLOT
  SGG <- ggplot() + theme_bw() +
    ## data
    geom_step(data = DTM, aes(x = FRACTION, y = MULT, alpha = TYPE, colour = SAMPLE), 
              lwd = LWD.STEP) +
    COLORS + Y.LAB +
    ## scales
    scale_alpha_discrete(range=c(0.5, 1)) +
    scale_y_continuous(limits = Y.LIMS, trans = Y.TRANS) +
    scale_x_continuous(limits = X.LIMS, breaks = seq(0,1,0.1)) +
    annotation_logticks(sides = "l", ) +
    ## theme
    theme(aspect.ratio = 1, legend.position = "bottom",
          panel.grid = element_blank(), panel.border = element_blank(),
          axis.title = element_text(family = FAM, colour = TCOL, size = XYT),
          axis.text = element_text(family = FAM, color = TCOL, size = XYT),
          axis.ticks.y = element_blank())
  
  ## results
  RES.L <- list("sampleList" = SAMPLE.L, "plotData" = DTM, "plot" = SGG)
  
  ## return
  if(isTRUE(RETURN.ALL)){ return(RES.L) } else { return(RES.L[["plot"]]) }
}
