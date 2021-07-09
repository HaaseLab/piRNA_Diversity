simpleSRViolin <- function(
                              ## INPUT
                              GRL = NULL,
                              
                              ## OPTIONS
                              Y.PPM = FALSE,
                              
                              ## SETTINGS
                              SAMPLE.ORDER = names(GRL),
                              SIZE.RANGE = c(18,50), 
                              NH.TAG = NULL,
                              TOTAL.SAMPLE="Total",
                              
                              ## OUTPUT
                              ADD.TABLE = TRUE,
                              SIG.FIGURES = 2,
                              SOURCE.DIR = NULL,
                              RETURN.ALL = FALSE,
                              
                              ## PLOT SETTINGS
                              Y.LIMS = NULL,
                              VIOLIN.ADJUST = 1, 
                              VIOLIN.SCALE = "width",
                              ASPECT.RATIO = 1,
                              LEGEND.POSITION = "right",
                              TCOL = "black",FAM = "Helvetica",XYT = 12,
                              TRANS = "log10",PLOT.COLORS = c("orange1","grey70")){
  
  ## AUTHOR: Pavol Genzor & Daniel Stoyko
  ## Use: Calculate and plot steps plot
  ## 06.28.21; Version 6; fixed and refined

  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("GenomicRanges")
    library("ggplot2"); library("tidyverse"); library("reshape"); 
    library("gridExtra");library("ggpubr")})
  
  ## INPUT CHECKING
  if(is.null(GRL)) stop("Please provide a named GRL object !")
  if(isFALSE(is.list(GRL))) stop("Input is not a GRL - Please proide a list of GRanges !")
  if(is.null(names(GRL))) stop("Please make sure that GRL is named !")
  if(is.null(SOURCE.DIR)) stop("Please provide a SOURCE.DIR with functions !")
  if(!"MULT" %in% colnames(mcols(GRL[[1]]))) stop("GRL needs to have MULT column !")
  if(!"NH" %in% colnames(mcols(GRL[[1]]))) stop("GRL needs to have NH column !")
  
  ## FUNCTION
  source(paste0(SOURCE.DIR,"simpleGRFilter.R"))
  
  ## PROCEED
  message("Processing...")
  
  ##
  ## DATA LOOP
  ##
  
  GRL.DT.L <- lapply(names(GRL), function(i){
    
    ## SUBSET
    GR <- GRL[[i]]
    
    ## FILTER
    GR <- simpleGRFilter(GR = GR, 
                         RANGE.NAME = i,
                         NH.TAG = NH.TAG, 
                         SIZE.RANGE = SIZE.RANGE)
    
    ## PREP THE TABLE
    GR.DT <- as.data.table(GR)
    NDT <- GR.DT[,.N, by = "MULT"]
    setorderv(x = NDT, cols = "MULT")
    colnames(NDT) <- c("MULT","SEQ")
    
    ## CALCULATE 
    NDT[, READ := MULT * SEQ]
    NDT[, SAMPLE := i]
    
    ## REARRANGE
    mNDT <- setDT(melt.data.table(data = NDT,
                                  id.vars = c("MULT","SAMPLE"), 
                                  variable.name = "TYPE", 
                                  value.name = "COUNT"))
    ## RETURN
    return(mNDT) })

  ## COMNBINE AND ORGANIZE FOR PLOTTING
  message(" calculating & plotting")
  GR.DT <- rbindlist(GRL.DT.L)
  GR.DT[["SAMPLE"]] <- factor(GR.DT[["SAMPLE"]], levels = SAMPLE.ORDER)
  GR.DT[["TYPE"]] <- factor(GR.DT[["TYPE"]], levels = c("SEQ","READ"))
  GR.DT[["GROUP"]] <- paste(GR.DT[["SAMPLE"]],GR.DT[["TYPE"]], sep = "_")
  for(i in unique(GR.DT[["GROUP"]])) set(x = setDT(GR.DT), i = which(GR.DT[["GROUP"]] %in% i), 
                                         j = "GROUP_SUM", value = sum(GR.DT[GROUP %in% i][["COUNT"]]))
  setDT(GR.DT)

    ##
  ## TABLE FOR PLOT
  ##
  
  TAB.DT <- GR.DT[, sum(.SD), by = c("SAMPLE","TYPE"), .SDcols = "COUNT"]
  for(t in unique(TAB.DT[["TYPE"]])) set(x = setDT(TAB.DT), 
                                         i = which(TAB.DT[["TYPE"]] %in% t), 
                                         j = "FRACTION", 
                                         value = round(TAB.DT[TYPE %in% t][["V1"]]/
                                                         TAB.DT[TYPE %in% t & SAMPLE %in% TOTAL.SAMPLE][["V1"]],
                                                       digits = SIG.FIGURES))
  TAB.W <- dcast.data.table(data = TAB.DT, formula = TYPE ~ SAMPLE, value.var = "FRACTION")
  
  ## CONVERT TO PPM
  if(isTRUE(Y.PPM)){
    message("... MULT to PPM conversion (adjust Y.LIMS) ")
    NOMINATOR <- GR.DT[["MULT"]]
    DENOMINATOR <- unique(GR.DT[GROUP %in% paste(TOTAL.SAMPLE,"READ",sep="_")][["GROUP_SUM"]])
    PPMs <- ((NOMINATOR/DENOMINATOR) * 1000000)
    GR.DT[["MULT"]] <- PPMs
    Y.LAB <- ylab("MULT.PPM")
  } else { Y.LAB <- ylab("MULT") }

  ##
  ## VIOLIN PLOT
  ##
  
  GG.VIOL <-  ggplot() + theme_pubclean() +
    ## DATA
    geom_violin(data = GR.DT, aes(x = SAMPLE, y = MULT, weight = COUNT/GROUP_SUM, fill = TYPE),
                color="black", lwd = 0.25, adjust= VIOLIN.ADJUST, scale = VIOLIN.SCALE, 
                position = position_dodge(width = 1)) +
    
    ## TITLE
    ggtitle(paste0("Violins settings\n",
                   paste0("size range: ", min(SIZE.RANGE)," - ", max(SIZE.RANGE)),"\n",
                   paste0("used NH tags: ", ifelse(is.null(NH.TAG),"all",NH.TAG) ))) + Y.LAB +
    
    ## SCALES
    scale_fill_manual(values = PLOT.COLORS) +
    scale_y_continuous(trans = TRANS, limits = Y.LIMS) +
    annotation_logticks(sides = "l") +
    theme(aspect.ratio = ASPECT.RATIO, 
          axis.ticks.y = element_blank(),
          legend.position = LEGEND.POSITION,
          panel.grid = element_blank(), 
          axis.text = element_text(family = FAM, color =TCOL, size = XYT))
  
  ## CONVERT INTO GROB
  if(isTRUE(ADD.TABLE)){
    TBL <- tableGrob(d = TAB.W, theme = ttheme_minimal(), rows = NULL)
    GG.VIOL <- grid.arrange(GG.VIOL,TBL, nrow = 2, heights = c(5,1)) }
  
  ## RETURN
  message("Done.")
  if(isTRUE(RETURN.ALL)){
    RES.L <- list("plotData" = GR.DT,"plot" = GG.VIOL)
    return(RES.L)
  } else { return(GG.VIOL) }
}
