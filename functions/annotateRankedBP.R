annotateRankedBP <- function(
                              ## INPUT
                              ANN.TAB.L = NULL,
                              
                              ## SETTINGS
                              ASPECT.RATIO = 2,
                              COORD.FLIP=FALSE,
                              Y.LIMS = c(-100,100),
                              Y.LAB = "Fraction (%)",
                              FAM = "Helvetica", TCOL = "black", XYT = 8,
                              RETURN.ALL = FALSE){
  
  ## Pavol Genzor
  ## To plot sense and antisense annotation categories
  ## 07.08.21; Version 2; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("ggplot2")})
  
  ## INPUT
  if(is.null(ANN.TAB.L)) stop("Please provide ANN.TAB.L !")
  if(!is.list(ANN.TAB.L)) stop("Please provide list of tables !")
  
  ## Combine samples into directional tables
  DT.UP <- setDT(ANN.TAB.L[["DTM.UP"]])
  #DT.UP[["SAMPLE"]] <- factor(DT.UP[["SAMPLE"]], levels = unique(DT.UP[["SAMPLE"]]))
  DT.UP[["NAME"]] <- factor(DT.UP[["NAME"]], levels = unique(DT.UP[["NAME"]]))
  
  DT.DOWN <- setDT(ANN.TAB.L[["DTM.DOWN"]])
  #DT.DOWN[["SAMPLE"]] <- factor(DT.DOWN[["SAMPLE"]], levels = unique(DT.DOWN[["SAMPLE"]]))
  DT.DOWN[["NAME"]] <- factor(DT.DOWN[["NAME"]], levels = unique(DT.DOWN[["NAME"]]))
  
  ## All data table
  DT <- rbindlist(list(DT.UP,DT.DOWN))

  ## Rotate
  if(isTRUE(COORD.FLIP)){COORD=coord_flip()
  } else {COORD=NULL}
  
  ## Plot
  GGBP <- ggplot() + theme_bw() +
    ## DATA
    geom_bar(data = DT.UP, aes(x = SAMPLE, y = UD.VALUE, fill = NAME), 
             stat = "identity", 
             #fill = DT.UP[["COLS"]], 
             colour = NA) + 
    geom_bar(data = DT.DOWN, aes(x = SAMPLE, y = UD.VALUE, fill = NAME), 
             stat = "identity", 
             #fill = DT.DOWN[["COLS"]], 
             colour = NA) + 
    
    ## LABELS
    #geom_text(data = DT.UP, aes(x = SAMPLE, y = UD.VALUE, label = NAME), 
    #          position = position_stack(vjust = 0.5), 
    #          colour=TCOL, family=FAM, size = 2) +
    #geom_text(data = DT.DOWN[NAME %in% "OTHER"], aes(x = SAMPLE, y = UD.VALUE, label = NAME), 
    #          position = position_stack(vjust = 0.5), 
    #          colour="white", family=FAM, size = 2) +
    
    ## LABS
    xlab("") + ylab(paste0(Y.LAB,"\n(-) antisense / (+) sense")) + 
    ggtitle(paste0("Annotation Barplot\n",
                   "used: ",unique(DT[["dataType"]]),"\n",
                   "filters: ",unique(DT[["FILTERS"]]))) +
    
    ## LINE
    geom_hline(yintercept = 0, size = 0.5, colour = "black" ) +
    
    ## SCALES
    scale_y_continuous(limits = Y.LIMS, breaks = seq(-100,100,10)) +
    scale_fill_manual(values = union(DT.UP[["COLS"]], DT.DOWN[["COLS"]])) +

    ## ROTATE    
    COORD +

    ## THEME
    theme(aspect.ratio = ASPECT.RATIO, 
          panel.border = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_text(family = FAM, color = TCOL, size = XYT, 
                                     angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(family = FAM, color = TCOL, size = XYT),
          axis.title = element_text(family = FAM, color = TCOL, size = XYT)) 

  
  ## RESULTS
  RES.L <- list("plotData" = DT, "plot" = GGBP)
  
  ## RETURN
  if(isTRUE(RETURN.ALL)) {
    return(RES.L) } else { return(RES.L[["plot"]]) }
}
