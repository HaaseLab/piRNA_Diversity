simpleMultBarPlot <- function(GR = NULL,
                              RANGE.NAME = NULL,
                              
                              ## SETTINGS
                              COLORS = c("orange","white"),
                              FAM = "Helvetica", 
                              TCOL = "black", 
                              XYT = 12,
                              ASPECT.RATIO = 0.1,
                              BAR.WIDTH = 0.75,
                              
                              ## RETURNS
                              RETURN.ALL = FALSE){
  
  
  ## AUTHOR: Pavol Genzor
  ## Use: Generate barplot from multiplicities (to supplement STEPS plot)
  ## 06.28.2021; Version 2; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr");
    library("ggplot2");library("ggpubr")})
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GRL !")
  if(is.null(RANGE.NAME)) stop("Please provide a RANGE.NAME !")
  
  ## Process data
  DT <- as.data.table(mcols(GR)[["MULT"]])[,.N, by = "V1"]
  colnames(DT) <- c("MULT","SEQ")
  DT[,READ := MULT*SEQ]
  DT[,fSEQ := (SEQ/sum(SEQ))]
  DT[,fREAD := (READ/sum(READ))]
  setorderv(x = DT, cols = "MULT", order = -1)
  DT[,afSEQ := Reduce('+', fSEQ, accumulate = T)]
  DT[,afREAD := Reduce('+', fREAD, accumulate = T)]
  
  ## Reduce to relevant and melt
  DTS <- rbindlist(l = list(DT[MULT %in% 1],DT[!MULT %in% 1][,lapply(.SD,sum)]))
  DTS[["MULT"]] <- factor(x = c("one","rest"), levels = c("rest","one"))
  DTM <- melt.data.table(DTS,id.vars = "MULT", 
                         measure.vars = c("fSEQ","fREAD"), 
                         variable.name = "type", 
                         value.name = "frequency")
  DTM[["type"]] <- factor(x = DTM[["type"]] , levels = c("fREAD","fSEQ"))
  
  ## Plot
  GGBAR <- ggplot() + theme_pubclean() +
    ## Data
    geom_bar(data = DTM, aes(x = type, y = frequency, fill = MULT),
             colour = "black",stat = "identity", width = BAR.WIDTH) + 
    ## Label
    geom_text(data = DTM, aes(x = type, y = frequency, label = 100*round(frequency,digits = 3)), 
              position = position_stack(vjust = 0.5), 
              colour=TCOL, family=FAM) +
    ## Colors
    scale_fill_manual(values = COLORS) +
    ## Coordinates & title
    coord_flip() + ggtitle(paste0(RANGE.NAME)) +
    ## Theme
    theme(aspect.ratio = ASPECT.RATIO, panel.grid = element_blank(),
          axis.text.x = element_text(family = FAM, colour = TCOL, size = XYT, 
                                     angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(family = FAM, colour = TCOL, size = XYT))
  
  ## Results
  RES.L <- list("plotData" = DTM, "ggplot" = GGBAR)
  if(isTRUE(RETURN.ALL)) { return(RES.L) } 
  else { return(RES.L[["ggplot"]]) }
}
