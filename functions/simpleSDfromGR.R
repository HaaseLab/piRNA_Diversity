simpleSDfromGR <- function(
                            ## INPUT
                            GR=NULL,
                            SAMPLE.NAME=NULL,
                            
                            ## OPTIONS
                            USE.READS=TRUE,
                            PLOT.FREQ=TRUE,
                            YLIMS=NULL,
                            
                            ## SETTINGS
                            ASPECT.RATIO=1,
                            X.BREAKS.BY=1,
                            BAR.FILL="grey80",
                            BAR.LINE="black",
                            XYT=8,
                            FAM="Helvetica",
                            TCOL="black",
                            RETURN.ALL=FALSE){
  
  
  ## AUTHOR: Pavol Genzor
  ## 06.30.21; Version 4; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({library(data.table);library(ggplot2);library(GenomicRanges)})

  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GR  made from .bam file !")
  if(is.null(SAMPLE.NAME)) stop("Please provide SAMPLE.NAME !")
  if(!"MULT" %in% colnames(mcols(GR))) stop("The GRanges must have MULT column!")
  
  ## OPTIONS
  PLOT.ON.Y=ifelse(isTRUE(PLOT.FREQ),"FREQ","COUNT")
  
  ## MAKE TABLE
  GR.DT  <- as.data.table(GR)
  
  ## COUNT SEQ OR READ
  if(isTRUE(USE.READS)){
    message("using reads")
    METRIC.USED = "Read"
    SD.DT <- GR.DT[, lapply(.SD,sum), by = "width", .SDcols = "MULT"] } 
  else {
    message("using sequences")
    METRIC.USED = "Sequence"
    SD.DT <- GR.DT[,.N, by = "width"] }
  
  colnames(SD.DT) <- c("SIZE","COUNT")
  SD.DT[, FREQ := lapply(.SD,function(i){(i/sum(.SD))*100 }),.SDcols = c("COUNT")]
  SD.DT[["SIZE"]] <- as.numeric(SD.DT[["SIZE"]])
  
  ## PLOT
  SD.GG <- ggplot() + theme_pubclean() +
    geom_bar(data = SD.DT, aes_string(x = "SIZE", y = PLOT.ON.Y), 
             stat = "identity", colour = BAR.LINE, fill = BAR.FILL) + 
    scale_x_continuous(breaks = seq(15,50,X.BREAKS.BY)) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE), 
                       limits = YLIMS, 
                       breaks = seq(0,100,10)) +
    ggtitle(paste(METRIC.USED,SAMPLE.NAME, sep = "; ")) +
    theme(panel.grid = element_blank(), 
          title = element_text(family = FAM, size = XYT, colour = TCOL),
          axis.text = element_text(family = FAM, size = XYT, colour = TCOL),
          axis.title = element_text(family = FAM, size = XYT, colour = TCOL),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          aspect.ratio = ASPECT.RATIO)
  
  ## RETURN
  if(isTRUE(RETURN.ALL)){
    RES.L <- list(SD.DT,SD.GG)
    names(RES.L) <- c(paste0(METRIC.USED,"_data_table"),"ggplot")
    return(RES.L)
  } else { return(SD.GG) }
}
