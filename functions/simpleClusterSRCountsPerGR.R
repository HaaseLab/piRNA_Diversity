simpleClusterSRCountsPerGR <- function(
                                      ## INPUT
                                      GR = NULL,
                                      RANGE.NAME = NULL,
                                      CLUSTER.NAME.COLUMN = "seqnames",
                                      SAMPLE.SEP = "_"){
  
  
  ## AUTHOR: Pavol Genzor
  ## Use: Retrieve sequence and read counts per cluster
  ## 07.07.2021; Version 2; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("plyr"); library("dplyr");.
    library("GenomicRanges"); library("reshape")})
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GR!")
  if(is.null(RANGE.NAME)) stop("Please provide a RANGE.NAME!")
  if(isFALSE("MULT" %in% colnames(mcols(GR)))) stop("GR must have MULT column!")

  message(" FUN: cluster counts")
  a.DT <- as.data.table(GR)
  
  ## sequences
  seq.count.DT <- a.DT[,.N, by = CLUSTER.NAME.COLUMN]
  colnames(seq.count.DT) <- c("CLUSTER", paste("seq.count", RANGE.NAME, sep = SAMPLE.SEP))
  seq.count.DT[,paste("seq.total",RANGE.NAME, sep = SAMPLE.SEP) := length(GR)]
  
  ## reads
  read.count.DT <- a.DT[, lapply(.SD, function(i){sum(i)}), by = CLUSTER.NAME.COLUMN, .SDcols = "MULT"]
  colnames(read.count.DT) <- c("CLUSTER",paste("read.count", RANGE.NAME, sep = SAMPLE.SEP))
  read.count.DT[,paste("read.total",RANGE.NAME, sep = SAMPLE.SEP) := sum(mcols(GR)[["MULT"]])]
  
  ## join
  sr.count.DT <- setDT(join_all(dfs = list(seq.count.DT,read.count.DT), by = "CLUSTER", type = "full"))
  for (j in names(sr.count.DT)) set(sr.count.DT,which(is.na(sr.count.DT[[j]])),j,0)
  
  ## return
  return(sr.count.DT)
}




