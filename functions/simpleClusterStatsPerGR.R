simpleClusterStatsPerGR <- function(
                                    ## INPUT      
                                    GR = NULL,
                                    RANGE.NAME = NULL,
                                    CLUSTER.NAME.COLUMN = "seqnames",
                                    SAMPLE.SEP = "_"){
  
  
  ## AUTHOR: Pavol Genzor
  ## Use: Calcualte mean and median stats per cluster
  ## 07.07.2021; Version 2; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("plyr"); library("dplyr");.
    library("GenomicRanges"); library("reshape")})
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GR!")
  if(is.null(RANGE.NAME)) stop("Please provide a RANGE.NAME!")
  if(isFALSE("MULT" %in% colnames(mcols(GR)))) stop("GR must have MULT column!")

  message(" FUN: cluster stats")
  a.DT <- as.data.table(GR)
  
  ## calculate mean
  cluster.mean.MULT <- a.DT[,lapply(.SD, function(i){mean(i)}), by = CLUSTER.NAME.COLUMN, .SDcols = "MULT"]
  colnames(cluster.mean.MULT) <- c("CLUSTER",paste("mean.mult",RANGE.NAME, sep = SAMPLE.SEP))
  
  ## calculate median
  cluster.median.MULT <- a.DT[,lapply(.SD, function(i){median(as.double(i))}), by = CLUSTER.NAME.COLUMN, .SDcols = "MULT"]
  colnames(cluster.median.MULT) <- c("CLUSTER",paste("median.mult",RANGE.NAME, sep = SAMPLE.SEP))
  
  ## compile
  cluster.stat.MULT <- setDT(join_all(dfs = list(cluster.mean.MULT, cluster.median.MULT), by = "CLUSTER", type = "full" ))
  
  ## return
  return(cluster.stat.MULT)
}




