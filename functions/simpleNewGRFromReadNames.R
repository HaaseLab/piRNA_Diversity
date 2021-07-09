simpleNewGRFromReadNames <- function(## INPUT
                                        GR = NULL,
                                        
                                        ## OPTIONS
                                        INCLUDE.N = FALSE,
                                        
                                        ## OPTIONS
                                        READ.NAMES.SPLIT = "_"){
  
  ## AUTHOR: Pavol Genzor
  ## Use: To make new GR from read names
  ## 07.07.21; Version 3; refined

  ## LIBRARIES
  suppressPackageStartupMessages({library(data.table);library(dplyr);library(plyr);library(GenomicRanges)})
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GR with read names !")
  if(!"NH" %in% colnames(mcols(GR)) || !"MULT" %in% colnames(mcols(GR))) 
    stop("Make sure mcols(GR) have NH and MULT columns !")
  
  message(" FUN: decoding names")
  NAMES.L <- tstrsplit(names(GR),split = READ.NAMES.SPLIT)
  AGR <- GRanges(seqnames = NAMES.L[[1]], 
                 ranges = IRanges(start = as.numeric(NAMES.L[[2]]), end = as.numeric(NAMES.L[[3]])), 
                 strand = NAMES.L[[4]])
  names(AGR) <- names(GR)
  mcols(AGR)[["CLUSTER"]] <- as.character(seqnames(GR))
  
  ## Option: Include sample count
  if(isTRUE(INCLUDE.N)){
    mcols(AGR)[["N"]] <- nth(tstrsplit(NAMES.L[[5]], split = "N"),2) }
  
  mcols(AGR)[["NH"]] <- mcols(GR)[["NH"]]
  mcols(AGR)[["MULT"]] <- mcols(GR)[["MULT"]]
  
  ## RETURN
  return(AGR)
}
