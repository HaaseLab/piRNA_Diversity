simpleGRFilter <- function(
                          ## INPUT
                          GR = NULL,
                          RANGE.NAME = NULL,
                          
                          ## FILTERS
                          NH.TAG = NULL,
                          SIZE.RANGE = NULL,
                          STRAND = NULL,
                          
                          # OPTIONS
                          SEQNAMES.SPLIT = ":",
                          L.STRAND.POSITION = 3){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Filter a single GR
  ## 06.28.21; Version 2, refined
  
  ## Libraries
  suppressPackageStartupMessages({library("data.table");library("plyr");
    library("Biostrings");library("GenomicRanges")})

  ## Check input
  if(is.null(GR)) stop("Please provide a GR !")
  if(ncol(mcols(GR)) == 0) stop("GR needs to have mcols() !")
  if(!"NH" %in% colnames(mcols(GR))) stop("mcols() need to have NH column !")
  if(!"MULT" %in% colnames(mcols(GR))) stop("mcols() need to have MULT column !")
  FILTER.LIST <- c("NH.TAG"=NH.TAG,"SIZE.RANGE" = SIZE.RANGE,"STRAND" = STRAND)
  if(is.null(FILTER.LIST)) stop("Please provide at least one filter criteria !")
  
  ## PROCEED
  message(paste0(" FUN: filtering..."))
  
  ## MAPPING
  if(!is.null(NH.TAG)){
    message(paste0("\tselect mappers: ",NH.TAG))
    GR <- GR[mcols(GR)[["NH"]] %in% NH.TAG]
  } else { message("\tselect mappers: all") }
  
  ## SIZE RANGE
  if(!is.null(SIZE.RANGE)){
    message(paste0("\tsize range used: ", min(SIZE.RANGE),"-", max(SIZE.RANGE)))
    SIZE.RANGE <- seq(min(SIZE.RANGE),max(SIZE.RANGE))
    GR <- GR[width(GR) %in% SIZE.RANGE]
  } else { message("\tall sizes") }
  
  ## STRAND
  if(!is.null(STRAND)){
    message("\tmatched strand: YES")
    mcols(GR)[["STRAND"]] <- nth(tstrsplit(seqnames(GR), split=SEQNAMES.SPLIT),L.STRAND.POSITION)
    GR <- GR[strand(GR) == mcols(GR)[["STRAND"]]]
  } else { message("\tmatched strand: NO") }
  
  ## RETURN
  return(GR)
}
