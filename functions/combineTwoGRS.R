combineTwoGRS <- function(
                            ## INPUT
                            GRL=NULL,
                            DUPLICATE.NAMES=NULL, 
                            
                            ## SETTINGS
                            MC.CORES=3
                            ){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Combine two Genomic Ranges
  ## 06.28.21; Version 3; refined
  ## NOTE: This function takes time
  
  ## LIBRARIES
  suppressPackageStartupMessages({library("data.table");library("plyr");
    library("dplyr");library("GenomicRanges");library("stringi");library("parallel")})

  ## INPUT CHECKING
  if(is.null(GRL)) stop("Provide named list of Genomic Ranges")
  if(isFALSE(class(GRL) %in% c("list"))) stop("Input must be a list of GRanges")
  if(isTRUE(is.null(names(GRL)))) stop("This list needs to be named")
  if(is.null(DUPLICATE.NAMES)) stop("Provide names for two samples")
  if(!length(DUPLICATE.NAMES)==2) stop("Provide exactly TWO ranges")
  
  message("Proceeding ...")
  
  message(" combining duplicates")
  D.REPS.GR <- c(GRL[[DUPLICATE.NAMES[1]]],
                 GRL[[DUPLICATE.NAMES[2]]])
  
  message(" making a data.table")
  D.REP.DT <- as.data.table(D.REPS.GR)
  
  message(" removing width")
  D.REP.DT[[c("width")]] <- NULL
  
  message(" creating a unique ID")
  D.REP.DT[["id"]] <- stringi::stri_c(D.REP.DT[["seqnames"]],
                                      D.REP.DT[["start"]], 
                                      D.REP.DT[["end"]], 
                                      D.REP.DT[["strand"]], sep = "__")
  
  message(" summarizing mcols")
  ## Sumarizing happens over UNIQUE ID here
  USQ.COUNT <- D.REP.DT[,.N, by = "id"]
  USQ.NH <- D.REP.DT[,lapply(.SD, max), by = "id", .SDcols = c("NH")]
  USQ.MULT <- D.REP.DT[,lapply(.SD, sum), by = "id", .SDcols = c("MULT")]
  
  message(" combining summary results") ## SLOW
  NF <- setDT(plyr::join_all(dfs = list(USQ.COUNT,USQ.NH,USQ.MULT), 
                             by = "id", type = "full"))
  
  message(" expanding the ID")
  EXPANDED.ID.L <- mclapply(seq(1,4), mc.cores = MC.CORES, function(i){
    dplyr::nth(data.table::tstrsplit(NF[["id"]],split="__",fixed=TRUE),i) })
  names(EXPANDED.ID.L) <- c("chr","start","end","strand")
  
  message(" adding expanded ID items to the table")
  NF[["chr"]] <- EXPANDED.ID.L[["chr"]]
  NF[["start"]] <- EXPANDED.ID.L[["start"]]
  NF[["end"]] <- EXPANDED.ID.L[["end"]]
  NF[["strand"]] <- EXPANDED.ID.L[["strand"]]
  
  message(" converting to GRanges")
  NF.GR <- GenomicRanges::makeGRangesFromDataFrame(df = NF, keep.extra.columns = TRUE)
  
  message(" cleaning up")
  mcols(NF.GR)[["id"]] <- NULL
  
  message("done!")
  return(NF.GR)
}
