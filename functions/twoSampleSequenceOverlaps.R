twoSampleSequenceOverlaps <- function(
                                      ## INPUT
                                      GRL=NULL,
                                      BSSPECIES=NULL,
                                      
                                      ## SETTINGS
                                      MC.CORES=4,
                                      
                                      ## OUTPUT
                                      REMOVE.SEQ=TRUE
                                      ){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Calculate two sample sequence overlaps
  ## 06.28.21; Version 3; refined

  ## load libraries
  suppressPackageStartupMessages({library("data.table");library("parallel");
    library("GenomicRanges");library("GenomicAlignments");library("BSgenome.Hsapiens.UCSC.hg38");
    library("BSgenome.Dmelanogaster.UCSC.dm6");library("BSgenome.Mmusculus.UCSC.mm10")})
  
  ## INPUT CHECKING
  if(is.null(GRL)) stop("Please provide named GRanges list !")
  if(is.null(BSSPECIES)) stop("Please provide BSSPECIES name !")
  if(!isTRUE(class(GRL) %in% "list")) stop("GRL must be a named LIST !")
  if(is.null(names(GRL))) stop("GRP must be a NAMED list !!!")
  if(!isTRUE(unique(unlist(lapply(GRL,class)) %in% "GRanges"))) stop("All the components must be GRanges !")
  
  message("Processing ...")
  
  ## GET SEQUENCES
  message("\textracting sequence")
  GRLS <- mclapply(names(GRL), mc.cores = MC.CORES, function(i){ 
    BSgenome::getSeq(eval(parse(text = BSSPECIES)),GRL[[i]]) })
  names(GRLS) <- names(GRL)
  
  message("\tchecking sequence duplication")
  for(i in names(GRLS)) {print(table(duplicated(GRLS[[i]]))  )}
  
  message("\t\nadding sequence to GRanges")
  for(i in names(GRLS)) {mcols(GRL[[i]])[["seq"]] <- GRLS[[i]]}
  
  message("\tmaking name map and assigning")
  TEMP.NAMES <- c("FA","FB")
  names(GRLS) <- TEMP.NAMES
  NAME.MAP <- data.table(tempName = TEMP.NAMES, 
                         realName = names(GRL),
                         Input_N = lapply(GRLS,length))
  ##
  ## ALL COMMON
  ##
  
  message(" OVERLAPPING: all common")
  SEQ.TWO.COMMON <- Reduce(BiocGenerics::intersect,GRLS)
  
  message("\tfiltering orginal GRL")
  TWO.COMMON.GRL <- lapply(names(GRL),function(i){ 
    A.GR <- GRL[[i]]
    N.GR <- A.GR[mcols(A.GR)[["seq"]] %in% SEQ.TWO.COMMON] 
    
    ## remove sequence
    if(isTRUE(REMOVE.SEQ)){mcols(N.GR)[["seq"]] <- NULL}
    
    ## return
    return(N.GR) })
  names(TWO.COMMON.GRL) <- paste0(names(GRL),"_COMMON")
  
  ##
  ## ALL EXCLUSIVE
  ##
  
  message(" OVERLAPPING: all exclusive")
  SEQ.TWO.EXCLUSIVE <- lapply(TEMP.NAMES,function(i){
    TEMP.NAME.ORDER <- c(i,TEMP.NAMES[!TEMP.NAMES %in% i])
    AN.EXCLUSIVE.SEQ <- Reduce(BiocGenerics::setdiff, GRLS[TEMP.NAME.ORDER])
    return(AN.EXCLUSIVE.SEQ)})
  names(SEQ.TWO.EXCLUSIVE)  <- paste0(names(GRL))
  
  message("\tfiltering orginal GRL")
  TWO.EXCLUSIVE.GRL <- lapply(names(GRL), function(i){
    A.GR <- GRL[[i]]
    N.GR <- A.GR[mcols(A.GR)[["seq"]] %in% SEQ.TWO.EXCLUSIVE[[i]]] 
    
    ## remove sequence
    if(isTRUE(REMOVE.SEQ)){mcols(N.GR)[["seq"]] <- NULL}
    
    ## return
    return(N.GR)})
  names(TWO.EXCLUSIVE.GRL)  <- paste0(names(GRL),"_EXCLUSIVE")
  
  ##
  ## COMPILE SUMMARIES
  ##
  
  message(" preparing summaries")
  NAME.MAP[["Common_N"]] <- length(SEQ.TWO.COMMON)
  NAME.MAP[["Exclusive_N"]] <- lapply(SEQ.TWO.EXCLUSIVE,length)
  
  ## CREATE A RESULTS LIST
  RES.L <- list(TWO.COMMON.GRL, TWO.EXCLUSIVE.GRL, NAME.MAP)
  names(RES.L) <- c("Common","Exclusive","SampleInfo"); message("Completed.")
  message(paste0("\n","Samples Info"))
  print(NAME.MAP)
  
  ## RETURN
  return(RES.L) 
}
