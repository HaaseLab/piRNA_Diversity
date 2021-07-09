threeSampleSequenceOverlaps <- function(
                                        ## INPUT
                                        GRL=NULL,
                                        BSSPECIES=NULL,
                                        
                                        ## SETTINGS
                                        MC.CORES=4,
                                        
                                        ## OUTPUT
                                        REMOVE.SEQ=TRUE){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Calculate three sample sequence overlaps
  ## 06.25.21; Version 3; refined
 
  ## load libraries
  suppressPackageStartupMessages({library("data.table");library("parallel");
    library("GenomicRanges");library("GenomicAlignments");library("BSgenome.Hsapiens.UCSC.hg38");
    library("BSgenome.Dmelanogaster.UCSC.dm6");library("BSgenome.Mmusculus.UCSC.mm10")})

  #source("/Users/genzorp/Documents/GITHUB/piDiversity/r/extractFastqMultiplicity.R")
  
  ## INPUT CHECKING
  if(is.null(GRL)) stop("Please provide named GRanges list !")
  if(is.null(BSSPECIES)) stop("Please provide BSSPECIES name !")
  if(!isTRUE(class(GRL) %in% "list")) stop("GRL must be a named LIST !")
  if(is.null(names(GRL))) stop("GRP must be a NAMED list !!!")
  if(!isTRUE(unique(unlist(lapply(GRL,class)) %in% "GRanges"))) stop("All the components must be GRanges !")
  
  message("Processing ...")
  
  ## get sequences
  message("\textracting sequences")
  GRLS <- mclapply(names(GRL), mc.cores = MC.CORES, function(i){
    BSgenome::getSeq(eval(parse(text = BSSPECIES)),GRL[[i]]) })
  names(GRLS) <- names(GRL)
  
  message("\tchecking sequence duplication")
  for(i in names(GRLS)) {print(table(duplicated(GRLS[[i]]))  )}
  
  message("\n\tadding sequence to GRanges")
  for(i in names(GRLS)) {mcols(GRL[[i]])[["seq"]] <- GRLS[[i]]}
  
  message("\tmaking name map and assigning\n")
  TEMP.NAMES <- c("FA","FB","FC")
  names(GRLS) <- TEMP.NAMES
  NAME.MAP <- data.table(tempName = TEMP.NAMES,
                         realName = names(GRL),
                         Input_N = lapply(GRLS,length))
  ##
  ## ALL COMMON
  ##
  
  message(" OVERLAPPING: all common")
  SEQ.THREE.COMMON <- Reduce(BiocGenerics::intersect,GRLS)
  
  message("\tfiltering orginal GRL")
  THREE.COMMON.GRL <- lapply(names(GRL),function(i){
    GRL[[i]][mcols(GRL[[i]])[["seq"]] %in% SEQ.THREE.COMMON] })
  names(THREE.COMMON.GRL) <- paste0(names(GRL),"_COMMON")
  
  # remove seq column
  if(isTRUE(REMOVE.SEQ)){
    THREE.COMMON.GRL <- lapply(THREE.COMMON.GRL,function(i){
      GR <- i; mcols(GR)[["seq"]] <- NULL;return(GR)}) }
  
  ##
  ## ALL EXCLUSIVE
  ##
  
  message(" OVERLAPPING: all exclusive")
  SEQ.THREE.EXCLUSIVE <- lapply(TEMP.NAMES,function(i){
    TEMP.NAME.ORDER <- c(i,TEMP.NAMES[!TEMP.NAMES %in% i])
    AN.EXCLUSIVE.SEQ <- Reduce(BiocGenerics::setdiff, GRLS[TEMP.NAME.ORDER])
    return(AN.EXCLUSIVE.SEQ)})
  
  names(SEQ.THREE.EXCLUSIVE)  <- paste0(names(GRL))
  
  message("\tfiltering orginal GRL")
  THREE.EXCLUSIVE.GRL <- lapply(names(GRL), function(i){
    GRL[[i]][mcols(GRL[[i]])[["seq"]] %in% SEQ.THREE.EXCLUSIVE[[i]]] })
  names(THREE.EXCLUSIVE.GRL)  <- paste0(names(GRL),"_EXCLUSIVE")
  
  # remove seq column
  if(isTRUE(REMOVE.SEQ)){
    THREE.EXCLUSIVE.GRL <- lapply(THREE.EXCLUSIVE.GRL,function(i){
      GR <- i; mcols(GR)[["seq"]] <- NULL;return(GR)}) }
  
  ##
  ## ALL PAIRWISE
  ##
  
  message(" OVERLAPPING: all pairwise")
  PW.ORDER <- list(c("FA","FB"),c("FB","FC"),c("FC","FA"))
  names(PW.ORDER) <- c("AB","BC","CA")
  SEQ.THREE.PAIRWISE <- lapply(PW.ORDER,function(i){
    A.PAIRWISE.SEQ <- Reduce(BiocGenerics::intersect, list(GRLS[[i[1]]],GRLS[[i[[2]]]]))
    return(A.PAIRWISE.SEQ)})
  names(SEQ.THREE.PAIRWISE)  <- names(PW.ORDER)
  
  message("\tfiltering orginal GRL")
  PAIRWISE.GRL <- lapply(names(PW.ORDER), function(i){
    PAIR.REAL.NAMES <- unlist(lapply(PW.ORDER[[i]],function(j){NAME.MAP[tempName %in% j][["realName"]]}))
    PAIR.GRL <- lapply(PAIR.REAL.NAMES, function(k){
      GRL[[k]][mcols(GRL[[k]])[["seq"]] %in% SEQ.THREE.PAIRWISE[[i]]] })
    names(PAIR.GRL) <- paste0(PAIR.REAL.NAMES,"_",i,"-PAIRWISE")
    
    # remove seq column
    if(isTRUE(REMOVE.SEQ)){
      PAIR.GRL <- lapply(PAIR.GRL,function(i){
        GR <- i; mcols(GR)[["seq"]] <- NULL;return(GR)}) }
    
    return(PAIR.GRL) })
  names(PAIRWISE.GRL) <- names(PW.ORDER)
  
  ## COMPILE SUMMARIES
  
  message("\tpreparing summaries")
  NAME.MAP[["Common_N"]] <- length(SEQ.THREE.COMMON)
  NAME.MAP[["Exclusive_N"]] <- lapply(SEQ.THREE.EXCLUSIVE,length)
  NAME.MAP[["Pairwise"]] <- names(PW.ORDER)
  NAME.MAP[["Pairwise_N"]] <- lapply(SEQ.THREE.PAIRWISE,length)
  
  ## COMPILE RESULTS & RETURN
  
  RES.L <- list("Common"=THREE.COMMON.GRL, 
                "Exclusive"=THREE.EXCLUSIVE.GRL,
                "Pairwise"=PAIRWISE.GRL, 
                "SampleInfo"=NAME.MAP)
  
  message("\treturning list of GRanges"); message("Completed.")
  message(paste0("\n","Samples Info"))
  print(NAME.MAP)
  return(RES.L)
}
