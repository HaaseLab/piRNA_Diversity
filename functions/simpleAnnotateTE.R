simpleAnnotateTE <- function(## INPUT
                            GR = NULL,
                            RMSK.GR = NULL,

                            ## SETTINGS
                            USE.CLASS = TRUE,
                            USE.FAMILY = TRUE,
                            ORDER.BY = "LENGTH",
                            OVERLAP.TYPE = "within",
                                
                            ## RETURN
                            RETURN.ALL = FALSE,
                            SHOW.PROGRESS = FALSE){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Generate annotation tables using RMSK
  ## 07.08.21; Version 3; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("GenomicRanges")})
  
  ## ERROR
  if(is.null(GR)) stop("Please provide GR to be annotated.")
  if(is.null(RMSK.GR)) stop("Provide RMSK.GR with family_id and class_id mcols().")
  
  message("Starting...\n")
  
  ## PREPARE RMSK
  
  message("\tPreparing RMSK")
  
  ## Split into GRL
  RMSK.GR.CLASS.L <- split(RMSK.GR, ~class_id)
  RMSK.GR.FAMILY.L <- split(RMSK.GR, ~family_id)
  
  ## Create element table and order it
  CLASS.INFO <- data.table("CLASS" = names(RMSK.GR.CLASS.L),
                           "ELEMEN_COUNT" = unlist(lapply(RMSK.GR.CLASS.L, length)),
                           "ELEMENT_BP_LENGTH" = unlist(lapply(RMSK.GR.CLASS.L, function(i){ sum(width(i)) })) )
  setorderv(x = CLASS.INFO, cols = grep(ORDER.BY, names(CLASS.INFO), value = TRUE),order = -1)
  FAMILY.INFO <- data.table("FAMILY" = names(RMSK.GR.FAMILY.L),
                            "ELEMEN_COUNT" = unlist(lapply(RMSK.GR.FAMILY.L, length)),
                            "ELEMENT_BP_LENGTH" = unlist(lapply(RMSK.GR.FAMILY.L, function(i){ sum(width(i)) })) )
  setorderv(x = FAMILY.INFO, cols = grep(ORDER.BY, names(CLASS.INFO), value = TRUE),order = -1)
  
  ## Set order 
  CLASS.ORDER = CLASS.INFO[["CLASS"]]
  FAMILY.ORDER = FAMILY.INFO[["FAMILY"]]
  
  ## ANNOTATE
  
  ## CLASS
  if(isTRUE(USE.CLASS)){
    
    message("\tAnnotating by CLASS")
    INPUT.GR = GR
    IN.CLASS.GRL.S = IN.CLASS.GRL.AS = GRangesList()
    
    ## Loop through classes
    for(c in CLASS.ORDER){
      if(isTRUE(SHOW.PROGRESS)) {message("\t class: ",c)}
      IN.CLASS.GRL.S[[c]] = subsetByOverlaps(x = INPUT.GR, ranges = RMSK.GR.CLASS.L[[c]], type = OVERLAP.TYPE)
      INPUT.GR = subsetByOverlaps(x = INPUT.GR, ranges = IN.CLASS.GRL.S[[c]], type = OVERLAP.TYPE, invert = TRUE)
      IN.CLASS.GRL.AS[[c]] = subsetByOverlaps(x = INPUT.GR, ranges = invertStrand(RMSK.GR.CLASS.L[[c]]), type = OVERLAP.TYPE)
      INPUT.GR = subsetByOverlaps(x = INPUT.GR, ranges = IN.CLASS.GRL.AS[[c]], type = OVERLAP.TYPE,invert = TRUE) }
    IN.CLASS.GRL.S[["OTHER"]] <- INPUT.GR
    
    ## Class tables
    CLASS.S.DT <- as.data.table(unlist(lapply(IN.CLASS.GRL.S,length)), keep.rownames = TRUE)
    colnames(CLASS.S.DT) <- c("CLASS","S")
    CLASS.AS.DT <- as.data.table(unlist(lapply(IN.CLASS.GRL.AS,length)), keep.rownames = TRUE)
    colnames(CLASS.AS.DT) <- c("CLASS","AS")
    CLASS.DT <- setDT(join_all(dfs = list(CLASS.S.DT, CLASS.AS.DT), by = "CLASS", type = "full",))
    for (j in names(CLASS.DT))set(CLASS.DT,which(is.na(CLASS.DT[[j]])),j,0)
    CLASS.DTm <- setDT(melt(data = CLASS.DT, id.vars = "CLASS", value.name = "COUNT",variable.name = "STRAND"))
    CLASS.DTm[["FREQ"]] = CLASS.DTm[["COUNT"]] / sum(CLASS.DTm[["COUNT"]])
    CLASS.DTm[["CLASS"]] = factor(x = CLASS.DTm[["CLASS"]], levels = c(CLASS.ORDER,"OTHER"))
    
  } else {IN.CLASS.GRL.S = IN.CLASS.GRL.AS = CLASS.DTm = "Not executed!"}
  
  ## FAMILY
  if(isTRUE(USE.FAMILY)){
    
    message("\tAnnotating by FAMILY")
    INPUT.GR = GR
    IN.FAMILY.GRL.S = IN.FAMILY.GRL.AS = GRangesList()
    
    ## Loop through families
    for(c in FAMILY.ORDER){
      if(isTRUE(SHOW.PROGRESS)) {message("\t family: ",c)}
      IN.FAMILY.GRL.S[[c]] = subsetByOverlaps(x = INPUT.GR, ranges = RMSK.GR.FAMILY.L[[c]], type = OVERLAP.TYPE)
      INPUT.GR = subsetByOverlaps(x = INPUT.GR, ranges = IN.FAMILY.GRL.S[[c]], type = OVERLAP.TYPE, invert = TRUE)
      IN.FAMILY.GRL.AS[[c]] = subsetByOverlaps(x = INPUT.GR, ranges = invertStrand(RMSK.GR.FAMILY.L[[c]]), type = OVERLAP.TYPE)
      INPUT.GR = subsetByOverlaps(x = INPUT.GR, ranges = IN.FAMILY.GRL.AS[[c]], type = OVERLAP.TYPE,invert = TRUE)}
    IN.FAMILY.GRL.S[["OTHER"]] <- INPUT.GR
    
    ## Family tables
    FAMILY.S.DT <- as.data.table(unlist(lapply(IN.FAMILY.GRL.S,length)), keep.rownames = TRUE)
    colnames(FAMILY.S.DT) <- c("FAMILY","S")
    FAMILY.AS.DT <- as.data.table(unlist(lapply(IN.FAMILY.GRL.AS,length)), keep.rownames = TRUE)
    colnames(FAMILY.AS.DT) <- c("FAMILY","AS")
    FAMILY.DT <- setDT(join_all(dfs = list(FAMILY.S.DT, FAMILY.AS.DT), by = "FAMILY", type = "full",))
    for (j in names(FAMILY.DT))set(FAMILY.DT,which(is.na(FAMILY.DT[[j]])),j,0)
    FAMILY.DTm <- setDT(melt(data = FAMILY.DT, id.vars = "FAMILY", value.name = "COUNT",variable.name = "STRAND"))
    FAMILY.DTm[["FREQ"]] = FAMILY.DTm[["COUNT"]] / sum(FAMILY.DTm[["COUNT"]])
    FAMILY.DTm[["FAMILY"]] = factor(x = FAMILY.DTm[["FAMILY"]], levels = c(FAMILY.ORDER,"OTHER"))
    
  } else {IN.FAMILY.GRL.S = IN.FAMILY.GRL.AS = FAMILY.DTm = "Not executed"}
  
  ## COMPILE RESULTS
  
  RES.L <- list("GRLS" = list("CLASS" = list("CLASS.S" = IN.CLASS.GRL.S, "CLASS.AS" = IN.CLASS.GRL.AS),
                              "FAMILY" = list("FAMILY.S" = IN.FAMILY.GRL.S), "FAMILY.AS" = IN.FAMILY.GRL.AS),
                "RMSK.DATA" = list("CLASS" = CLASS.INFO, "FAMILY" = FAMILY.INFO),
                "RESULTS" = list("CLASS" = CLASS.DTm, "FAMILY" = FAMILY.DTm))
  
  ## RETURN
  
  message("Done.")
  if(isTRUE(RETURN.ALL)) {return(RES.L)} 
  else {return(RES.L[["RESULTS"]])}
  
}
