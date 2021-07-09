annotateGRanked <- function(
                              ## INPUTS
                              GR=NULL,
                              
                              ## ANNOTATIONS
                              CATEGORY.1.GR=NULL,
                              CATEGORY.2.GR=NULL,
                              CATEGORY.3.GR=NULL,
                              CATEGORY.4.GR=NULL,
                              CATEGORY.5.GR=NULL,
                              CATEGORY.NAMES=NULL,
                              
                              ## OPTIONS
                              SAMPLE.NAME=NULL,
                              USE.READS=TRUE,
                              
                              ## FILTER
                              NH.TAG = NULL, SIZE.RANGE = c(18,50),
                              
                              ## SETTINGS
                              REPORT.METRIC="FRACTION",
                              OVERLAP.TYPE="within",
                              SOURCE.DIR = NULL,
                              RETURN.ALL=FALSE){
  
  ## Pavol Genzor
  ## Use: Annotate GR with multiple categories in order
  ## 07.08.21; Version 5; Refined

  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("GenomicRanges")})
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide the GR to annotate !")
  CATS <- c("CAT1" = CATEGORY.1.GR,"CAT2" = CATEGORY.2.GR,
            "CAT3" = CATEGORY.3.GR,"CAT4" = CATEGORY.4.GR,"CAT5" = CATEGORY.5.GR)
  if(is.null(unlist(CATS))) stop("Please provide at least one annotation category in order !")
  if(is.null(CATEGORY.NAMES)) stop("Please provide vector of ordered (!) CATEGORY.NAMES !")
  if(is.null(SOURCE.DIR)) stop("Please provide SOURCE.DIR with filter function !")
  
  ## FUNCTIONS
  source(paste0(SOURCE.DIR,"simpleGRFilter.R"))
  
  ## PROGRAM
  message("Annotating ...")
  message(" order:\n\tCAT1 > CAT2 > CAT3 > CAT4 > CAT5 > OTHER")
  
  ## FILTER
  
  FINFO <- paste0("NH",NH.TAG,"SR",min(SIZE.RANGE),max(SIZE.RANGE))
  GR <- simpleGRFilter(GR = GR, RANGE.NAME = SAMPLE.NAME,
                       NH.TAG = NH.TAG, 
                       SIZE.RANGE = SIZE.RANGE)
  message("")
  
  ## OVERLAPPING
  
  ## INITIATE LIST
  FINAL.GRL <- GRangesList()
  
  message(" calculating overlaps")
  
  if(!is.null(CATEGORY.1.GR)){
    message("\twith Category 1")
    IN.CAT1.S <- subsetByOverlaps(x = GR, ranges = CATEGORY.1.GR, type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT1.S, type = OVERLAP.TYPE, invert = TRUE)
    IN.CAT1.AS <- subsetByOverlaps(x = GR, ranges = invertStrand(CATEGORY.1.GR), type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT1.AS, type = OVERLAP.TYPE, invert = TRUE)
    
    FINAL.GRL[["CAT1_S"]] <- IN.CAT1.S
    FINAL.GRL[["CAT1_AS"]] <- IN.CAT1.AS 
  } else {
    IN.CAT1.S <- GRanges()
    IN.CAT1.AS <- GRanges() } 
  
  if(!is.null(CATEGORY.2.GR)){
    message("\twith Category 2")
    IN.CAT2.S <- subsetByOverlaps(x = GR, ranges = CATEGORY.2.GR, type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT2.S, type = OVERLAP.TYPE, invert = TRUE)
    IN.CAT2.AS <- subsetByOverlaps(x = GR, ranges = invertStrand(CATEGORY.2.GR), type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT2.AS, type = OVERLAP.TYPE, invert = TRUE)
    
    FINAL.GRL[["CAT2_S"]] <- IN.CAT2.S
    FINAL.GRL[["CAT2_AS"]] <- IN.CAT2.AS 
  } else {
    IN.CAT2.S <- GRanges()
    IN.CAT2.AS <- GRanges() } 
  
  if(!is.null(CATEGORY.3.GR)){
    message("\twith Category 3")
    IN.CAT3.S <- subsetByOverlaps(x = GR, ranges = CATEGORY.3.GR, type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT3.S, type = OVERLAP.TYPE, invert = TRUE)
    IN.CAT3.AS <- subsetByOverlaps(x = GR, ranges = invertStrand(CATEGORY.3.GR), type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT3.AS, type = OVERLAP.TYPE, invert = TRUE)
    
    FINAL.GRL[["CAT3_S"]] <- IN.CAT3.S
    FINAL.GRL[["CAT3_AS"]] <- IN.CAT3.AS 
  } else {
    IN.CAT3.S <- GRanges()
    IN.CAT3.AS <- GRanges() } 
  
  if(!is.null(CATEGORY.4.GR)){
    message("\twith Category 4")
    IN.CAT4.S <- subsetByOverlaps(x = GR, ranges = CATEGORY.4.GR, type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT4.S, type = OVERLAP.TYPE, invert = TRUE)
    IN.CAT4.AS <- subsetByOverlaps(x = GR, ranges = invertStrand(CATEGORY.4.GR), type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT4.AS, type = OVERLAP.TYPE, invert = TRUE)
    
    FINAL.GRL[["CAT4_S"]] <- IN.CAT4.S
    FINAL.GRL[["CAT4_AS"]] <- IN.CAT4.AS 
  } else {
    IN.CAT4.S <- GRanges()
    IN.CAT4.AS <- GRanges() } 
  
  if(!is.null(CATEGORY.5.GR)){
    message("\twith Category 5")
    IN.CAT5.S <- subsetByOverlaps(x = GR, ranges = CATEGORY.5.GR, type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT5.S, type = OVERLAP.TYPE, invert = TRUE)
    IN.CAT5.AS <- subsetByOverlaps(x = GR, ranges = invertStrand(CATEGORY.5.GR), type = OVERLAP.TYPE)
    GR <- subsetByOverlaps(x = GR, ranges = IN.CAT5.AS, type = OVERLAP.TYPE, invert = TRUE)
    
    FINAL.GRL[["CAT5_S"]] <- IN.CAT5.S
    FINAL.GRL[["CAT5_AS"]] <- IN.CAT5.AS 
  } else {
    IN.CAT5.S <- GRanges()
    IN.CAT5.AS <- GRanges() } 
  
  message("\tremainder is in OTHER")
  GR <- subsetByOverlaps(x = GR, ranges = c(unlist(FINAL.GRL)), type = OVERLAP.TYPE, invert = TRUE)
  FINAL.GRL[["OTHER"]] <- GR
  
  ## SEQUENCE OR READ TABLE
  
  if(isTRUE(USE.READS)){
    message(" extracting READS info")
    DT <- as.data.table(unlist(lapply(FINAL.GRL, function(i){
      sum(mcols(i)[["MULT"]])})), keep.rownames = TRUE)
    DT[["dataType"]] = "reads"}
  else{
    message(" extracting SEQUENCE info")
    DT <- as.data.table(unlist(lapply(FINAL.GRL, function(i){
      length(i)})), keep.rownames = TRUE)
    DT[["dataType"]] = "sequences"}
  
  ## PROCESSING THE RESULS
  
  colnames(DT) <- c("ANNOTATION", paste0(SAMPLE.NAME,"__COUNT"),"dataType")
  DT[["STRAND"]] <- nth(tstrsplit(DT[["ANNOTATION"]],split="_"),2)
  COUNT.COL.NAME <- grep("COUNT",colnames(DT),value = TRUE)
  DT[["FRACTION"]] <- lapply(DT[[COUNT.COL.NAME]],function(i){( i/sum(DT[[COUNT.COL.NAME]]))*100 })
  colnames(DT) <- gsub("FRACTION",paste0(SAMPLE.NAME,"__FRACTION"), colnames(DT))
  
  ##  MELT AND PREPARE DATA

  DTM <- suppressWarnings(melt.data.table(DT, id.vars = c("ANNOTATION","STRAND","dataType"), value.name = "VALUE"))
  DTM[["CATS"]] <-  nth(tstrsplit(DTM[["ANNOTATION"]],split="_"),1)
  DTM[["MEASURE"]] <- nth(tstrsplit(DTM[["variable"]],split="__"),2)
  DTM[["UPDOWN"]] <- ifelse(DTM[["STRAND"]] %in% "S","+","-")
  DTM[["UD.VALUE"]] <- as.double(paste0(DTM[["UPDOWN"]],DTM[["VALUE"]]))
  DTM[["SAMPLE"]] <- nth(tstrsplit(DTM[["variable"]],split="__"),1)
  
  ## ADD COLORS
  COL.DT <- data.table(CATS = c("CAT1","CAT2","CAT3","CAT4","CAT5","OTHER"),
                       COLS = c("brown2","skyblue2","goldenrod2","gold","chartreuse4","grey30"))
  for(c in DTM[["CATS"]]) set(x = setDT(DTM), 
                              i = which(DTM[["CATS"]] %in% c), 
                              j = "COLS", 
                              value = COL.DT[CATS %in% c][["COLS"]])
  
  ## ADD NAMES
  CATEGORIES = unique(DTM[["CATS"]])[!unique(DTM[["CATS"]]) %in% "OTHER"]
  for(n in 1:length(CATEGORIES)) set(x = DTM,
                                     i = which(DTM[["CATS"]] %in% CATEGORIES[n]), 
                                     j = "NAME", 
                                     value = CATEGORY.NAMES[n])
  for(i in 1:nrow(DTM)) set(x = DTM, i = which(is.na(DTM[["NAME"]])), j = "NAME","OTHER")
  
  ## SUBSET BY REPORT.METRIC
  DTM <- DTM[MEASURE %in% REPORT.METRIC]
  DTM[["FILTERS"]] <- FINFO
  
  ## SPLIT INTO TWO TABLES
  DTM.UP <- DTM[UPDOWN %in% "+"]
  DTM.UP[["CATS"]] <- factor(DTM.UP[["CATS"]], levels = rev(DTM.UP[["CATS"]]))
  DTM.DOWN <- DTM[UPDOWN %in% "-"]
  DTM.DOWN[["CATS"]] <- factor(DTM.DOWN[["CATS"]], levels = rev(DTM.DOWN[["CATS"]]))
  
  ## RETURN
  RES.L <- list("DTM.UP" = DTM.UP, "DTM.DOWN" = DTM.DOWN, "DTM" = DTM, "COLORS" = COL.DT)
  if(isTRUE(RETURN.ALL)){
    return(RES.L)
  } else { return(RES.L[c("DTM.UP","DTM.DOWN")])}
  
}
