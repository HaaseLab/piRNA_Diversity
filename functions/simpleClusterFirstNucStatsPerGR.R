simpleClusterFirstNucStatsPerGR <- function(
                                            ## INPUT
                                            GR = NULL,
                                            RANGE.NAME = NULL,
                                            CLUSTER.NAME.COLUMN = "seqnames",
                                            
                                            BSSPECIES = NULL,
                                            BSSPECIES.VERSION = "mm10",
                                            
                                            SAMPLE.SEP = "_",
                                            READ.NAMES.SPLIT = "_",
                                            ALPHABET = c("A","C","G","T"),
                                            SOURCE.DIR = NULL){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Calculate statistics surrounding the first nucleotide of piRNAs
  ## 07.07.2021; Version 2; refined

  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("plyr"); library("dplyr"); library("reshape"); library("BSgenome");
    library("GenomicRanges"); library("GenomeInfoDb"); library("XVector"); library("Biostrings")  })
  
  ## CONDITIONAL LIBRARIES
  if(BSSPECIES.VERSION %in% "mm9"){
    if(any(grepl("package:BSgenome.Mmusculus.UCSC.mm10", search()))) detach("package:BSgenome.Mmusculus.UCSC.mm10") 
    suppressWarnings(suppressMessages(library("BSgenome.Mmusculus.UCSC.mm9")))
  } else { 
    if(any(grepl("package:BSgenome.Mmusculus.UCSC.mm9", search()))) detach("package:BSgenome.Mmusculus.UCSC.mm9") 
    suppressWarnings(suppressMessages(library("BSgenome.Mmusculus.UCSC.mm10"))) }
  
  ## HOMEMADE FUNCTIONS
  source(paste0(SOURCE.DIR,"simpleNewGRFromReadNames.R"))
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GR!")
  if(is.null(BSSPECIES)) stop("Please provide BSSPECIES!")
  if(is.null(RANGE.NAME)) stop("Please provide a RANGE.NAME!")
  if(isFALSE("MULT" %in% colnames(mcols(GR)))) stop("GR must have MULT column!")
  if(is.null(SOURCE.DIR)) stop("Please provide SOURCE.DIR with functions!")

  ## transpose with names
  N.GR <- simpleNewGRFromReadNames(GR = GR, 
                                  READ.NAMES.SPLIT = READ.NAMES.SPLIT)
  
  message(" FUN: cluster counts")
  
  ## fix the levels
  seqlevels(N.GR) <- sortSeqlevels(seqlevels(N.GR))
  seqinfo(N.GR) <- keepStandardChromosomes(seqinfo(x = eval(parse(text = BSSPECIES))))
  
  ## get sequences
  message("\tgetting sequences ... slow")
  N.SEQ <- getSeq(eval(parse(text = BSSPECIES)), N.GR)
  
  ## add first name & nucleotide
  mcols(N.GR)[["seq"]] <- subseq(N.SEQ, start = 1, end = 1)
  
  ## make data.table
  N.DT <- as.data.table(N.GR)
  N.DT[["sample"]] <- RANGE.NAME
  
  ## message
  message(" FUN: calculate stats")
  
  ## total
  TOTAL.SAMPLE.READS <- sum(N.DT[["MULT"]])
  
  ## counts
  message("\t calculating counts")
  COUNTS <- N.DT[,sum(.SD), by = c("CLUSTER","seq"), .SDcols = "MULT"]
  COUNT <- dcast.data.table(data = COUNTS, CLUSTER ~ seq, value.var = "V1")
  colnames(COUNT) <- c("CLUSTER",paste(paste("count",ALPHABET, sep = "."),RANGE.NAME,sep = SAMPLE.SEP))
  COUNT[, paste("count.SUM",RANGE.NAME,sep = SAMPLE.SEP) := sum(.SD), by = "CLUSTER"]
  COUNT <- COUNT[!is.na(COUNT[[paste("count.SUM",RANGE.NAME,sep = SAMPLE.SEP)]])]
  COUNT[, paste("total.reads",RANGE.NAME,sep = SAMPLE.SEP) := TOTAL.SAMPLE.READS]
  
  # per cluster mean
  message("\t calculating means")
  P.CL.MEANS <- N.DT[,lapply(.SD, mean), by = c("CLUSTER"), .SDcols = "MULT"]
  colnames(P.CL.MEANS) <- c("CLUSTER",paste("cluster.mean",RANGE.NAME,sep = "_"))
  
  # per nucleotide mean
  P.NT.MEANS <- N.DT[,lapply(.SD, mean), by = c("CLUSTER","seq"), .SDcols = "MULT"]
  MEAN <- dcast.data.table(data = P.NT.MEANS, CLUSTER ~ seq, value.var = "MULT")
  colnames(MEAN) <- c("CLUSTER",paste(paste("mean",ALPHABET, sep = "."),RANGE.NAME,sep = SAMPLE.SEP))
  MEAN[, paste("mean.SUM",RANGE.NAME,sep = SAMPLE.SEP) := sum(.SD), by = "CLUSTER"]
  nMEAN <- setDT(join_all(dfs = list(MEAN,P.CL.MEANS),  by = "CLUSTER", type = "full"))
  nMEAN <- nMEAN[!is.na(nMEAN[[paste("mean.SUM",RANGE.NAME,sep = SAMPLE.SEP)]])]
  nMEAN[, paste("mean.AGC",RANGE.NAME,sep = SAMPLE.SEP) := rowMeans(.SD), by = "CLUSTER",.SDcols = 
          paste(paste("mean",c("A","C","G"), sep = "."),RANGE.NAME,sep = SAMPLE.SEP)]
  
  # per cluster median
  message("\t calculating medians")
  P.CL.MEDIANS <- N.DT[,lapply(.SD, median), by = c("CLUSTER"), .SDcols = "MULT"]
  colnames(P.CL.MEDIANS) <- c("CLUSTER",paste("cluster.median",RANGE.NAME,sep = "_"))
  
  # per nucleotide median
  P.NT.MEDIANS <- N.DT[,lapply(.SD, median), by = c("CLUSTER","seq"), .SDcols = "MULT"]
  MEDIAN <- dcast.data.table(data = P.NT.MEDIANS, CLUSTER ~ seq, value.var = "MULT")
  colnames(MEDIAN) <- c("CLUSTER",paste(paste("median",ALPHABET, sep = "."),RANGE.NAME,sep = SAMPLE.SEP))
  MEDIAN[, paste("median.SUM",RANGE.NAME,sep = SAMPLE.SEP) := sum(.SD), by = "CLUSTER"]
  nMEDIAN <- setDT(join_all(dfs = list(MEDIAN,P.CL.MEDIANS),  by = "CLUSTER", type = "full"))
  nMEDIAN <- nMEDIAN[!is.na(nMEDIAN[[paste("median.SUM",RANGE.NAME,sep = SAMPLE.SEP)]])]
  nMEDIAN[, paste("median.AGC",RANGE.NAME,sep = SAMPLE.SEP) := lapply(.SD,median), by = "CLUSTER",.SDcols = 
            paste(paste("median",c("A","C","G"), sep = "."),RANGE.NAME,sep = SAMPLE.SEP)]
  
  # join stat tables
  N1J <- setDT(join_all(dfs = list(COUNT,nMEAN,nMEDIAN), by = "CLUSTER", type = "full"))
  
  ## return
  return(N1J)
}





