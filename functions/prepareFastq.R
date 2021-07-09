prepareFastq <- function(
                        ## IO
                        FASTQ.FILE=NULL,
                        OUTPUT.DIR=NULL,
                        
                        ## FUNCTIONS (UNUSED)
                        #SOURCE.DIR = NULL
                        
                        ## OPTION: PCR DUPLICATE REMOVAL
                        REMOVE.PCR.DUPLICATES=FALSE,
                        
                        ## OPTION: UMI REMOVAL
                        REMOVE.UMI.N=FALSE,
                        FIVE.PRIME.N.NUMBER=8,
                        THREE.PRIME.N.NUMBER=2,
                        
                        ## OPTION: SIZE FILTERING
                        FILTER.BY.SIZE=FALSE,
                        SIZE.RANGE=NULL,
                        
                        ## OUTPUT
                        RETURN.TABLE=FALSE){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Preprocess fastq file with UMIs for alignment
  ## 06.24.2021; Version 3; refined
  
  ## load libraries
  suppressPackageStartupMessages({library(data.table);library(plyr);library(dbplyr);library(ShortRead)})
  
  ## check inputs
  if(is.null(FASTQ.FILE)) stop("Please provide fastq file !")
  if(isTRUE(FILTER.BY.SIZE)){if(is.null(SIZE.RANGE)) stop("Please provide read size range !")}
  
  message("Processing ...")
  
  ## extract file names and set output.dir
  FILE.NAME <- dplyr::nth(data.table::tstrsplit(basename(FASTQ.FILE),split = "\\."),1)
  if(!is.null(OUTPUT.DIR)){OUT.DIR <- OUTPUT.DIR} else { OUT.DIR <- dirname(FASTQ.FILE) }
  
  message("\treading in the .fastq file")
  message(paste0("  ",FILE.NAME))
  FASTQ.F <- ShortRead::readFastq(dirPath = FASTQ.FILE)
  
  message("\textracting sequence information")
  FQSEQ <- ShortRead::sread(FASTQ.F)
  LENGTH.FQSEQ <- length(FQSEQ)
  
  if(isTRUE(REMOVE.PCR.DUPLICATES)){
    message("\tremoving PCR duplicates")
    NOPCRDUPS <- BiocGenerics::unique(FQSEQ)
    
    # for stat report
    LENGTH.NOPCRDUPS = length(NOPCRDUPS)
    PCNTDUPS = 100 - (length(NOPCRDUPS)/length(FQSEQ))*100
    message(paste0("\t\tduplication percentage: ",round(PCNTDUPS, digits = 3)))
  } else {
    message("\tNOT removing PCR duplicates")
    NOPCRDUPS <- FQSEQ; LENGTH.NOPCRDUPS = NA; PCNTDUPS = NA}
  
  if(isTRUE(REMOVE.UMI.N)){
    message("\tremowing UMI sequences")
    message(paste0("  5`N: ",FIVE.PRIME.N.NUMBER," 3`N: ",THREE.PRIME.N.NUMBER))
    NOUMI <- XVector::subseq(x = NOPCRDUPS,
                             start = FIVE.PRIME.N.NUMBER+1,
                             end = -(THREE.PRIME.N.NUMBER+1))
    
    # for stat report
    LENGTH.NOUMI = length(NOUMI)
  } else {
    message("\tNOT removing UMI sequences")
    NOUMI <- NOPCRDUP; LENGTH.NOUMI = NA}
  
  if(isTRUE(FILTER.BY.SIZE)){
    message("\tremoving reads outside of SIZE.RANGE")
    message(paste0(" ",paste(seq(SIZE.RANGE[1],SIZE.RANGE[2]),collapse = ", ")))
    NOUMI.ALL <- NOUMI[width(NOUMI) %in% seq(SIZE.RANGE[1],SIZE.RANGE[2])]
    NOUMI.ALL <- ShortRead::srsort(NOUMI.ALL)
    
    # for stat report
    LENGTH.NOUMI.ALL <- length(NOUMI.ALL)
  } else {
    message("\tNOT removing reads outside of SIZE.RANGE")
    NOUMI.ALL <- NOUMI; NOUMI.ALL <- ShortRead::srsort(NOUMI.ALL)
    LENGTH.NOUMI.ALL = NA}
  
  message("\tcreating a data.table of all reads")
  DTT <- as.data.table(NOUMI.ALL)
  DTT[["rname"]] <- paste0(FILE.NAME,"-","R",seq(1,nrow(DTT)))
  colnames(DTT) <- gsub("x","seq",colnames(DTT))
  
  # for stat report
  LENGTH.DTT <- nrow(DTT)
  
  message("\tcollapsing into unique sequences")
  DT <- as.data.table(NOUMI.ALL)[,.N, by=c("x")]
  DT[["sname"]] <- paste0("S",seq(1,nrow(DT)),"M",DT[["N"]])
  DT[["lsname"]] <- paste0(FILE.NAME,"-","S",seq(1,nrow(DT)),"M",DT[["N"]])
  colnames(DT) <- gsub("x","seq",colnames(DT))
  
  # for stat report
  LENGTH.DT <- nrow(DT)
  
  message("\tjoining ALL and UNIQUE")
  LDT <- setDT(plyr::join_all(dfs = list(DTT,DT), by = "seq", type = "full"))
  LDT[["rsname"]] <- paste0(LDT[["rname"]],LDT[["sname"]])
  
  message("\tcreating new DNAStringSet objects")
  UNIQ.SEQ <- DNAStringSet(DT[["seq"]])
  names(UNIQ.SEQ) <- DT[["lsname"]]
  ALL.SEQ <- DNAStringSet(LDT[["seq"]])
  names(ALL.SEQ) <- LDT[["rsname"]]
  
  # the stat report
  message("\tcreating stat report")
  LIBSTAT <- data.table(metric = c("initialFastqReads","withoutPCRDuplicates","withoutUMI",
                                   "selectedSizeRange", "allReads", "uniqueSequences"),
                        values = c(LENGTH.FQSEQ, LENGTH.NOPCRDUPS, LENGTH.NOUMI,
                                   LENGTH.NOUMI.ALL, LENGTH.DTT, LENGTH.DT) )
  
  message("\twritting .fasta files and saving reports")
  OUT.DIR.NAMES <- c("totalReads/","uniqueSequences/","libraryStats/")
  OUT.DIR.PATHS <- paste(OUT.DIR,OUT.DIR.NAMES, sep = "/")
  for(i in OUT.DIR.PATHS){ifelse(dir.exists(path = i),TRUE,dir.create(path = i))}
  
  message(paste0(" output directory: ",OUT.DIR))
  writeFasta(object = ALL.SEQ, file = paste0(OUT.DIR,"/totalReads/",FILE.NAME,".ALLREADS.fa"))
  writeFasta(object = UNIQ.SEQ, file = paste0(OUT.DIR,"/uniqueSequences/",FILE.NAME,".UNIQSEQS.fa"))
  write.csv(x = LIBSTAT, file = paste0(OUT.DIR,"/libraryStats/",FILE.NAME,".stats.txt"))
  
  if(isTRUE(RETURN.TABLE)){
    message(" returning table");message("Done."); message("")
    return(LDT)
  } else {message("Done.");message("") }
  
}
