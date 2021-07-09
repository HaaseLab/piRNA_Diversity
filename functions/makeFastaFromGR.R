makeFastaFromGR <- function(
                            ## INPUT
                            GRL=NULL,
                            BSSPECIES=NULL,
                            BSSPECIES.VERSION="mm10",
                            
                            ## OPTIONS
                            CHECK.MAPPING=TRUE, 
                            
                            ## OUTPUT
                            SAVE.TO.FILE=TRUE,
                            FA.DIR=NULL,
                            MC.CORES = 3){
  
  ## Pavol Genzor
  ## Use: Make fasta from Genomic Range
  ## 06.30.21; Version 4; refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({library("data.table");library("ShortRead");library("Biostrings");
    library("GenomicRanges");library("parallel");library("BSgenome");
    library("BSgenome.Dmelanogaster.UCSC.dm6");library("BSgenome.Hsapiens.UCSC.hg38")})

  if(BSSPECIES.VERSION %in% "mm9"){
    suppressWarnings(suppressMessages(library("BSgenome.Mmusculus.UCSC.mm9")))
  } else { suppressWarnings(suppressMessages(library("BSgenome.Mmusculus.UCSC.mm10"))) }

  ## INPUT CHECKING
  if(is.null(GRL)) stop("Please provide a GR or GRL !")
  if(is.null(BSSPECIES)) stop("Please provide BSSPECIES name !")
    if(isFALSE(is.list(GRL))) {
    message("NOTE: single GR. Making a list ...")
    GRL <- list("singleGR"=GRL)}
  if(is.null(names(GRL))) stop("The GRL list has to be named !")
  
  if(isTRUE(SAVE.TO.FILE)){
    if(is.null(FA.DIR)) stop("You need to set FA.DIR where to save a file !")
    else{ message("Processing ...") }
  } else { message("Processing without saving ...") }
  
  ##
  ## PROCESSING
  ##
  
  FA.L <- mclapply(names(GRL), mc.cores = MC.CORES, function(i){
    
    message(paste0("Processing: ",i))
    GR <- GRL[[i]]
    
    message("  trimming out-of-bound")
    seqinfo(GR) <- keepStandardChromosomes(seqinfo(eval(parse(text = BSSPECIES))))
    GR.TRIMMED <- trim(GR)
    GR <- GR[which(GR == GR.TRIMMED)]
    
    ## MAPPING
    NH <- unique(mcols(GR)[["NH"]])
    SEC.NAME <- ifelse(length(NH) == 1,"unique","all")
    message(paste0("  sample mapping: ",SEC.NAME))
    
    message("  getting sequences ...")
    SEQ <- getSeq(eval(parse(text = BSSPECIES)),GR)
    GR.DT <- as.data.table(GR)
    
    message("  preparing IDs")
    ID.MCOLS <- paste0(paste0("N",GR.DT[["N"]]),
                       paste0("NH",GR.DT[["NH"]]),
                       paste0("M",GR.DT[["MULT"]]))
    
    ID.GR <- paste(GR.DT[["seqnames"]],GR.DT[["start"]],GR.DT[["end"]],
                   GR.DT[["strand"]],ID.MCOLS,sep="_")
    
    message("  adding sequence names")
    names(SEQ) <- ID.GR
    
    message("  writting to file")
    writeFasta(object = SEQ, file = paste0(FA.DIR,i,"_",SEC.NAME,".fa"))
    return(SEQ) })
  names(FA.L) <- names(GRL)
  
  ## RETURN
  message("Done.")
  if(length(GRL) == 1){return(unlist(FA.L))}
  else { return(FA.L)}
}
