filterBam <- function(
                      ## INPUTS
                      BAMFILE=NULL,
                      BSSPECIES=NULL,
                      EXTENTION=".Aligned.sortedByCoord.out.bam",

                      ## OPTIONS
                      SIMPLECIGAR=TRUE,
                      INCLUDE.SECONDARY.ALIGNEMNT=FALSE,
                      
                      GET.ORIGINAL.SEQUENCE=FALSE,
                      STANDARD.CONTIGS.ONLY=TRUE,
                      PERFECT.MATCH.ONLY=TRUE,
                      
                      FILTER.BY.FLAG=TRUE,
                      SELECTFLAG=c(0,16),
                      
                      USE.SIZE.FILTER=TRUE,
                      READ.SIZE.RANGE=c(18,50),
                      
                      TAGS=c("NH","NM","MD"),
                      WHAT=c("flag"),
                      
                      ## name specific
                      SPLIT.NAME.BY="-"
                      ){
  
  
  ## AUTHOR: Pavol Genzor
  ## Use: Load bam file into R
  ## 06.25.21; Version 7; refined

  ## NOTE ON BAM INFO
  ## Flag: 256 = not primary alignment; 272 = reverse strand not primary alignment; 
  ## Flag: 0 = forward unpaired unique alignment 16 = reverse unpaired unique alignment
  ## Tags: NH:i:1 = unique alignment
  ## Tags: NM = edit distance to the reference 
  ##
  ## NOTE: This program requires prepareFastq.R to be able to extract multiplicity information

  ## libraries
  suppressPackageStartupMessages({library("data.table");library("dplyr");library("Rsamtools");
    library("GenomicAlignments");library("BSgenome.Hsapiens.UCSC.hg38");
    library("BSgenome.Dmelanogaster.UCSC.dm6"); library("BSgenome.Mmusculus.UCSC.mm10")})

  ## check input
  if(is.null(BAMFILE)) stop("Please provide full path to a .bam file !!!")
  if(is.null(BSSPECIES)) stop("Please provide BSSPECIES name !!!")
  if(isTRUE(GET.ORIGINAL.SEQUENCE)){WHAT=c("flag","seq")}

  ## for report
  PROGRESS.L <- list()
  
  FILE.NAME <- gsub(EXTENTION,"",basename(BAMFILE))
  message("Processing ...")
  
  ## PARAMETERS FOR LOADING BAM FILE
  PARAM = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, 
                                                                isSecondaryAlignment = INCLUDE.SECONDARY.ALIGNEMNT),
                                  tag = TAGS, simpleCigar = SIMPLECIGAR, what = WHAT)
  message(" prepared loading parameters")
  message(paste0("\tTAGS:\t",paste(TAGS, collapse = ", ")))
  message(paste0("\tCIGAR:\t",ifelse(isTRUE(SIMPLECIGAR),"simple cigar","all cigar")))
  message(paste0("\tWHAT:\t",paste0(WHAT,collapse = ", ")))
  message(" loading .bam file into GAlignemnts")
  GA <- GenomicAlignments::readGAlignments(file = BAMFILE, use.names = TRUE, param = PARAM)
  
  ## ***
  GA.IN <- length(GA)
  PROGRESS.L[["INPUT"]] <- GA.IN
  message(paste0("\tIMPORTED: ", GA.IN))
  
  if(isTRUE(USE.SIZE.FILTER)){
    message(" filtering by read size")
    message(paste0("\tRANGE:\t",paste0(READ.SIZE.RANGE, collapse = "-")))
    GA <- GA[width(GA) %in% seq(READ.SIZE.RANGE[1],READ.SIZE.RANGE[2],by = 1)] 
    
    ## ***
    REMAINDER = (length(GA)/GA.IN)*100
    PROGRESS.L[["SIZE_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2) )) }
  
  if(isTRUE(STANDARD.CONTIGS.ONLY)){
    message(" removing non-standard contigs")
    REG.CHR <- standardChromosomes(eval(parse(text = BSSPECIES)))
    GAR <- GA[seqnames(GA) %in% REG.CHR]
    GenomeInfoDb::seqlevels(GAR) <- REG.CHR
    
    ## ***
    REMAINDER = (length(GAR)/GA.IN)*100
    PROGRESS.L[["CONTIG_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2) )) }
  else{ GAR <- GA }

  if(isTRUE(FILTER.BY.FLAG)){
    message(" selecting ONLY primary alignemnts") ## Adjust to use different flags
    GARP <- GAR[mcols(GAR)[["flag"]] %in% SELECTFLAG]
    message(" removing flag column")
    mcols(GARP)[["flag"]] <- NULL
    
    ## ***
    REMAINDER = (length(GARP)/GA.IN)*100
    PROGRESS.L[["FLAG_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2) )) }
  else { GARP <- GAR }
  
  if(isTRUE(PERFECT.MATCH.ONLY)){
    message(" removing reads with mismatches")
    GARP <- GARP[mcols(GARP)[["NM"]] %in% c(0)] 
    mcols(GARP)[["NM"]] <- NULL
    mcols(GARP)[["MD"]] <- NULL 
    
    ## ***
    REMAINDER = (length(GARP)/GA.IN)*100
    PROGRESS.L[["MISMATCH_FILTER"]] <- REMAINDER
    message(paste0("\tREMAINDER: ", round(REMAINDER, digits = 2) )) }
  
  if(isTRUE(GET.ORIGINAL.SEQUENCE)){
    message(" retrieving original read sequences")
    BAMSEQ <- mcols(GARP)[["seq"]]
    ISONMINUS <- as.logical(GenomicAlignments::strand(GARP) == "-")
    BAMSEQ[ISONMINUS] <- Biostrings::reverseComplement(BAMSEQ[ISONMINUS])
    mcols(GARP)[["seq"]] <- BAMSEQ }
  
  message(" converting to GRanges")
  GARP.GR <- GenomicRanges::granges(GARP, use.names = TRUE, use.mcols = TRUE)
  PROGRESS.L[["FINAL"]] <- length(GARP.GR)
  
  message(" adding multiplicity column [MULT]")
  mcols(GARP.GR)[["MULT"]] <- as.integer(dplyr::nth(
    data.table::tstrsplit(
      dplyr::nth(
        data.table::tstrsplit(names(GARP.GR),split = SPLIT.NAME.BY),-1), split = "M"),-1))
  
  ## RETURN
  message("Done!")
  message("")
  print(as.data.table(PROGRESS.L))
  return(GARP.GR)
}
