miRbase2BED <- function(
                        ## INPUT
                        miRBASEFILE = NULL,
                        
                        ## OPTIONS
                        miRNAOnly = TRUE){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Import miRbase gff/gtf file into R
  ## 06.28.21; Version 3; refined

  ## LIBRARIES
  suppressPackageStartupMessages({library("data.table")})
  
  ## check for input
  if(is.null(miRBASEFILE)){stop("Please provide .gtf2 or .gff3 file from miRbase")}
  
  ## which file type
  FILE.TYPE <- tstrsplit(tail(unlist(tstrsplit(miRBASEFILE, split = "/")), n=1), split = "\\.")[[2]]
  
  ## GFF3
  
  if(FILE.TYPE %in% "gff3"){
    message("Processing .gff3 file")
    GFF <- fread(miRBASEFILE, header = FALSE)
    GFF.COLNAMES <- c("chr","score0","feature","start","end","score1","strand","score2","description")
    colnames(GFF) <- GFF.COLNAMES
    
    ## split descriptions column
    ID = tstrsplit(tstrsplit(GFF$description, split = ";")[[1]], split = "=")[[2]]
    ALIAS = tstrsplit(tstrsplit(GFF$description, split = ";")[[2]], split = "=")[[2]]
    NAME = tstrsplit(tstrsplit(GFF$description, split = ";")[[3]], split = "=")[[2]]
    ORIGIN = tstrsplit(tstrsplit(GFF$description, split = ";")[[4]], split = "=")[[2]]
    
    ## make bed
    miRBED <- data.table(chr = GFF$chr, start = GFF$start, end = GFF$end,
                         name = NAME, feature = GFF$feature, strand = GFF$strand,
                         gene_id = ID, alias = ALIAS, origin =  ORIGIN)
    
    ## only mature miRNAs
    if(isTRUE(miRNAOnly)){miRBED <- miRBED[miRBED$feature %in% "miRNA",]} }
  
  ## GTF2
  
  if(FILE.TYPE %in% "gtf2"){
    message("Processing .gtf2 file")
    GTF2 <- fread(miRBASEFILE, header = FALSE)
    GTF2.COLNAMES <- c("chr","database","feature","start","end","score1","strand","score2","description")
    colnames(GTF2) <- GTF2.COLNAMES
    
    ## split descriptions column
    GID <- tstrsplit(tstrsplit(tstrsplit(GTF2$description, split = ";")[[1]], split =" ")[[2]], split ="\"")[[2]]
    TID <- tstrsplit(tstrsplit(tstrsplit(GTF2$description, split = ";")[[2]], split =" ")[[3]], split ="\"")[[2]]
    
    ## MAKE BED TABLE
    miRBED <- data.table(chr = GTF2$chr, start = GTF2$start,end = GTF2$end,
                         gene_id = GID, transript_id = TID, strand = GTF2$strand,
                         database = GTF2$database) 
    
    ## only mature miRNAs
    if(isTRUE(miRNAOnly)){miRBED <- miRBED[miRBED$feature %in% "miRNA",]}}
  
  message("Done.")
  return(miRBED)
}
