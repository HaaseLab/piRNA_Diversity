rmskGTF2BED <- function(## INPUT
                        RMSK.GTF = NULL,
                        
                        ## OUTPUT OPTIONS
                        RETURN.GR = TRUE,
                        KEEP.ALL.CONTIGS = FALSE,
                        
                        ## EXPORT OPTIONS
                        EXPORT.BED = FALSE,
                        
                        ## OTHER
                        MC.CORES = 4,
                        GTF.COLUMN.NAMES = c("chr","database","feature","start","end","score1","strand","score2","description")
                        ){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Load repeatmasker gtf annotation into R
  ## 07.08.21; Verions 6, refined
  
  ## LIBRARIES
  suppressPackageStartupMessages({library("data.table");library("stringr");
    library("GenomicRanges") })

  ## INPUT CHECKING
  if(is.null(RMSK.GTF)){stop("Please provide rmsk (TE toolkit) GTF file !")}
  
  message("Processing ...")
  message(" loading")
  GTF <- data.table::fread(input = RMSK.GTF, nThread = MC.CORES)
  colnames(GTF) <- GTF.COLUMN.NAMES
  
  ## EXTRACT FROM DESCRIPTION COLLUMN
  DESCRIPTION <- stringr::str_split(GTF[["description"]], pattern = "\"", simplify = TRUE)

  message(" creating bed")
  TE.BED <- data.table(chr = GTF[["chr"]],
                       start = GTF[["start"]],
                       end = GTF[["end"]],
                       name = paste(GTF[["database"]], seq(1,nrow(GTF)), sep = "_"),
                       strand = GTF[["strand"]],
                       gene_id = DESCRIPTION[,2],
                       transcript_id = DESCRIPTION[,4],
                       family_id = DESCRIPTION[,6],
                       class_id = DESCRIPTION[,8])
  
  if(isFALSE(KEEP.ALL.CONTIGS)){
    # note: not the best way
    message(" removing unused contigs")
    TE.BED <- TE.BED[grep("_",TE.BED[["chr"]],invert = TRUE),] }
  
  if(isTRUE(EXPORT.BED)){
    FILE.NAME <- data.table::tstrsplit(utils::tail(unlist(tstrsplit(GTF, split = "/")), n=1), split = "\\.")[[1]]
    utils::write.table(x = TE.BED, file = paste(FILE.NAME,".bed", sep = ""), 
                       quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    message(paste("Saved",paste(FILE.NAME,".bed", sep = ""),"file in the working directory", sep = " "))}
  
  message("Done.")
  if(isTRUE(RETURN.GR)){
    TE.GR <- makeGRangesFromDataFrame(df = TE.BED, keep.extra.columns = TRUE)
    return(TE.GR) }
  else{return(TE.BED) }
}
