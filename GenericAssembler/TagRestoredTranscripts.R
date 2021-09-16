# libraries
require(dplyr)



# #### Hadcoded arguments ####
# gtfFileName <- "/home/schakraborty/Documents/HomeOffice/LadderSeq/refBasedAnalysis/30Mill/sim_1/restoringAnalysis/test/singleExonsRemovedInitialAssembly.restored.combined.gtf"
# outFileName <- ""




#### Command line arguments
args = commandArgs(trailingOnly=TRUE)
## gtf
gtfFileName <- args[1]
## output file name
outFileName <- args[2]




tagGtf <- function(gtfFileName,outFileName){
  gtfFile <- read.delim(gtfFileName, header = FALSE, comment.char="#",quote = "")
  
  transcriptList <- lapply(gtfFile$V9 , function(x) strsplit(as.character(x),'transcript_id')[[1]][2])
  transcriptList<- as.data.frame(unlist(transcriptList))
  transcriptList <- lapply(transcriptList$`unlist(transcriptList)` , function(x) strsplit(as.character(x),';')[[1]][1])
  transcriptList<- as.data.frame(unlist(transcriptList))
  gtfFile$TranscriptId <- transcriptList$`unlist(transcriptList)`
  
  
  oldTranscriptIdList <- lapply(gtfFile$V9 , function(x) strsplit(as.character(x),'oId')[[1]][2])
  oldTranscriptIdList<- as.data.frame(unlist(oldTranscriptIdList))
  oldTranscriptIdList <- lapply(oldTranscriptIdList$`unlist(oldTranscriptIdList)` , function(x) strsplit(as.character(x),';')[[1]][1])
  oldTranscriptIdList<- as.data.frame(unlist(oldTranscriptIdList))
  gtfFile$oId <- oldTranscriptIdList$`unlist(oldTranscriptIdList)`
  
  
 # restoredTranscriptIds <- gtfFile %>% subset(grep("restored",NA)) 
  restoredTranscriptIds <- gtfFile[grep("restored", gtfFile$oId), ]
  restoredTranscriptIds <- as.data.frame(restoredTranscriptIds[,10])
  names(restoredTranscriptIds) <- c("TranscriptId")
  
  
  unrestoredTranscripts <- gtfFile %>% subset(!(TranscriptId %in% restoredTranscriptIds$TranscriptId))
  restoredTranscripts <- gtfFile %>% subset(TranscriptId %in% restoredTranscriptIds$TranscriptId)
  
  
  restoredTranscriptsV9 <- lapply(restoredTranscripts$V9 , function(x) gsub("transcript_id \"TCONS_", "transcript_id \"TCONS_restored_", x))
  restoredTranscriptsV9 <- as.data.frame(unlist(restoredTranscriptsV9))
  restoredTranscripts$V9 <- restoredTranscriptsV9$`unlist(restoredTranscriptsV9)`
  
  
  restoredTranscripts <- restoredTranscripts[,c(-10,-11)]
  unrestoredTranscripts <- unrestoredTranscripts[,c(-10,-11)]
  
  ## Writing the tagged gtf
  taggedGtf <- rbind(unrestoredTranscripts,restoredTranscripts)
  write.table(taggedGtf,file = outFileName, sep = "\t",row.names = F, col.names = F, quote = F)
  
}


tagGtf(gtfFileName,outFileName)
