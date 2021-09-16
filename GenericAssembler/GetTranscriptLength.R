### Gets length of transcripts from BioMart ###

### Author: Shounak Chakraborty###
### 20.08.2020 ###



require(biomaRt)


## Command line arguments ##

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

## human
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'www.ensembl.org')


## writing down the mart stuff to file so that its faster to access!!
groundTruthFile <- read.table(args[1],header = T)
transcriptNames <- lapply(groundTruthFile$target_id, function(x) strsplit(as.character(x),'\\.')[[1]][1])
transcriptNames <- as.data.frame(unlist(transcriptNames))
names(transcriptNames) <- c("transcriptNames")
groundTruthFile$transcriptNames <- transcriptNames$transcriptNames
tLength <- biomaRt::getBM(attributes = c("transcript_length","ensembl_transcript_id","ensembl_gene_id","transcript_biotype"),
                          filters='ensembl_transcript_id', values = groundTruthFile$transcriptNames ,
                          mart = mart)
names(tLength) <- c("TranscriptLength","transcriptNames","GeneName","transcriptBiotype")
mergedTranscripts <- merge(groundTruthFile,tLength, by= "transcriptNames")
mergedTranscripts <- mergedTranscripts[,c(2,7,8,9)]
names(mergedTranscripts) <- c("TranscriptName","TranscriptLength","GeneName","transcriptBiotype")
write.table(mergedTranscripts, file =args[2], col.names = T, row.names = F, quote = F)
############## filnished writing mart stuff to file
