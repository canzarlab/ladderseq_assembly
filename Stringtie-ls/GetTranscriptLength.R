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
tLength <- biomaRt::getBM(attributes = c("transcript_length","ensembl_transcript_id","ensembl_gene_id","transcript_biotype"),
                          filters='ensembl_transcript_id', values = groundTruthFile$target_id , 
                          mart = mart)
names(tLength) <- c("TranscriptLength","TranscriptName","GeneName","transcriptBiotype")
write.table(tLength, file =args[2], col.names = T, row.names = F, quote = F)
############## filnished writing mart stuff to file

