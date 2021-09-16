# Reads a set of transcripts assembled (transcript_names containing the band where they were assembled) using Ladder Seq reads.
# Filters out the transcripts which are assembled in a particular band but whose length does not agree with the length of the band
# Inputs:
# 1.



# Libraries
require(dplyr)

# Hardcoded Arguments
#initialGffFile <- "/algbio1/shounak/raw/simulated/human/RefAssemblyAnalysis/stringtie_out_individual_noFilter.gtf"
#resolution <- "HIGH"
#filterType <- "CORRECT_PLUS_ONE"
#outputPath <- "/home/schakraborty/Documents/HomeOffice/LadderSeq/refBasedAnalysis/30Mill/sim_1/assemblies/normal/band6/test"
#bandProbabilitiesFile <- "/home/schakraborty/Documents/HomeOffice/LadderSeq/kallistoAnalysis/GoodExamples/estimatedBeta_high_text.txt"
#threshold <- 2
## INDIVIDUAL or COMBINED
#bandType <- "INDIVIDUAL"
## Poly A tail length
#polyA <- 200
#bandName <- "comb23"


#### Command line arguments
args = commandArgs(trailingOnly=TRUE)
## annotation gtf
initialGffFile <- args[1]
## Resolution
resolution <- args[2]
## Filter Type : EQUAL, LARGER, SHORTER
filterType <- args[3]
## OutputPath
outputPath <- args[4]
## High resoulution band probabilities file
bandProbabilitiesFile <- args[5]
## Threshold
threshold <- args[6]
# INDIVIDUAL or COMBINED
bandType <- args[7]
# Poly A tail length
polyA <- args[8]
# bandName
bandName <- args[9]




getBand <- function(transcriptLength,resolution,bandProbabilitiesFile,threshold){
  if(resolution=="LOW"){

    band_list <- c(0,1000,1500,2000,3000,4000,6000)
    correct_band <- findInterval(transcriptLength,band_list)
    bands_to_include <- correct_band

  }else if(resolution=="HIGH"){
    #transcriptLength <-
    band_list_high <- c(0,910,1083,1268,1440,1590,1757,1920,2070,2273,2462,2701,2932,3230,3554,3893,4525,5324,6378)
    band_probabilities <- read.delim(bandProbabilitiesFile, header = F, comment.char="#",quote = "")
    band_probabilities <- band_probabilities[, colSums(is.na(band_probabilities)) == 0]
    correct_sub_band <- findInterval(transcriptLength,band_list_high)
    correctSubBandProb <- band_probabilities[correct_sub_band,]

    n <- length(correctSubBandProb)
    max_prob <- sort(correctSubBandProb,partial=n-1)[1,n]
    max_pos <- match(max_prob,correctSubBandProb)

    second_max_prob <- sort(correctSubBandProb,partial=n-1)[1,n-1]
    second_max_pos <- match(second_max_prob,correctSubBandProb)

    if(((max_prob/second_max_prob) >= 0.5) && ((max_prob/second_max_prob) <= as.numeric(threshold))){
      bands_to_include <- c(max_pos,second_max_pos)
    }else{
      bands_to_include <- c(max_pos)
    }
  }

  return(bands_to_include)
}

## High resolution filtering is always correct and correct plus one
filterHighRes <- function(assembledBand,correctBand,bandType){
  if(bandType == "INDIVIDUAL"){
    if(as.numeric(assembledBand) %in% as.vector(correctBand) || ((as.numeric(assembledBand) - 1) %in% as.vector(correctBand))){
      returnVal = TRUE
    }else{
      returnVal = FALSE
    }
  }else if(bandType == "COMBINED"){
    if((as.numeric(assembledBand) %in% as.vector(correctBand)) || ((as.numeric(assembledBand) - 1) %in% as.vector(correctBand)) || ((as.numeric(assembledBand) - 2) %in% as.vector(correctBand))){
      returnVal = TRUE
    }else{
      returnVal = FALSE
    }
  }

  return(returnVal)
}


## High resolution restoring is always correct and correct plus one
restoreHighRes <- function(assembledBand,correctBand,bandType){

  if(bandType == "INDIVIDUAL"){
    if(!((as.numeric(assembledBand) %in% as.vector(correctBand)) ||
        ((as.numeric(assembledBand) - 1) %in% as.vector(correctBand))) &&
       (as.numeric(assembledBand) < min(as.vector(correctBand)))){
      returnVal = TRUE
    }else{
      returnVal = FALSE
    }
  }else if(bandType == "COMBINED"){
    if(!((as.numeric(assembledBand) %in% as.vector(correctBand)) ||
        ((as.numeric(assembledBand) - 1) %in% as.vector(correctBand)) ||
        ((as.numeric(assembledBand) - 2) %in% as.vector(correctBand))) &&
       ( (as.numeric(assembledBand) - 1) < min(as.vector(correctBand)))){
      returnVal = TRUE
    }else{
      returnVal = FALSE
    }
  }

  return(returnVal)
}

ladderFilter <- function(initialGffFile,resolution,filterType,outputPath,bandProbabilitiesFile,threshold,bandType,polyA,bandName){

  ### Reading the gff containing the assembled transcripts and calculating their lengths
  initialGff <- read.delim(initialGffFile, header = F, comment.char="#",quote = "")
  transcriptList <- lapply(initialGff$V9 , function(x) strsplit(as.character(x),' ')[[1]][4])
  transcriptList<- as.data.frame(unlist(transcriptList))
  initialGff$TranscriptId <- transcriptList$`unlist(transcriptList)`
  initialGff$AssembledBand <- substr(bandName,start = 5,stop = 5)

  onlyExons <- initialGff %>% subset(V3=="exon")
  onlyExons$Length <- abs(onlyExons$V4-onlyExons$V5)
  TranscriptLengths <- onlyExons %>% group_by(TranscriptId) %>% summarise(TotalLength=sum(Length))

  initialGff <- merge(initialGff,TranscriptLengths,by = "TranscriptId")

  ## Correcting for the poly a tails
  initialGff$TotalLength <- initialGff$TotalLength + as.numeric(polyA)

  ## Assigning the correct bands
  initialGff$CorrectBand <- lapply(initialGff$TotalLength,function(x) getBand(x,resolution,bandProbabilitiesFile,threshold))

  if(bandType == "COMBINED"){
    initialGff$AssembledBand <- as.numeric(initialGff$AssembledBand) + 1
  }


  ## Filtering for Low resolution
  if(resolution=="LOW"){
    if(bandType == "INDIVIDUAL"){
      if(filterType=="EQUAL"){
        filteredGtf <- initialGff %>% subset(AssembledBand==CorrectBand)
      }else if(filterType=="LARGER"){
        filteredGtf <- initialGff %>% subset(AssembledBand>=CorrectBand)
      }else if(filterType=="SHORTER"){
        filteredGtf <- initialGff %>% subset(AssembledBand<=CorrectBand)
      }else if(filterType=="CORRECT_PLUS_ONE"){
        filteredGtf <- initialGff %>% subset(AssembledBand==CorrectBand | AssembledBand==(as.numeric(CorrectBand)+1))
      }else if(filterType=="CORRECT_PLUS_TWO"){
        filteredGtf <- initialGff %>% subset(AssembledBand==CorrectBand | AssembledBand==(as.numeric(CorrectBand)+1) | AssembledBand==(as.numeric(CorrectBand)+2))
      }
    }else if(bandType == "COMBINED"){
      ## the calculations are different for the combined band since the assembled band is the higher of the ombined bands
      if(filterType=="EQUAL"){
        filteredGtf <- initialGff %>% subset(AssembledBand==CorrectBand)
      }else if(filterType=="LARGER"){
        filteredGtf <- initialGff %>% subset(AssembledBand>=CorrectBand)
      }else if(filterType=="SHORTER"){
        filteredGtf <- initialGff %>% subset(AssembledBand<=CorrectBand)
      }else if(filterType=="CORRECT_PLUS_ONE"){
        filteredGtf <- initialGff %>% subset(AssembledBand==CorrectBand | AssembledBand==(as.numeric(CorrectBand)+1) | AssembledBand==(as.numeric(CorrectBand)+2))
      }else if(filterType=="CORRECT_PLUS_TWO"){
        filteredGtf <- initialGff %>% subset(AssembledBand==CorrectBand | AssembledBand==(as.numeric(CorrectBand)+1) | AssembledBand==(as.numeric(CorrectBand)+2) | AssembledBand==(as.numeric(CorrectBand)+3))
      }
    }
  }else if(resolution=="HIGH"){
    row <- 1
    for (row in 1:nrow(initialGff)){
      assembledBand <- as.numeric(initialGff[row,11])
      correctBand <- as.vector(initialGff[row,13][[1]])
      toBeKept <- filterHighRes(assembledBand,correctBand,bandType)
      initialGff[row,"ToBeKept"] <- toBeKept
      toBeRestored <- restoreHighRes(assembledBand,correctBand,bandType)
      initialGff[row,"ToBeRestored"] <- toBeRestored
    }

    filteredGtf <- initialGff %>% subset(ToBeKept)
    toBerestoredGtf <- initialGff %>% subset(ToBeRestored)
  }

  # cc <- initialGff %>% subset(initialGff$TranscriptId == "\"STRGcomb12.1445.1\";")
  # ccfiltered <- filteredGtf %>% subset(filteredGtf$TranscriptId == "\"STRGcomb12.1445.1\";")



  ## Writing the filtered gtf
  filteredGtf <- filteredGtf[,c(-1,-11,-12,-13,-14)]
  outFileName <- paste(outputPath,".gtf",sep="")
  write.table(filteredGtf,file = outFileName, sep = "\t",row.names = F, col.names = F, quote = F)


  ## Writing the transcripts which have been filtered out
  toBerestoredGtf <- toBerestoredGtf[,c(-1,-11,-12,-13,-14)]
  outFileName <- paste(outputPath,"_filteredOut.gtf",sep="")
  write.table(toBerestoredGtf,file = outFileName, sep = "\t",row.names = F, col.names = F, quote = F)

}

ladderFilter(initialGffFile,resolution,filterType,outputPath,bandProbabilitiesFile,threshold,bandType,polyA,bandName)
