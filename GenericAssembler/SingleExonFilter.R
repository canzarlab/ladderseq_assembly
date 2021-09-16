## libraries
require(dplyr)




# setSFile <- "/algbio1/shounak/raw/simulated/human/ImperfectSep100HighRes/AssemblyByCoverage/tpm0.1/substituted/onlyIndStringLowKallistoFilter0.1/test/setS.txt"
# setGFile <- "/algbio1/shounak/raw/simulated/human/ImperfectSep100HighRes/AssemblyByCoverage/tpm0.1/substituted/onlyIndStringLowKallistoFilter0.1/test/setG.txt"
# mappingFile <- "/algbio1/shounak/raw/simulated/human/ImperfectSep100HighRes/AssemblyByCoverage/tpm0.1/substituted/onlyIndStringLowKallistoFilter0.1/test/oldId_concatId_map.txt"
# outFileName <- "/algbio1/shounak/raw/simulated/human/ImperfectSep100HighRes/AssemblyByCoverage/tpm0.1/substituted/onlyIndStringLowKallistoFilter0.1/test/setA.txt"



#### Command line arguments
args = commandArgs(trailingOnly=TRUE)

setSFile <- args[1]
setGFile <- args[2]
mappingFile <- args[3]
## Output path
outFileName <- args[4]



subtractSets <- function(setSFile,setGFile,mappingFile,outFileName){
  setS <- read.delim(setSFile, header = F, comment.char="#",quote = "")
  names(setS) <- "newId"
  setG <- read.delim(setGFile, header = F, comment.char="#",quote = "")
  names(setG) <- "oldId"
  mapping <- read.delim(mappingFile, header = F, comment.char="#",quote = "",sep = " ")
  names(mapping) <- c("newId","oldId")


  setG <- merge(setG,mapping,by = "oldId")


  setG_minus_S <- setG %>% subset(!newId %in% setS$newId )

  outputSet <- as.data.frame(unique(setG_minus_S$newId))
  names(outputSet) <- "newId"


  write.table(outputSet,file=outFileName, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}



subtractSets(setSFile,setGFile,mappingFile,outFileName)
