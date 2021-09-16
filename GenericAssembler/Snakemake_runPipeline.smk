### Snakemake file for ###
### 1. Running STAR on different bands###
### 2. Running Stringtie ###

### Author: Shounak ###
### 30.07.2018 ###

include: "config.py"
include: "Snakefile_prepareData.smk"
include: "Snakefile_runAssembler.smk"
include: "Snakefile_runLadderAssembly.smk"




rule all:
	input:
		expand(DATA_BASE_PATH+"/LadderSeqAssembly.gtf"),
