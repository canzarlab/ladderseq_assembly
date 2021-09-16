### Snakemake file for ###
### 1. Single e ###

### Author: Shounak ###
### 23.07.2018 ###


# MOD_KALLISTO_FILTER = 0.1
# SUBEXON_LENGTH_THRESHOLD = 2
# RESTORE = 'TRUE'
# RESTORED_TPM_FILTER = 1

def determineBands(wildcards):
    bands = wildcards.combBand

    if bands=="comb12":
        returnList = [str(DATA_BASE_PATH+"/band1/individual_assembly.gtf"),str(DATA_BASE_PATH+"/band2/individual_assembly.gtf")]
    elif bands=="comb23":
        returnList = [str(DATA_BASE_PATH+"/band2/individual_assembly.gtf"),str(DATA_BASE_PATH+"/band3/individual_assembly.gtf")]
    elif bands=="comb34":
        returnList = [str(DATA_BASE_PATH+"/band3/individual_assembly.gtf"),str(DATA_BASE_PATH+"/band4/individual_assembly.gtf")]
    elif bands=="comb45":
        returnList = [str(DATA_BASE_PATH+"/band4/individual_assembly.gtf"),str(DATA_BASE_PATH+"/band5/individual_assembly.gtf")]
    elif bands=="comb56":
        returnList = [str(DATA_BASE_PATH+"/band5/individual_assembly.gtf"),str(DATA_BASE_PATH+"/band6/individual_assembly.gtf")]
    elif bands=="comb67":
        returnList = [str(DATA_BASE_PATH+"/band6/individual_assembly.gtf"),str(DATA_BASE_PATH+"/band7/individual_assembly.gtf")]

    returnList.append(str(DATA_BASE_PATH+"/"+wildcards.combBand+"/combined_assembly.gtf"))

    return returnList


rule combineAssemblies:
    input:
        determineBands,
    output:
        DATA_BASE_PATH+"/{combBand}/pair.combined.gtf"
    shell:
        """
            mkdir -p {DATA_BASE_PATH}/{wildcards.combBand}
            {GFFCOMPARE_BINARY} -e 0 -o {DATA_BASE_PATH}/{wildcards.combBand}/pair -C {input}
        """


rule run_gffCompare_merge:
    input:
        expand(DATA_BASE_PATH+"/{combBand}/pair.combined.gtf", combBand=COMBINED_BAND)
    output:
        temp(DATA_BASE_PATH+"/gffMerged.combined.gtf"),
        DATA_BASE_PATH+"/gffMerged/gffMergedTranscripts.gtf"
    shell:
        """
        mkdir -p {DATA_BASE_PATH}/gffMerged
        {GFFCOMPARE_BINARY} -o {DATA_BASE_PATH}/gffMerged {input}
        grep -v '^#' {output[0]} > {output[1]}
        """


rule run_gffCompare_mergeRejectedTranscriptsToOneFile:
    input:
        expand(DATA_BASE_PATH+"/{band}/combined_assembly_filteredOut.gtf", band=COMBINED_BAND),
        expand(DATA_BASE_PATH+"/{band}/individual_assembly_filteredOut.gtf", band=INDIVIDUAL_BAND)
    output:
        temp(DATA_BASE_PATH+"/rejectedTranscripts.combined.gtf"),
        DATA_BASE_PATH+"/rejectedTranscripts_combined.gtf",

    shell:
        """
        {GFFCOMPARE_BINARY} -D -o {DATA_BASE_PATH}/rejectedTranscripts {input}
        grep -v '^#' {output[0]} > {output[1]}
        """


rule run_stringtieMerge_gffMerge:
	input:
		DATA_BASE_PATH+"/gffMerged/gffMergedTranscripts.gtf"
	output:
		temp(DATA_BASE_PATH+"/concatStringtie.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/concatStringtieMerged.gtf")
	shell:
		"""
		mkdir -p {DATA_BASE_PATH}/LadderSeqAssembly
		{STRINGTIE2_BINARY} --merge -o {output[0]} {input}
		grep -v '^#' {output[0]} > {output[1]}
		"""


rule run_gffCompare_compareGffandStringtieMerge:
	input:
		DATA_BASE_PATH+"/gffMerged/gffMergedTranscripts.gtf",
		DATA_BASE_PATH+"/LadderSeqAssembly/concatStringtieMerged.gtf"
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/compareGff_str_merge.concatStringtieMerged.gtf.refmap"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/compareGff_str_merge.concatStringtieMerged.gtf.tmap")
	shell:
		"""
		{GFFCOMPARE_BINARY} -r {input[0]} -o {DATA_BASE_PATH}/compareGff_str_merge {input[1]}
		"""


rule run_generate_sets:
	input:
		DATA_BASE_PATH+"/gffMerged/gffMergedTranscripts.gtf",
		DATA_BASE_PATH+"/LadderSeqAssembly/compareGff_str_merge.concatStringtieMerged.gtf.tmap"
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/setG.txt"), ## single exon transcripts in gffMerged assembly
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/setS.txt"), ## single exon transcripts in strMerge which are equal to single exon transcripts in gffMerged
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/setU.txt"), ## set of unknown single exons in stringtie merge
	shell:
		"""
		sed 's/\s/\t/g' {input[0]} | awk '$3 == "exon"'| datamash -s  -g 10 count 3 | awk '$2 == 1 {{print $1}}' | sed 's/\"//g' | sed 's/\;//g' > {output[0]}
		awk '($3 == "=" && $6 == 1) {{print $2}}' {input[1]} | sort > {output[1]}
		awk '($3 == "u" && $6 == 1) {{print $5}}' {input[1]} | sort > {output[2]}
		"""

rule run_setOperations_SetDifference:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/setG.txt",
		DATA_BASE_PATH+"/LadderSeqAssembly/setS.txt",
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/setA.txt"),
	shell:
		"""
		comm -23 {input[0]} {input[1]} > {output[0]}
		"""


rule run_filterSingleExons_gffFilter:
	input:
		DATA_BASE_PATH+"/gffMerged/gffMergedTranscripts.gtf",
		DATA_BASE_PATH+"/LadderSeqAssembly/setA.txt",
		DATA_BASE_PATH+"/LadderSeqAssembly/concatStringtieMerged.gtf",
		DATA_BASE_PATH+"/LadderSeqAssembly/setU.txt"
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_withSubfrags.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_C.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_X.gtf")
	shell:
		"""
			{GTFFILTER_BINARY} --list {input[1]} -m blackT {input[0]} {output[0]}
			{GFFCOMPARE_BINARY} -e 0 -o {DATA_BASE_PATH}/LadderSeqAssembly/discardAllSubfrags_C -C {output[0]}
			mv {DATA_BASE_PATH}/LadderSeqAssembly/discardAllSubfrags_C.combined.gtf {output[1]}
			{GFFCOMPARE_BINARY} -e 0 -o {DATA_BASE_PATH}/LadderSeqAssembly/discardAllSubfrags_X -X {output[0]}
			mv {DATA_BASE_PATH}/LadderSeqAssembly/discardAllSubfrags_X.combined.gtf {output[2]}
			awk '{{for (i=1; i<=NF; ++i) {{ if ($i ~ "transcript_id") print $(i+1)}} }}' {output[1]} | sed 's/\"//g' | sed 's/\;//g' | sort | uniq > {DATA_BASE_PATH}/LadderSeqAssembly/transcriptList_C.txt
			awk '{{for (i=1; i<=NF; ++i) {{ if ($i ~ "transcript_id") print $(i+1)}} }}' {output[2]} | sed 's/\"//g' | sed 's/\;//g' | sort | uniq > {DATA_BASE_PATH}/LadderSeqAssembly/transcriptList_X.txt
			comm -23 {DATA_BASE_PATH}/LadderSeqAssembly/transcriptList_C.txt {DATA_BASE_PATH}/LadderSeqAssembly/transcriptList_X.txt > {DATA_BASE_PATH}/LadderSeqAssembly/presentInCNotX.txt
		"""

rule run_subExonFilter:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_C.gtf"
	threads:
	        N_THREADS
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly.gtf")
	shell:
		"""
		bash FindUniqSubExons.sh {input[0]} {DATA_BASE_PATH}/LadderSeqAssembly/ 2 {EXONREFINE_BINARY}
		comm -12 {DATA_BASE_PATH}/LadderSeqAssembly/transcripts_uniqueSubexonLength_threshold.txt {DATA_BASE_PATH}/LadderSeqAssembly/presentInCNotX.txt > {DATA_BASE_PATH}/LadderSeqAssembly/uniqueNotInX.txt
		{GTFFILTER_BINARY} --list {DATA_BASE_PATH}/LadderSeqAssembly/uniqueNotInX.txt -m blackT {input[0]} {output[0]}
		"""



rule make_listOfTranscriptsToBeKept:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly.gtf",
		DATA_BASE_PATH+"/rejectedTranscripts_combined.gtf"
	threads:
		N_THREADS
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/restored_rejectedTrans_tagged.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/TranscriptsToBePutBack.txt"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/TranscriptsToBePutBack_sorted.txt"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/TranscriptsToBePutBack.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly.restored.combined.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_tagged.gtf"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_restored.gtf")
	shell:
		"""
			sed 's/TCONS_/TCONS_rejected_/g' {input[1]} > {output[0]}
			python MergeRejectedTranscripts.py {input[0]} {output[0]} {output[1]}
			sort {output[1]} > {output[2]}
			{GTFFILTER_BINARY} --list {output[2]} -m whiteT {output[0]} {output[3]}
			sed -i 's/TCONS_/TCONS_restored/g' {output[3]}
			{GFFCOMPARE_BINARY} -o {DATA_BASE_PATH}/LadderSeqAssembly/singleExonsRemovedInitialAssembly.restored {input[0]} {output[3]}
			Rscript TagRestoredTranscripts.R {output[4]} {output[5]}
			grep -v '^#' {output[5]} > {output[6]}
		"""





rule run_gffRead_SingleExon:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_restored.gtf",
		ANNOTATION_FASTA_FILE
	threads:
	        N_THREADS
	output:
		DATA_BASE_PATH+"/LadderSeqAssembly/firstRoundAssembly_singleExons.fa"
	shell:
		"""
		{GFFREAD_BINARY}  {input[0]} -g {input[1]} -w {output[0]}
		"""



rule run_kallistoMod_SingleExon_index:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/firstRoundAssembly_singleExons.fa"
	threads:
	        N_THREADS
	output:
		DATA_BASE_PATH+"/LadderSeqAssembly/modKallistoIndex_singleExons.idx"
	shell:
		"""
		{KALLISTO_LS_BINARY} index -i {output[0]} {input[0]}
		"""


rule run_kallisto_SingleExon_modified:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/modKallistoIndex_singleExons.idx",
		ladderReads = expand(DATA_BASE_PATH+"/{band}/R{pair}_sim.fq",band=INDIVIDUAL_BAND,pair=[1,2])
	threads:
	        N_THREADS
	output:
		DATA_BASE_PATH+"/LadderSeqAssembly/abundance.tsv"
	shell:
		"""
		{KALLISTO_LS_BINARY} quant -t {threads} -o {DATA_BASE_PATH}/LadderSeqAssembly/ -i {input[0]} {input.ladderReads}
		"""






rule run_filterModKallistoGff_SingleExon:
	input:
		unfilteredGff = DATA_BASE_PATH+"/LadderSeqAssembly/singleExonsRemovedInitialAssembly_restored.gtf",
		modKallistoQuantResult = DATA_BASE_PATH+"/LadderSeqAssembly/abundance.tsv"
	threads:
	        N_THREADS
	output:
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/restoredTranscripts_tobeKept.txt"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/initialAssemblyTranscripts_tobeKept.txt"),
		temp(DATA_BASE_PATH+"/LadderSeqAssembly/transcripts_tobeKept.txt"),
		DATA_BASE_PATH+"/LadderSeqAssembly/LadderSeqAssembly.gtf"
	shell:
		"""
		grep "TCONS_restored" {input[1]} | awk '$5 >= 1 {{print $1}}' > {output[0]}
		grep -v "TCONS_restored" {input[1]} | awk '$5 >= 0.1 {{print $1}}' > {output[1]}
		cat {output[0]} {output[1]} > {output[2]}
		{GTFFILTER_BINARY} --list {output[2]} -m whiteT {input[0]} {output[3]}
		"""


rule cleanup:
	input:
		DATA_BASE_PATH+"/LadderSeqAssembly/LadderSeqAssembly.gtf"
	threads:
	        N_THREADS
	output:
		DATA_BASE_PATH+"/LadderSeqAssembly.gtf"
	shell:
		"""
        mv {DATA_BASE_PATH}/LadderSeqAssembly/LadderSeqAssembly.gtf {DATA_BASE_PATH}/
        rm -r {DATA_BASE_PATH}/LadderSeqAssembly
        rm -r {DATA_BASE_PATH}/band*
        rm -r {DATA_BASE_PATH}/comb*
        rm -r {DATA_BASE_PATH}/gffMerged
        rm {DATA_BASE_PATH}/gffMerged*
        rm {DATA_BASE_PATH}/compare*
        rm {DATA_BASE_PATH}/rejected*
		"""
