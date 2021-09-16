### Snakemake file for preparing/downloading data for the ladder-assembly pipeline  ###

### Author: Shounak Chakraborty###
### 20.08.2020 ###


rule runSTARGenomeGenerate:
    input:
        ANNOTATION_FASTA_FILE
    output:
        STAR_GENOME_DIR
    shell:
        "mkdir -p {STAR_GENOME_DIR};"
        "{STAR_BINARY} --runThreadN 2 --runMode genomeGenerate --genomeDir {STAR_GENOME_DIR} --genomeFastaFiles {ANNOTATION_FASTA_FILE}"


rule kallisto_index:
    input:
        CDNA_FASTA_FILE
    output:
        KAL_IDX
    threads: 1
    shell:
        KALLISTO_LS_BINARY + ' index -i ' + KAL_IDX + ' ' + CDNA_FASTA_FILE


rule organizeReads:
    input:
        DATA_BASE_PATH+"/R1_{individualBand}.fq",
        DATA_BASE_PATH+"/R2_{individualBand}.fq"
    threads:
        N_THREADS
    output:
        DATA_BASE_PATH+"/{individualBand}/R1_sim.fq",
        DATA_BASE_PATH+"/{individualBand}/R2_sim.fq"
    run:
        shell("mkdir -p {DATA_BASE_PATH}/{wildcards.individualBand}/")
        shell("cp {input[0]} {DATA_BASE_PATH}/{wildcards.individualBand}/R1_sim.fq")
        shell("cp {input[1]} {DATA_BASE_PATH}/{wildcards.individualBand}/R2_sim.fq")



rule estimateMigProbs:
    input:
        KAL_IDX,
        ladderReads = expand(DATA_BASE_PATH+"/{band}/R{pair}_sim.fq",band=INDIVIDUAL_BAND,pair=[1,2])
    threads:
        N_THREADS
    output:
        DATA_BASE_PATH+"/migProb.tsv"
    run:
        shell("{KALLISTO_LS_BINARY} migrate -o {DATA_BASE_PATH}/ -i {input}")


def determineBands(wildcards):
    bandCombination = wildcards.combBand
    pair = wildcards.pair
    if bandCombination=="comb12":
        returnList = [str(DATA_BASE_PATH+"/band1/R"+pair+"_sim.fq"),str(DATA_BASE_PATH+"/band2/R"+pair+"_sim.fq")]
    elif bandCombination=="comb23":
        returnList = [str(DATA_BASE_PATH+"/band2/R"+pair+"_sim.fq"),str(DATA_BASE_PATH+"/band3/R"+pair+"_sim.fq")]
    elif bandCombination=="comb34":
        returnList = [str(DATA_BASE_PATH+"/band3/R"+pair+"_sim.fq"),str(DATA_BASE_PATH+"/band4/R"+pair+"_sim.fq")]
    elif bandCombination=="comb45":
        returnList = [str(DATA_BASE_PATH+"/band4/R"+pair+"_sim.fq"),str(DATA_BASE_PATH+"/band5/R"+pair+"_sim.fq")]
    elif bandCombination=="comb56":
        returnList = [str(DATA_BASE_PATH+"/band5/R"+pair+"_sim.fq"),str(DATA_BASE_PATH+"/band6/R"+pair+"_sim.fq")]
    elif bandCombination=="comb67":
        returnList = [str(DATA_BASE_PATH+"/band6/R"+pair+"_sim.fq"),str(DATA_BASE_PATH+"/band7/R"+pair+"_sim.fq")]

    return returnList



rule combineReads:
    input:
        determineBands
    output:
        DATA_BASE_PATH+"/{combBand}/Comb_R{pair}_sim.fq"
    shell:
        """
            mkdir -p {DATA_BASE_PATH}/{wildcards.combBand}
            cat {input} > {DATA_BASE_PATH}/{wildcards.combBand}/Comb_R{wildcards.pair}_sim.fq
        """



rule run_STAR_Band_INDIVIDUAL:
    input:
        genomeDir = STAR_GENOME_DIR,
        inputFilePair1 = DATA_BASE_PATH+"/{band}/R1_sim.fq",
        inputFilePair2 = DATA_BASE_PATH+"/{band}/R2_sim.fq"
    threads:
        N_THREADS
    output:
        DATA_BASE_PATH+"/{band}/Aligned.sortedByCoord.out.bam"
    shell:
        """
        mkdir -p {DATA_BASE_PATH}/{wildcards.band}
        {STAR_BINARY} --runThreadN 2 --genomeDir {input.genomeDir} --readFilesIn {input.inputFilePair1} {input.inputFilePair2} --outSAMstrandField intronMotif --outFileNamePrefix {DATA_BASE_PATH}/{wildcards.band}/ --outSAMtype BAM SortedByCoordinate
        """



rule run_STAR_Band_COMBINED:
        input:
                genomeDir = STAR_GENOME_DIR,
                inputFilePair1 = DATA_BASE_PATH+"/{combBand}/Comb_R1_sim.fq",
                inputFilePair2 = DATA_BASE_PATH+"/{combBand}/Comb_R2_sim.fq"
        threads:
                N_THREADS
        output:
                DATA_BASE_PATH+"/{combBand}/Comb.Aligned.sortedByCoord.out.bam"
        shell:
                """
                {STAR_BINARY} --runThreadN 2 --genomeDir {input.genomeDir} --readFilesIn {input.inputFilePair1} {input.inputFilePair2} --outSAMstrandField intronMotif --outFileNamePrefix {DATA_BASE_PATH}/{wildcards.combBand}/Comb. --outSAMtype BAM SortedByCoordinate
                """
