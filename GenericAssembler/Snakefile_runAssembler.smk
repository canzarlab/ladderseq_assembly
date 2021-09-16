### Snakemake file for ###
### 1. Single e ###

### Author: Shounak ###
### 23.07.2018 ###

rule run_assembler_individual:
  input:
    alignmentFile = DATA_BASE_PATH+"/{individualBand}/Aligned.sortedByCoord.out.bam",
    bandProb = DATA_BASE_PATH+"/migProb.tsv"
  threads:
    N_THREADS
  output:
    DATA_BASE_PATH+"/{individualBand}/individual_assembly.gtf",
    DATA_BASE_PATH+"/{individualBand}/individual_assembly_filteredOut.gtf"
  run:
    shell("mkdir -p {DATA_BASE_PATH}/{wildcards.individualBand}")
    # Enter the command for running the assembler of your choice in the line below
    # The input file should be {input.alignmentFile} and the output file, {output[0]}
    # E.G. for Scallop : scallop -i {input.alignmentFile} -o {output[0]}
    # E.G. for Stringtie2 : stringtie2 {input.alignmentFile} -o {output[0]}
    shell("<YOUR CODE GOES HERE>")
    shell("sed -i '/^#/ d' {output[0]}")
    shell("Rscript GffLadderFilter.R {output[0]} HIGH CORRECT_PLUS_ONE {DATA_BASE_PATH}/{wildcards.individualBand}/filteredAssembly_individual {input.bandProb} 1 INDIVIDUAL 200 {wildcards.individualBand}")
    shell("mv {DATA_BASE_PATH}/{wildcards.individualBand}/filteredAssembly_individual.gtf {output[0]}")
    shell("mv {DATA_BASE_PATH}/{wildcards.individualBand}/filteredAssembly_individual_filteredOut.gtf {output[1]}")


rule run_assembler_combined:
  input:
    alignmentFile = DATA_BASE_PATH+"/{combBand}/Comb.Aligned.sortedByCoord.out.bam",
    bandProb = DATA_BASE_PATH+"/migProb.tsv"
  threads:
    N_THREADS
  output:
    DATA_BASE_PATH+"/{combBand}/combined_assembly.gtf",
    DATA_BASE_PATH+"/{combBand}/combined_assembly_filteredOut.gtf"
  run:
    shell("mkdir -p {DATA_BASE_PATH}/{wildcards.combBand}")
    # Enter the command for running the assembler of your choice in the line below
    # The input file should be {input.alignmentFile} and the output file, {output[0]}
    # E.G. for Scallop : scallop -i {input.alignmentFile} -o {output[0]}
    # E.G. for Stringtie2 : stringtie2 {input.alignmentFile} -o {output[0]}
    shell("<YOUR CODE GOES HERE>")
    shell("sed -i '/^#/ d' {output[0]}")
    shell("Rscript GffLadderFilter.R {output[0]} HIGH CORRECT_PLUS_ONE {DATA_BASE_PATH}/{wildcards.combBand}/filteredAssembly_combined {input.bandProb} 1 COMBINED 200 {wildcards.combBand}")
    shell("mv {DATA_BASE_PATH}/{wildcards.combBand}/filteredAssembly_combined.gtf {output[0]}")
    shell("mv {DATA_BASE_PATH}/{wildcards.combBand}/filteredAssembly_combined_filteredOut.gtf {output[1]}")
