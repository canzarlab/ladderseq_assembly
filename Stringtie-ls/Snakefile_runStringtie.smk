### Snakemake file for ###
### 1. Single e ###

### Author: Shounak ###
### 23.07.2018 ###

# FILTER_TYPE_INDIVIDUAL="CORRECT_PLUS_ONE"
# FILTER_TYPE_COMBINED="CORRECT_PLUS_ONE"
# LADDER_RESOLUTION = "HIGH"
# LADDER_HIGH_FILTER_THRESHOLD = 1
# POLY_A = 200


rule run_stringtie_individual:
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
    shell("{STRINGTIE2_BINARY} {input.alignmentFile} -o {output[0]}")
    shell("sed -i '/^#/ d' {output[0]}")
    shell("Rscript GffLadderFilter.R {output[0]} HIGH CORRECT_PLUS_ONE {DATA_BASE_PATH}/{wildcards.individualBand}/filteredAssembly_individual {input.bandProb} 1 INDIVIDUAL 200 {wildcards.individualBand}")
    shell("mv {DATA_BASE_PATH}/{wildcards.individualBand}/filteredAssembly_individual.gtf {output[0]}")
    shell("mv {DATA_BASE_PATH}/{wildcards.individualBand}/filteredAssembly_individual_filteredOut.gtf {output[1]}")


rule run_stringtie_combined:
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
    shell("{STRINGTIE2_BINARY} {input.alignmentFile} -o {output[0]}")
    shell("sed -i '/^#/ d' {output[0]}")
    shell("Rscript GffLadderFilter.R {output[0]} HIGH CORRECT_PLUS_ONE {DATA_BASE_PATH}/{wildcards.combBand}/filteredAssembly_combined {input.bandProb} 1 COMBINED 200 {wildcards.combBand}")
    shell("mv {DATA_BASE_PATH}/{wildcards.combBand}/filteredAssembly_combined.gtf {output[0]}")
    shell("mv {DATA_BASE_PATH}/{wildcards.combBand}/filteredAssembly_combined_filteredOut.gtf {output[1]}")
