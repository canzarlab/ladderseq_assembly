### Snakemake configuration file for running ladder-assembly pipeline using any RNA-Seq assembler###
### Author: Shounak Chakraborty###
### 20.08.2020 ###


# The path to your data folder
# The folder where the reads are stored. e.g. '/reads' (if the reads are in the folder called 'reads')
DATA_BASE_PATH = ''


# Annotation Paths
#  The name of the downloaded and unzipped cdna file e.g. 'Homo_sapiens.GRCh38.cdna.all.fa'
CDNA_FASTA_FILE = ''
#  The name of the downloaded and unzipped dna sequence file e.g. 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
ANNOTATION_FASTA_FILE = ''

# the path where a STAR index would be stored For running the pipeline, the reads need to be aligned to a reference genome. we use STAR for that.
# STAR needs to create an index. This variable provides the path where the index is to be stored.
# This index needs around 30 GB of space for the human genome. Please provide a path with at least 30 GB of space.
STAR_GENOME_DIR = DATA_BASE_PATH+"/starGenomeDirectory"




N_THREADS = 8


# software configurations
# The following binaries are provided in the directory ext. If the user wants to use a binary different from the ones provided, then these vriables need to be changed.
STRINGTIE2_BINARY = '../ext/linuxBinaries/stringtie'
STAR_BINARY = '../ext/linuxBinaries/STAR'
GFFCOMPARE_BINARY = '../ext/linuxBinaries/gffcompare'
GTFFILTER_BINARY = '../ext/linuxBinaries/gtfFilter'
GFFREAD_BINARY = '../ext/linuxBinaries/gffread'
KALLISTO_LS_BINARY = '../ext/linuxBinaries/kallisto-ls'
SAMTOOLS_BINARY = '../ext/linuxBinaries/samtools'
EXONREFINE_BINARY = '../ext/linuxBinaries/exonRefine'


# Parameters for Ladder seq assembly
INDIVIDUAL_BAND = ["band1","band2","band3","band4","band5","band6","band7"]
COMBINED_BAND = ["comb12","comb23","comb34","comb45","comb56","comb67"]
KAL_IDX = DATA_BASE_PATH+'/kallisto_index.idx'
