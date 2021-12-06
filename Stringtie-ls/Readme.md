
Pipeline for assembling transcripts using Ladder-seq data.


Step1: Prepare Data
---------------------------------------------------------------------------------------
Paired end ladder-seq reads need to be organised in the following format in a directory.
R1_band1.fq, R2_band1.fq
R1_band2.fq, R2_band2.fq
.
.
.
R1_band7.fq, R2_band7.fq
(Right now the pipeline is implemented for only 7 bands).


Step2: Download software
---------------------------------------------------------------------------------------
The following software needs to be installed in order to run this pipeline.

a. kallisto-ls - https://github.com/canzarlab/kallisto_ladder
b. stringtie2 - https://github.com/skovaka/stringtie2
c. snakemake - https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
d. STAR - https://github.com/alexdobin/STAR
e. R (latest version)
f. Essential R packages:
   dplyr : https://cran.r-project.org/web/packages/dplyr/readme/README.html
   tidyr : https://www.r-project.org/nosvn/pandoc/tidyr.html
   reshape2 : https://www.r-project.org/nosvn/pandoc/reshape2.html
   biomaRt : https://bioconductor.org/packages/release/bioc/html/biomaRt.html
g. Python 3.0
h. Essential python packages:
   biopython : https://biopython.org/wiki/Download
   SeqIo : https://biopython.org/wiki/SeqIO
   sys : https://docs.python.org/3/library/sys.html
   HTSeq : https://pypi.org/project/HTSeq/




Step3: Get annotation and cdna files and prepare snakemake config.py file
---------------------------------------------------------------------------------------
The snakemake pipeline has a config file called config.py where the user needs to specify the directory where the data is stored and also needs to specify the paths and filenames of the cdna and genome fasta files corresponding to the data. The following variables need to be set in the config.py file

a. DATA_BASE_PATH : The folder where the data is organised as in step 1 is stored.
   e.g. '/user/reads' (if the reads are in the folder called 'reads')

b. CDNA_FILE_NAME : the path and name of the downloaded and unzipped cdna file
   e.g. 'path/Homo_sapiens.GRCh38.cdna.all.fa'

c. ANNOTATION_FASTA_FILE_NAME : the name of the downloaded and unzipped dna sequence file
   e.g. 'path/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

d. Path to the respective binaries are initially set to the '../ext/linuxBinaries/' folder. The user can change these to a binary of choice.



Step4: Run the snakemake pipeline
---------------------------------------------------------------------------------------
Run the snakemake pipeline using the following command

'snakemake -s Snakemake_runPipeline.smk '

Check snakemake documentation for further details on running snakemake using multiple cores ( -j <number of cores>)

After the pipeline finishes, the assembled transcripts would be in a file called 'LadderSeqAssembly.gtf' in the location 'DATA_BASE_PATH'.
