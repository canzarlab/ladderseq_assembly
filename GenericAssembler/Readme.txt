

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
The following softwares need to be installed in order to run this pipeline.
 
a. kallisto - https://pachterlab.github.io/kallisto/
b. kallisto-ls - https://github.com/canzarlab/kallisto_ladder
c. stringtie2 - https://github.com/skovaka/stringtie2
d. snakemake - https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
e. STAR - https://github.com/alexdobin/STAR
f. R (latest version)
g. Essential R packages:
   dplyr : https://cran.r-project.org/web/packages/dplyr/readme/README.html
   tidyr : https://www.r-project.org/nosvn/pandoc/tidyr.html
   reshape2 : https://www.r-project.org/nosvn/pandoc/reshape2.html
   biomaRt : https://bioconductor.org/packages/release/bioc/html/biomaRt.html
h. Python 3.0
i. Essential python packages:
   biopython : https://biopython.org/wiki/Download
   SeqIo : https://biopython.org/wiki/SeqIO
   sys : https://docs.python.org/3/library/sys.html
   HTSeq : https://pypi.org/project/HTSeq/
j. RNA-seq assembler of your choice
   
   


Step3: Get annotation and cdna ftp paths and prepare snakemake config.py file
---------------------------------------------------------------------------------------
The snakemake pipeline has a config file called config.py where the user needs to specify the directory where the data is stored and also needs to specify the ftp paths for downloading the cdna and genome fasta files corresponding to the data. The following variables need to be set in the config.py file

a. DATA_BASE_PATH : The folder where the data organised as in step 1 is stored.
   e.g. '/user/reads' (if the reads are in the folder called 'reads')

b. CDNA_FTP : the ftp path to download the cdna file from
   e.g. 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' 

c. CDNA_FILE_NAME : the name of the downloaded and unzipped cdna file
   e.g. 'Homo_sapiens.GRCh38.cdna.all.fa'
   
d. DNA_PRIMARY_ASSEMBLY_FTP : the ftp path to download the dna sequence
   e.g. 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

e. ANNOTATION_FASTA_FILE_NAME : the name of the downloaded and unzipped dna sequence file
   e.g. 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'

f. STAR_GENOME_DIR : the path where a STAR index would be stored
   For running the pipeline, the reads need to be aligned to a reference genome. we use STAR for that. STAR needs to create an index. This variable provides the      path where the index is to be stored. This index needs around 30 GB of space for the human genome. Please provide a path with at least 30 GB of space.
   
g. KALLISTO_ORIGINAL : path to the location of the original kallist binary

h. KALLISTO_LADDER : path to kallisto-ls binary


Step5: Fill in commands for the RNA-seq assembler of your choice
---------------------------------------------------------------------------------------
This pipeline works with any program that takes in RNA-seq reads as input and assembled transcripts from them.
The user need to modify the file "Snakefile_runAssembler.smk" and insert the exact coomand for the assembler of choice.
Insert the command for the assembler in the part <YOUR CODE GOES HERE> in the "Snakefile_runAssembler.smk". 
It is essential to name the input and the output files correctly. The input to the assembler should be "{input.alignmentFile}" and the output should be called "{output[0]}".
This modification is there in two places in the "Snakefile_runAssembler.smk" file, once for the individual bands and once for the combined bands.

E.g. If the choice of assembler is scallop then substitute the line <YOUR CODE GOES HERE> with the following:

scallop -i {input.alignmentFile} -o {output[0]}


Step6: Run the snakemake pipeline
---------------------------------------------------------------------------------------
Run the snakemake pipeline using the following command

'snakemake -s Snakemake_runPipeline.smk '

Check snakemake documentation for further details on running snakemake using multiple cores ( -j <number of cores>)

After the pipeline finishes, the assembled transcripts would be in a file called 'LadderSeqAssembly.gtf' in the location 'DATA_BASE_PATH'.

