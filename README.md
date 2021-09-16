# Reference-based assembly from Ladder-seq

![Ladder-Seq Assembly pipeline](/AssemblyImage400.png)

**Ladder-seq**  is the concerted advancement of the RNA-seq protocol and its computational methods. It experimentally separates transcripts according to their length prior to sequencing to achieve a "coloring" of reads that connects them along transcript isoforms.

In the **reference-based transcript assembly** from Ladder-seq reads we use [StringTie2](https://ccb.jhu.edu/software/stringtie/) to build a separate splicing graph in each band and apply transcript length constraints derived from read colorings to break too long erroneous fusions and to eliminate too short transcript fragments. Finally, reads are assigned to assembled transcripts guided by their coloring using [kallisto-ls](https://github.com/canzarlab/kallisto-ls), our extension of [kallisto](https://pachterlab.github.io/kallisto/about) to Ladder-seq.

 <br />



## Installation

This repository includes precompiled StringTie 2.1.4, read aligner [STAR 2.7.5c](https://github.com/alexdobin/STAR), [samtools 1.10]() and kallisto-ls binaries for Linux x86_64. It includes kallisto-ls also as a submodule that is automatically included when obtaining this repository through command:
```shell
git clone --recursive https://github.com/canzarlab/LadderSeq-Assembly.git
```
Then, kallisto-ls can be re-built from source as described [here](https://github.com/canzarlab/kallisto-ls) and the path to the kallist-ls binary adjusted in the `config.py` file (see below).
All required steps, including read alignment (using STAR), estimation of migration patterns, length-contrained transcript assembly, integration, and quantification are run by the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system, which needs to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on your system.

## Transcript assembly using StringTie2

We refer to the StringTie2 based assembly workflow for Ladder-seq as __StringTie-ls__. It assumes that Ladder-seq input reads for each band are stored in a common directory with the following file name convention: ```R1_band<i>.fq``` and ```R2_band<i>.fq``` for paired-end reads in band ```i```. __StringTie-ls__ currently assumes a separation of reads into 7 bands.

In file `config.py` in the Stringtie-ls subdirectory, the user needs to set variables for the directory path to the input read files, and to cDNA and genome sequences in FASTA format. More detailed instructions are provided in file `config.py`. To run __StringTie-ls__, simply execute the Snakemake workflow:


```shell
cd Stringtie-ls
snakemake -s Snakemake_runPipeline.smk
```
The number of cores to be used by Snakemake can be specified using `-j <number of cores>`, which is described in more detailed in the Snakemake [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options).

Assembled transcripts are reported in file `LadderSeqAssembly.gtf` in the directory specified by variable `DATA_BASE_PATH` in `config.py`.

## Transcript assembly using an alternative short-read assembler

In directory `GenericAssembler` we provide a more generic workflow in which StringTie2 can easily be replaced by a short-read RNA-seq assembly method of the user's choice.
In file `Snakefile_runAssembler.smk` the user needs to replace two occurrences of 
```shell
<YOUR CODE GOES HERE>
``` 
with the command used to call the preferred assemby method, using
`{input.alignmentFile}` and  `{output[0]}` to specify input read alignments in `.bam` format and output file containing transcripts in `.gtf` format, respectively.

As described above, set variables in `config.py` in directory `GenericAssembler` and run the Snakemake workflow using

```shell
cd GenericAssembler
snakemake -s Snakemake_runPipeline.smk
```
