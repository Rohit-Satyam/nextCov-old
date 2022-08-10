# nextCov
nextCov is a Nextflow based short read variant calling and assembly pipeline for COVID short read paired-end sequencing data. It contains currently 7 independent modules that can be flexibely used in any order. The guidelines to run the pipeline on `ibex` or on `local` computer are as follows:

# Organising your data
For nextflow workflow to run without problem you need to organise your data first. This is a common step for both IBEX or your local system. Make two directories/folders with the following command:
```bash
mkdir data resources
```
The command `mkdir` will make directories `data` and `resources`. The `data` directory contains your sequencing data (paired end reads) and resources directory contains the **reference genome** fasta file, adapter sequence fasta file used by trimmomatic and other additional files that you might require for your downstream analysis. An example of the contents of `resources` directory is shown below:

```bash
tree resources/
```

```bash
resources/
├── exclude.bed
├── sequence.fasta
├── sequence.fasta.fai
└── universal_adapter.fasta
```
Now clone this repository using `git clone `

## Running on IBEX
Running this pipeline on IBEX is super easy since all the modules required are already installed there and you can load all necessary modules using a single liner:

```bash
module load nextflow fastqc multiqc gatk bcftools bwa bedtools trimmomatic tabix covtobed samtools
```



