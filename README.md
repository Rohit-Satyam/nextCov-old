# nextCov
nextCov is a Nextflow based short read variant calling and assembly pipeline for COVID short read paired-end sequencing data. It contains currently 7 independent modules that can be flexibely used in any order. The guidelines to run the pipeline on `ibex` or on `local` computer are as follows:

# UPDATE:
### Date `25-Aug-2022`
1. Included variant summarization step in `08_gatk.nf`.
2. Update README. Added primer file in resources.
### Date `23-Aug-2022`
1. `Resolved `MarkDdduplicateSpark`failing on large BAM files due to memory constraints using `--jobs` option.

### Date `22-Aug-2022`
1. Added the `--help` argument to view help for running each module (`.nf` file).
2. Included `DRAGMAP` indexing option to test `DRAGEN-GATK` pipeline when the `DRAGMAP` issue is resolved.
3. Added script to calculate horizontal and vertical coverage stats post alignment. 
4. Markduplicates is now an independent step. `Multiqc` will be run for steps `FASTQC` and `MARKDUPLICATES` by default. We included `BQSR` step and we save `BEFORE` and `AFTER` recalibration plots. However, `BQSR` might be removed or made optional in the pipeline.
5. User can now specify threads in several steps.
6. Fixed somme minor bugs in `fastqc` and `trimming`steps.


# Organising your data
For nextflow workflow to run without problem you need to organise your data first. This is a common step for both IBEX or your local system. First clone this repository using `git clone https://github.com/Rohit-Satyam/nextCov.git` to download all the scripts. 


```
cd nextCov
```

Make two directories/folders with the following command inside nextCov directory. We already provide these directories and an examplary data in this repository so you can skip the command below. Replace this data/files with your actual data

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
├── primers.bed
├── sequence.fasta
├── sequence.fasta.fai
└── universal_adapter.fasta
```
Once everything is organised your working directory will look like this:
```
nextCov/
├── 00_index.nf
├── 01_fastqc.nf
├── 02_trimming.nf
├── 04_alignbwa.nf
├── 04_ivar.nf
├── 05_markduplicates.nf
├── 06_lowCoverage.nf
├── 07_coverage.nf
├── 08_gatk.nf
├── 09_assembly.nf
├── data
├── multiqc.nf
├── resources
```

> For the moment, organise your data as shown above. The scripts (`.nf` files) must be in same folder as `data` and `resources`

## Running on IBEX

#### Loading Necessary Tools
Running this pipeline on IBEX is super easy since all the modules required are already installed there and you can load all necessary modules using a single liner:

```bash
module load nextflow fastqc multiqc gatk bcftools bwa bedtools trimmomatic tabix covtobed samtools
```
Each `.nf` file is like a module as well  and the instructions and command syntax can be obtained using `--help` flag. Eg:

```
nextflow run 00_index.nf --help
```
> NOTE: mosdepth is not installed on IBEX. A request has been filed for it's installation. `07_coverage.nf` step might therefore fail.

#### Installing tools locally using Conda

```bash
## Install mamba to make download faster
conda install mamba -n base -c conda-forge
## Create separate environment for nextCov installations and activate it
conda create -n nextCov
conda activate nextCov
## Install necessary packages
mamba install -c conda-forge -c bioconda nextflow fastqc multiqc trimmomatic dragmap covtobed mosdepth gatk4 ivar samtools openjdk==8.0.332=h166bdaf_0
```
The commands to run
```
# To index reference fasta file
nextflow run 00_indexbwa.nf --ref "resources/sequence.fasta" --outdir="results/00_indexes"

# Running fastqc on multiple files
nextflow run 01_fastqc.nf --raw 'data/*R{1,2}_001.fastq.gz' --outdir="results/01_rawfastqc"

# Trimming the data (optional) based on fastqc report
nextflow run 02_trimming.nf --raw "data/*R{1,2}_001.fastq.gz" --run "ibex"  

# Fastqc on trimmed reads
nextflow run 01_fastqc.nf --raw "results/02_trimmed/*R{1,2}_P.fastq.gz" --outdir "results/03_trimmedfastqc" 

# Aligning reads
nextflow run 04_alignbwa.nf --raw "results/02_trimmed/*R{1,2}_P.fastq.gz" --cpus 20

# Trimming ARTIC Primers
nextflow run 04_ivar.nf --outdir results/04_ivar_Primertrimmed

# Mark Duplicates
nextflow run 05_markduplicates.nf --raw 'results/04_ivar_Primertrimmed/*.bam' --outdir results/05_markduplicates

# Coverage estimation of regions with less than 30X coverage. These regions will be hard masked (replaced by NNN repeates) later during assembly step
nextflow run 06_lowCoverage.nf --raw 'results/05_markduplicates/*.bam' --outdir results/06_lowcoverageBed

# Horizontal and vertical Coverage estimation (Post Alignment QC)
nextflow run 07_coverage.nf --raw 'results/05_markduplicates/*.bam' --outdir results/07_coverage

# Running variant calling and filtering
nextflow run 08_gatk.nf --raw 'results/05_markduplicates/*.bam' --outdir results/08_gatkresults

# Summarising GATK QC stats
nextflow run multiqc.nf --raw results/08_gatkresults/

# Producing assembly with variants of high confidence
nextflow run 09_assembly.nf --raw 'results/08_gatkresults/vcfs/*.vcf.gz' --outdir results/09_assembly --bed /home/subudhak/Documents/Illumina_ONT_ARTIC_Comparison/illumina/results/06_lowcoverageBed/


```
The resulting directories will be as follows:

```
results/
├── 00_indexes
├── 01_rawfastqc
├── 02_trimmed
├── 03_trimmedfastqc
├── 04_alignments
├── 04_ivar_Primertrimmed
├── 05_markduplicates
├── 06_lowcoverageBed
├── 07_coverage
├── 08_gatkresults
└── 09_assembly
```
> **NOTE:** `*R{1,2}_001.fastq.gz` is a regular expression (a string common to all fastq files). R1_001.fastq.gz and R2_001.fastq.gz are common to all paired end files and therefore can be reduced to `R{1,2}_001.fastq.gz`

### Submitting as a job on ibex
Below given is a ready to use `sbatch` script that can be used to submit the jobs on cluster. Replace my email with yours to receive the job completion update on your email.

```
#SBATCH --nodes=1
#SBATCH --mem=500GB
#SBATCH --partition=batch
#SBATCH --cpus-per-task 30
#SBATCH -J nextCov
#SBATCH -o nextcov.out
#SBATCH -e nextcov.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=rohit.satyam@kaust.edu.sa
#SBATCH --mail-type=ALL

module load nextflow fastqc multiqc gatk bcftools bwa bedtools trimmomatic tabix covtobed samtools

## Copy the commands from above here
```
Save the above file as `nextcov.sbatch` and submit job on ibex using command

```
sbatch nextcov.sbatch
```
