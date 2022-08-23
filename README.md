# nextCov
nextCov is a Nextflow based short read variant calling and assembly pipeline for COVID short read paired-end sequencing data. It contains currently 7 independent modules that can be flexibely used in any order. The guidelines to run the pipeline on `ibex` or on `local` computer are as follows:

# UPDATE:
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
├── sequence.fasta
├── sequence.fasta.fai
└── universal_adapter.fasta
```
Once everything is organised your working directory will look like this:
```
nextCov/
├── 00_indexbwa.nf
├── 01_fastqc.nf
├── 02_trimming.nf
├── 04_alignbwa.nf
├── 05_lowCoverage.nf
├── 06_gatk.nf
├── 07_assembly.nf
├── data
├── resources

```

> For the moment, organise your data as shown above. The scripts (`.nf` files) must be in same folder as `data` and `resources`

## Running on IBEX
Running this pipeline on IBEX is super easy since all the modules required are already installed there and you can load all necessary modules using a single liner:

```bash
module load nextflow fastqc multiqc gatk bcftools bwa bedtools trimmomatic tabix covtobed samtools
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
nextflow run 04_alignbwa.nf --raw "results/02_trimmed/*R{1,2}_P.fastq.gz" 

# Identifying regions of low coverage to mask in the final assembly
nextflow run 05_lowCoverage.nf --raw 'results/04_alignments/*bam'

# Running GATK variant calling
nextflow run 06_gatk.nf

# Making final assembly
nextflow run 07_assembly.nf
```
The resulting directories will be as follows:

```
results/
├── 00_indexes
├── 01_rawfastqc
├── 02_trimmed
├── 03_trimmedfastqc
├── 04_alignments
├── 05_lowcoverageBed
├── 06_gatkresults
└── 07_assembly
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

# To index reference fasta file
nextflow run 00_indexbwa.nf --ref "resources/sequence.fasta" --outdir="results/00_indexes"

# Running fastqc on multiple files
nextflow run 01_fastqc.nf --raw 'data/*R{1,2}_001.fastq.gz' --outdir="results/01_rawfastqc"

# Trimming the data (optional) based on fastqc report
nextflow run 02_trimming.nf --raw "data/*R{1,2}_001.fastq.gz" --run "ibex"  

# Fastqc on trimmed reads
nextflow run 01_fastqc.nf --raw "results/02_trimmed/*R{1,2}_P.fastq.gz" --outdir "results/03_trimmedfastqc" 

# Aligning reads
nextflow run 04_alignbwa.nf --raw "results/02_trimmed/*R{1,2}_P.fastq.gz" 

# Identifying regions of low coverage to mask in the final assembly
nextflow run 05_lowCoverage.nf --raw 'results/04_alignments/*bam'

# Running GATK variant calling
nextflow run 06_gatk.nf

# Making final assembly
nextflow run 07_assembly.nf
```
Save the above file as `nextcov.sbatch` and submit job on ibex using command

```
sbatch nextcov.sbatch
```
