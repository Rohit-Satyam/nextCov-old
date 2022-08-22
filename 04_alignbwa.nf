nextflow.enable.dsl=2
params.raw = "results/02_trimmed/*R{1,2}_P.fastq.gz"
params.ref = "$baseDir/results/00_indexes/bwaidx/sequence.fasta"
params.outdir="results/04_alignments"
params.help = false
params.cpus = 8


// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 4: Alignment
=============================================
Usage:
	nextflow run 04_alignbwa.nf --raw ${params.raw} --ref ${params.ref} --outdir ${params.outdir} --cpus ${params.cpu}
Input:
	* --raw: Path of fastq files. Defult [${params.raw}]
	* --outdir: name of output directory. Default [${params.outdir}]
	* --ref Absolute path to reference FASTA file. Default [${params.ref}]
	* --cpus No. of threads to use. Default [${params.cpus}]
"""

exit 0
}


process BWAMEM{
	publishDir "$params.outdir", mode: 'copy'
	cpus params.cpus
	input:
	tuple val(sample_id), path(reads)
	output:
	file "*.sorted.bam"
	file "*.bai" 
	shell:
	'''
	id=$(zcat !{reads[0]} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//')
	bwa mem -M -R "$(echo "@RG\\tID:${id}\\tSM:!{sample_id}\\tPL:ILLUMINA")" -t !{task.cpus} !{params.ref} !{reads[0]} !{reads[1]} | samtools sort -@8 -o !{sample_id}.sorted.bam - 
	samtools index -@8 !{sample_id}.sorted.bam
	'''
}

reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )

workflow {
  BWAMEM(reads_ch)
}
