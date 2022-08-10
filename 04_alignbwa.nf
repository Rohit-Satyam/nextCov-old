nextflow.enable.dsl=2
params.raw = "data/*{1,2}.fastq.gz"
params.ref = "$baseDir/results/00_indexes/bwaidx/sequence.fasta"
params.outdir="results/04_alignments"

process BWAMEM{
	publishDir "$params.outdir", mode: 'copy'
	input:
	tuple val(sample_id), path(reads)
	output:
	file "*.sorted.bam"
	file "*.bai" 
	shell:
	'''
	id=$(zcat !{reads[0]} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//')
	echo "$id"
	bwa mem -M -R "$(echo "@RG\\tID:${id}\\tSM:!{sample_id}\\tPL:ILLUMINA")" -t 8 !{params.ref} !{reads[0]} !{reads[1]} | samtools sort -@8 -o !{sample_id}.sorted.bam -
	samtools index -@8 !{sample_id}.sorted.bam
	'''
}

reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )

workflow {
  BWAMEM(reads_ch)
}
