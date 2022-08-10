nextflow.enable.dsl=2
params.ref = "resources/sequence.fasta"
params.outdir="results/00_indexes"

process BWAINDEX{
	publishDir "$params.outdir", mode: 'copy'
	input:
	path (fasta)
	output:
	path "*"
	script:
	"""
	mkdir bwaidx
	cp $fasta bwaidx/
	bwa index bwaidx/$fasta
	gatk CreateSequenceDictionary -R bwaidx/$fasta
	samtools faidx bwaidx/$fasta
	"""
}

fa_ch=Channel.fromPath(params.ref, checkIfExists: true)
workflow {
  BWAINDEX(fa_ch)
}
