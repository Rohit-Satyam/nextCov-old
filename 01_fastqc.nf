//process_exercise_input_answer.nf

nextflow.enable.dsl=2
params.raw = "data/*{1,2}.fastq.gz"
params.outdir="results/01_rawfastqc"
process FASTQC {
	publishDir "$params.outdir", mode: 'copy'
	input:
	tuple val(sample_id), path(reads)
	output:
	file("*.{html,zip}")
	script:
	"""
	fastqc -t 10 ${reads[0]} ${reads[1]}
	"""
}


reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )

workflow {
  FASTQC(reads_ch)
}
