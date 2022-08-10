nextflow.enable.dsl=2
params.raw = "results/04_alignments/*.bam"
params.exclude = "$baseDir/resources/exclude.bed"
params.outdir="results/05_lowcoverageBed"

process COV{
	publishDir "$params.outdir", mode: 'copy'
	input:
	tuple val(sid), path(bam)
	output:
	file "*.bed"
	script:
	"""	
	covtobed -x 30 ${bam} > ${sid}_lowcoverage.bed
	subtractBed -a ${sid}_lowcoverage.bed -b ${params.exclude} | mergeBed  > ${sid}.bed
	"""
}

reads_ch = Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) }

workflow {
  COV(reads_ch)
}
