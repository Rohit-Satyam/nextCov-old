nextflow.enable.dsl=2
params.raw = "results/06_gatkresults/vcfs/*.vcf.gz"
params.ref = "$baseDir/results/00_indexes/bwaidx/sequence.fasta"
params.outdir="results/07_assembly"
params.bed="$baseDir/results/05_lowcoverageBed"


process COV{
	publishDir "$params.outdir", mode: 'copy'
	input:
	tuple val(sid), path(vcf)
	output:
	path "*.fa"
	script:
	file= sid + '.fa'
	"""
	bcftools consensus -f ${params.ref} -m ${params.bed}/${sid}.bed ${vcf.toRealPath()} | sed 's/NC_045512.2/${sid}/' > $file
	"""
}

reads_ch = Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) }

workflow {
  COV(reads_ch)
}

