nextflow.enable.dsl=2
params.raw = "results/08_gatkresults/vcfs/*.vcf.gz"
params.ref = "$baseDir/results/00_indexes/bwaidx/sequence.fasta"
params.outdir="results/09_assembly"
params.bed="$baseDir/results/05_lowcoverageBed"

params.help = false
// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 9: Substitution Based Assembly
=============================================
Usage:
	nextflow run 09_assembly.nf --raw ${params.raw} --outdir ${params.outdir} --bed ${params.bed} --ref ${params.ref}
Input:
	* --raw: Path of VCF files. Defult [${params.raw}]
	* --ref Absolute path of Reference FASTAfile. Default [${params.ref}]
	* --outdir: Name of output directory. Default [${params.outdir}]
	* --bed Absolute path of directory where BED files containing regions to be Hard Masked (with 'N') are present.
	Default [${params.bed}]
"""

exit 0
}

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

