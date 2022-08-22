nextflow.enable.dsl=2
params.raw = "results/05_markduplicates/*.bam"
params.exclude = "$baseDir/resources/exclude.bed"
params.outdir="results/06_lowcoverageBed"

params.help = false
// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 5: MarkDuplicates (Optical/PCR)
=============================================
Usage:
	nextflow run 06_lowCoverage.nf --raw ${params.raw} --outdir ${params.outdir} --exclude ${params.exclude}
Input:
	* --raw: Path to .bam files. Defult [${params.raw}]
	* --outdir: Name of output directory. Default [${params.outdir}]
	* --exclude: BED file for the known regions to be excluded from resulting low-coverage BED file (eg: Poly-A Tail).
"""

exit 0
}


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
