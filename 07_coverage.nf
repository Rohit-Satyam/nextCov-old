nextflow.enable.dsl=2
params.raw = "results/05_markduplicates/*.bam"
params.outdir="results/07_coverage"
params.cpus = 20
params.help = false

// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 7: BAM Coverage
=============================================
Usage:
	nextflow run 07_coverage.nf --raw ${params.raw} --outdir ${params.outdir} --cpus ${params.cpu}
Input:
	* --raw: Path of .bam files from MarkDuplicate Step. Defult [${params.raw}]
	* --outdir: name of output directory. Default [${params.outdir}]
	* --cpus No. of threads to use. Default [${params.cpus}]
"""

exit 0
}


process SAMCOV{
	publishDir "$params.outdir", mode: 'copy'
	cpus params.cpus
	input:
	tuple val(sid), path(bam)
	output:
	file "*"
	script:
	"""	
	samtools coverage ${bam} > ${sid}.coverage.tsv
	samtools depth -@ ${task.cpus} -s -d 0 -H ${bam} > ${sid}.depth.tsv
	mosdepth -t ${task.cpus} -x ${sid} ${bam.toRealPath()}
	"""
}

process MULTIQC {
	publishDir "$params.outdir", mode: 'copy'
	input:
	file x
	output:
	file "multiqc_report.html"
	file "multiqc_data"

    script:
    """
    multiqc -f $x
    """
}



reads_ch = Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) }

workflow {
  SAMCOV(reads_ch) | collect | MULTIQC
}
