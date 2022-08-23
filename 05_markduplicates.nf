nextflow.enable.dsl=2
params.raw = "results/04_alignments/*.bam"
params.outdir="results/05_markduplicates"
params.jobs = 5
params.help = false
// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 5: MarkDuplicates (Optical/PCR)
=============================================
Usage:
	nextflow run 05_markduplicates.nf --raw ${params.raw} --outdir ${params.outdir} --jobs ${params.jobs}
Input:
	* --raw: Path to .bam files. Defult [${params.raw}]
	* --outdir: Name of output directory. Default [${params.outdir}]
	* --jobs: Number of jobs to be run in parallel. Decrease this when getting error. Default [$params.jobs]
"""

exit 0
}


process MARKDUP{
publishDir "$params.outdir", mode: 'copy'
maxForks params.jobs
	input:
	tuple val(sid), path(bam)
	output:
	file "*"
	script:
	"""	
	gatk MarkDuplicatesSpark -I ${bam} -O ${sid}.dedup.bam -M ${sid}_markdup_metrics.txt --tmp-dir . 
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


//reads_ch = Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) }

workflow {
 Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) } | MARKDUP | collect | MULTIQC
}

