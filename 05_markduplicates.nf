nextflow.enable.dsl=2
params.raw = "results/04_alignments/*.bam"
params.outdir="results/05_markduplicates"

params.help = false
// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 5: MarkDuplicates (Optical/PCR)
=============================================
Usage:
	nextflow run 05_markduplicates.nf --raw ${params.raw} --outdir ${params.outdir}
Input:
	* --raw: Path to .bam files. Defult [${params.raw}]
	* --outdir: Name of output directory. Default [${params.outdir}]
"""

exit 0
}


process MARKDUP{
publishDir "$params.outdir", mode: 'copy'
	input:
	tuple val(sid), path(bam)
	output:
	path "*"
	script:
	"""	
	gatk MarkDuplicatesSpark -I ${bam} -O ${sid}.dedup.bam -M ${sid}_markdup_metrics.txt
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
  MARKDUP(reads_ch) | collect | MULTIQC
}

