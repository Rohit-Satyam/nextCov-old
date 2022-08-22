nextflow.enable.dsl=2
params.raw = "data/*.fastq.gz"
params.outdir="results/01_rawfastqc"
params.help = false
// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 1: FASTQC
=============================================
Usage:
	nextflow run 01_fastqc.nf --raw ${params.raw} --outdir ${params.outdir}
Input:
	* --raw: Path of fastq files. Defult [${params.raw}]
	* --outdir: name of output directory. Default [${params.outdir}]
"""

exit 0
}


process FASTQC {
	publishDir "$params.outdir", mode: 'copy'
	input:
	path reads
	output:
	file "*"
	script:
	"""
	fastqc -t 10 ${reads}
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



workflow {
 Channel.fromPath(params.raw, checkIfExists: true )| FASTQC | collect | MULTIQC
  
}


