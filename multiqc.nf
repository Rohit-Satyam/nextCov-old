nextflow.enable.dsl=2
params.raw = "results/08_gatkresults"
params.help = false
// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step *: MULTIQC
=============================================
Usage:
	nextflow run 01_fastqc.nf --raw ${params.raw}
Input:
	* --raw: Path of folder to summarize files.
"""

exit 0
}



process MULTIQC {
	publishDir "$params.raw", mode: 'copy'
	input:
	path x
	output:
	file "multiqc_report.html"
	file "multiqc_data"

    script:
    """
    multiqc -f $x
    """
}



workflow {
 Channel.fromPath(params.raw, checkIfExists: true )| MULTIQC
  
}


