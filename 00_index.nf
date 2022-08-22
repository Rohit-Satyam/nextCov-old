nextflow.enable.dsl=2
params.ref = "resources/sequence.fasta"
params.outdir="results/00_indexes"
params.runidx="bwa"
params.help = false    

// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 0: Indexing
=============================================
Usage:
	nextflow run 00_index.nf --ref ${params.ref} --outdir ${params.outdir} --runidx ${params.runidx}
Input:
	* --ref: Path of reference file. Defult [${params.ref}]
	* --outdir: name of output directory. Default [${params.outdir}]
	* --runidx: Name of tool to run indexing. Valid values are "bwa" and "dragmap". Default [${params.runidx}]
"""

exit 0
}


process DRAGMAPINDEX{
	publishDir "$params.outdir", mode: 'copy'
	input:
	path (fasta)
	output:
	path "*"
	script:
	"""
	mkdir dragmapidx
	cp $fasta dragmapidx/
	samtools faidx dragmapidx/$fasta 
	gatk CreateSequenceDictionary -R dragmapidx/$fasta	
	dragen-os --build-hash-table true --ht-reference dragmapidx/$fasta  --output-directory dragmapidx --ht-num-threads 20
	gatk ComposeSTRTableFile -R dragmapidx/$fasta -O dragmapidx/str_table.tsv
	"""
}


process BWAINDEX{
	publishDir "$params.outdir", mode: 'copy'
	input:
	path (fasta)
	output:
	path "*"
	script:
	"""
	mkdir bwaidx
	cp $fasta bwaidx/
	bwa index bwaidx/$fasta
	gatk CreateSequenceDictionary -R bwaidx/$fasta
	samtools faidx bwaidx/$fasta
	"""
}

fa_ch=Channel.fromPath(params.ref, checkIfExists: true)



workflow {
if ("${params.runidx}" == "bwa"){
 BWAINDEX(fa_ch)
} else if ("${params.runidx}" == "dragmap"){
  DRAGMAPINDEX(fa_ch)
  } else {
  exit 1,  "Invalid argument passed to --runidx"
  }
}


