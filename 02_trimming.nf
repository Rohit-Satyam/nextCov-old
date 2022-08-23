nextflow.enable.dsl=2

params.raw = "data/*R{1,2}_001.fastq.gz"
params.adapter = "$baseDir/resources/universal_adapter.fasta"
params.mode="local"
params.outdir="results/02_trimmed"
params.help = false


// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 2: Trimming
=============================================
Usage:
	nextflow run 02_trimming.nf --raw ${params.raw} --adapter ${params.adapter} --mode ${params.mode} --outdir ${params.outdir}
Input:
	* --raw: Path of fastq files. Defult [${params.raw}]
	* --outdir: name of output directory. Default [${params.outdir}]
	* --mode Run mode to define. If running on your system choose "local". For cluster, use "ibex", Default [${params.mode}]
	* --adapter Specify the path of adapter file. Default [${params.adapter}]
"""

exit 0
}


process TRIM {
	publishDir "$params.outdir", mode: 'copy'
	input:
	tuple val(sample_id), file(reads)
	output:
         tuple file(fq_1_paired), file(fq_2_paired)

	script:
	
    fq_1_paired = sample_id + '_R1_P.fastq.gz'
    fq_1_unpaired = sample_id + '_R1_UP.fastq.gz'
    fq_2_paired = sample_id + '_R2_P.fastq.gz'
    fq_2_unpaired = sample_id + '_R2_UP.fastq.gz'
if ("${params.mode}" == "local")
	"""
	trimmomatic PE \
	-threads 20 \
	${reads[0]} \
	${reads[1]}\
	$fq_1_paired \
	$fq_1_unpaired \
	$fq_2_paired \
	$fq_2_unpaired \
	ILLUMINACLIP:${params.adapter}:2:30:10 \
	MINLEN:30
	"""
else if ("${params.mode}" == "ibex")
	"""
java -jar $TRIMMOMATIC_JAR PE -threads 20 ${reads[0]} ${reads[1]} \
        $fq_1_paired \
        $fq_1_unpaired \
        $fq_2_paired \
        $fq_2_unpaired \
        ILLUMINACLIP:${params.adapter}:2:30:10 \
        MINLEN:30
	"""

}

workflow {
  reads_ch=channel.fromFilePairs(params.raw, checkIfExists: true ) 
TRIM(reads_ch)
}
