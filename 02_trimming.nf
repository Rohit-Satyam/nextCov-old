//process_exercise_input_answer.nf

nextflow.enable.dsl=2
params.raw = "data/*{1,2}.fastq.gz"
params.adapter = "$baseDir/resources/universal_adapter.fasta"
params.run="ibex"
process TRIM {
	publishDir "results/02_trimmed", mode: 'copy'
	input:
	tuple val(sample_id), file(reads)
	output:
          tuple val(sample_id), file(fq_1_paired), file(fq_2_paired)
	script:
    fq_1_paired = sample_id + '_R1_P.fastq.gz'
    fq_1_unpaired = sample_id + '_R1_UP.fastq.gz'
    fq_2_paired = sample_id + '_R2_P.fastq.gz'
    fq_2_unpaired = sample_id + '_R2_UP.fastq.gz'
	if ("${params.run}" == "local")
	"""
	trimmomatic PE -threads 20 ${reads[0]} ${reads[1]} \
	$fq_1_paired \
	$fq_1_unpaired \
	$fq_2_paired \
	$fq_2_unpaired \
	ILLUMINACLIP:${params.adapter}:2:30:10 \
	MINLEN:30
	"""
	else if ("${params.run}" == "ibex")
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

reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )
workflow {
  TRIM(reads_ch)
}
