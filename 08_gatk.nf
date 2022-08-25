nextflow.enable.dsl=2
params.raw = "results/05_markduplicates/*.bam"
params.ref = "$baseDir/results/00_indexes/bwaidx/sequence.fasta"
params.outdir="results/08_gatkresults"
params.help = false



// Help Section

if( params.help ) {
log.info """
nextCov@KAUST Step 8: GATK Variant Calling
=============================================
Usage:
	nextflow run 08_gatk.nf --raw ${params.raw} --ref ${params.ref} --outdir ${params.outdir}
Input:
	* --raw: Path of fastq files. Defult [${params.raw}]
	* --outdir: name of output directory. Default [${params.outdir}]
	* --ref Absolute path of reference FASTA file. Default [${params.ref}]
"""

exit 0
}

process GATK{
publishDir "$params.outdir/intermediate", pattern: '*.vcf', mode: 'copy'
publishDir "$params.outdir/vcfs", pattern: '*.vcf.gz*', mode: 'copy'
publishDir "$params.outdir/vcfs", pattern: '*.txt*', mode: 'copy'
publishDir "$params.outdir/bqsrfiles",pattern: "*.table", mode: 'copy'
publishDir "$params.outdir/bqsrfiles",pattern: "*.pdf", mode: 'copy'

	input:
	tuple val(sid), path(bam)
	output:
	file "*.vcf"
	file "*.vcf.gz*"
	file "*.pdf"
	file "*.table"
	file "*.txt"
	script:
	"""	
	gatk HaplotypeCaller --sample-name ${sid} -ploidy 1 --native-pair-hmm-threads 8 -R ${params.ref} -I ${bam} -O ${sid}.raw.vcf
        gatk VariantFiltration -R ${params.ref} -V ${sid}.raw.vcf --filter-expression \"QD<2.0 || FS > 60.0 || SOR > 4.0 || ReadPosRankSum < -8.0\"  -filter-name  \"good_filter\" -O ${sid}.filtered_variants.vcf
        gatk SelectVariants -V ${sid}.filtered_variants.vcf --select-type-to-include SNP   -O ${sid}_known.vcf
        gatk BaseRecalibrator -R ${params.ref} -I ${sid}.dedup.bam --known-sites  ${sid}_known.vcf  -O ${sid}_beforerecal_data.table
	gatk ApplyBQSR -R ${params.ref} -I ${sid}.dedup.bam --bqsr-recal-file ${sid}_beforerecal_data.table -O ${sid}_recal.bam
	gatk BaseRecalibrator -R ${params.ref} -I ${sid}_recal.bam --known-sites ${sid}_known.vcf -O ${sid}_afterrecal_data.table
	gatk AnalyzeCovariates -before ${sid}_beforerecal_data.table -after ${sid}_afterrecal_data.table -plots ${sid}_AnalyzeCovariates.pdf 
        gatk HaplotypeCaller --sample-name ${sid} -ploidy 1 --native-pair-hmm-threads 8 -R ${params.ref} -I ${sid}_recal.bam -O ${sid}_recal.raw.vcf
        gatk SelectVariants -V ${sid}_recal.raw.vcf --select-type-to-include SNP  -O ${sid}_snps.vcf
        gatk SelectVariants -V ${sid}_recal.raw.vcf --select-type-to-include INDEL  -O ${sid}_indels.vcf
        gatk VariantFiltration -R ${params.ref} -V ${sid}_snps.vcf --filter-expression \"QD<2.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"  -filter-name  \"bad_SNP\" -O ${sid}_filtered.snps.vcf
        gatk VariantFiltration -R ${params.ref} -V ${sid}_indels.vcf --filter-expression \"QD<2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0\"  -filter-name  \"bad_INDEL\" -O ${sid}_filtered.indels.vcf
        gatk SelectVariants --exclude-filtered true -V ${sid}_filtered.indels.vcf -O ${sid}_good_indels.vcf
        gatk SelectVariants --exclude-filtered true -V ${sid}_filtered.snps.vcf -O ${sid}_good_snps.vcf
        gatk MergeVcfs -I ${sid}_good_snps.vcf -I ${sid}_good_indels.vcf -O ${sid}.mergeSNPsIndels.vcf
	bgzip ${sid}.mergeSNPsIndels.vcf
	bcftools index ${sid}.mergeSNPsIndels.vcf.gz
	bcftools stats ${sid}.mergeSNPsIndels.vcf.gz > ${sid}.bcftoolstats.txt
	"""
}


workflow {
 Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) } | GATK 
}
