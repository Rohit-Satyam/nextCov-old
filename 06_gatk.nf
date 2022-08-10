nextflow.enable.dsl=2
params.raw = "results/04_alignments/*.bam"
params.ref = "$baseDir/results/00_indexes/bwaidx/sequence.fasta"
params.outdir="results/06_gatkresults"

process COV{
    publishDir "$params.outdir/bams", pattern: '*.bam', mode: 'copy'
	publishDir "$params.outdir/intermediate", pattern: '*.vcf', mode: 'copy'
    publishDir "$params.outdir/vcfs", pattern: '*.vcf.gz*', mode: 'copy'

	input:
	tuple val(sid), path(bam)
	output:
	file "*.bam"
	file "*.vcf"
	file "*.vcf.gz*"
	script:
	"""	
	gatk MarkDuplicatesSpark -I ${bam} -O ${sid}.dedup.bam
	gatk HaplotypeCaller --sample-name ${sid} -ploidy 1 --native-pair-hmm-threads 8 -R ${params.ref} -I ${sid}.dedup.bam -O ${sid}.raw.vcf
        gatk VariantFiltration -R ${params.ref} -V ${sid}.raw.vcf --filter-expression \"QD<2.0 || FS > 60.0 || SOR > 4.0 || ReadPosRankSum < -8.0\"  -filter-name  \"good_filter\" -O ${sid}.filtered_variants.vcf
        gatk SelectVariants -V ${sid}.filtered_variants.vcf --select-type-to-include SNP   -O ${sid}_known.vcf
        gatk BaseRecalibrator -R ${params.ref} -I ${sid}.dedup.bam --known-sites  ${sid}_known.vcf  -O ${sid}_recal_data.table
        gatk ApplyBQSR -R ${params.ref} -I ${sid}.dedup.bam --bqsr-recal-file ${sid}_recal_data.table -O ${sid}_recal.bam
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
	"""
}

reads_ch = Channel.fromPath(params.raw).map { file -> tuple(file.simpleName, file) }

workflow {
  COV(reads_ch)
}
