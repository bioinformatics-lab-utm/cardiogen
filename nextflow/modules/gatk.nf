process GATK_HAPLOTYPECALLER_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    // No publishDir: the unfiltered calls stay in work/ only. The published GATK results
    // are the hard-filtered ones (GATK_HARDFILTER_HG38_BWAMEM).

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path exome_bed
    path exome_bed_tbi

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    def reference = "/reference/hg38/hg38.fa"
    def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
    def intervals = params.exome_mode && exome_bed.name != 'NO_FILE' ? "-L ${exome_bed} --interval-padding ${params.exome_padding}" : ''
    """
    gatk ${avail_mem} HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        -O ${prefix}.vcf.gz \\
        --native-pair-hmm-threads $task.cpus \\
        ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
    END_VERSIONS
    """
}

// Hard filtering of HaplotypeCaller calls (GATK Best Practices recommended thresholds).
// A single-sample panel has too few variants for VQSR, so hard filters are used instead.
// SNPs and INDELs (+ MIXED) are filtered separately, because the thresholds differ, then merged.
// Kept as a separate process so changing the filters does not re-run HaplotypeCaller (-resume).
// Outputs:
//   *.hardfiltered.vcf.gz - all calls, failing ones marked in the FILTER column (QD2, FS60, ...);
//                           intermediate, kept in work/ only
//   *.pass.vcf.gz         - only FILTER=PASS calls; the only GATK result published to outdir
process GATK_HARDFILTER_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_low'

    // Only the PASS callset is published; *.hardfiltered.vcf.gz stays in work/
    publishDir "${params.outdir}/05_variant_calling/${meta.run}/gatk/bwamem/${meta.qc_tool}/hg38/hard_filtered", mode: 'copy', pattern: '*.pass.vcf.gz*'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.hardfiltered.vcf.gz"), path("*.hardfiltered.vcf.gz.tbi"), emit: filtered
    tuple val(meta), path("*.pass.vcf.gz"),         path("*.pass.vcf.gz.tbi"),         emit: pass
    path "versions.yml",                                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    def reference = "/reference/hg38/hg38.fa"
    def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
    """
    # SNPs
    gatk ${avail_mem} SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        -O snps.vcf.gz

    gatk ${avail_mem} VariantFiltration \\
        -R ${reference} \\
        -V snps.vcf.gz \\
        -filter "QD < 2.0"              --filter-name "QD2" \\
        -filter "QUAL < 30.0"           --filter-name "QUAL30" \\
        -filter "SOR > 3.0"             --filter-name "SOR3" \\
        -filter "FS > 60.0"             --filter-name "FS60" \\
        -filter "MQ < 40.0"             --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -O snps.filtered.vcf.gz

    # INDELs (+ MIXED sites)
    gatk ${avail_mem} SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        --select-type-to-include MIXED \\
        -O indels.vcf.gz

    gatk ${avail_mem} VariantFiltration \\
        -R ${reference} \\
        -V indels.vcf.gz \\
        -filter "QD < 2.0"               --filter-name "QD2" \\
        -filter "QUAL < 30.0"            --filter-name "QUAL30" \\
        -filter "FS > 200.0"             --filter-name "FS200" \\
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        -O indels.filtered.vcf.gz

    # Merge back into one sorted, indexed VCF
    gatk ${avail_mem} MergeVcfs \\
        -I snps.filtered.vcf.gz \\
        -I indels.filtered.vcf.gz \\
        -O ${prefix}.hardfiltered.vcf.gz

    # PASS-only callset
    gatk ${avail_mem} SelectVariants \\
        -R ${reference} \\
        -V ${prefix}.hardfiltered.vcf.gz \\
        --exclude-filtered \\
        -O ${prefix}.pass.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
    END_VERSIONS
    """
}
