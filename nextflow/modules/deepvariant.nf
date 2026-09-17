process DEEPVARIANT_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    // No publishDir: the raw VCF / gVCF (which contain RefCall/LowQual/NoCall records) stay in
    // work/ only. The published DeepVariant result is the PASS-only VCF (DV_PASS_FILTER_HG38_BWAMEM).

    container 'google/deepvariant:1.6.1'

    input:
    tuple val(meta), path(bam), path(bai)
    path exome_bed
    path exome_bed_tbi

    output:
    tuple val(meta), path("*[!g].vcf.gz"),     emit: vcf
    tuple val(meta), path("*[!g].vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.g.vcf.gz"),       emit: gvcf
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    def reference = "/reference/hg38/hg38.fa"
    def model_type = params.seq_platform ?: 'WGS'
    def regions = params.exome_mode && exome_bed.name != 'NO_FILE' ? "--regions exome_regions.bed" : ''
    """
    export TMPDIR=\$(pwd)/tmp && mkdir -p \$TMPDIR

    # DeepVariant requires uncompressed BED files for --regions
    if [ -f "${exome_bed}" ] && [ "${exome_bed}" != "NO_FILE" ]; then
        gunzip -c ${exome_bed} > exome_regions.bed
    fi

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=${model_type} \\
        --ref=${reference} \\
        --reads=${bam} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        --num_shards=$task.cpus \\
        ${regions} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version 2>&1) | sed 's/^.*DeepVariant version //; s/ .*\$//')
    END_VERSIONS
    """
}

// Keeps only DeepVariant FILTER=PASS records. Dropped:
//   RefCall - candidate site genotyped as reference (GT 0/0 or ./.), not a variant
//   LowQual - variant confidence below the calling threshold
//   NoCall  - no genotype could be determined (GT ./.)
// Kept as a separate process so it does not re-run DeepVariant (-resume).
process DV_PASS_FILTER_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_low'

    publishDir "${params.outdir}/05_variant_calling/${meta.run}/deepvariant/bwamem/${meta.qc_tool}/hg38", mode: 'copy', pattern: '*.pass.vcf.gz*'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.pass.vcf.gz"), path("*.pass.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    """
    bcftools view \\
        --apply-filters PASS \\
        --output-type z \\
        --output ${prefix}.pass.vcf.gz \\
        ${vcf}
    bcftools index --tbi ${prefix}.pass.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}