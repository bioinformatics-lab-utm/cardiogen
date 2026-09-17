process DELLY_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/05_variant_calling/${meta.run}/delly/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bcf"),     emit: vcf
    tuple val(meta), path("*.bcf.csi"), emit: tbi
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    def reference = "/reference/hg38/hg38.fa"
    // Note: Delly doesn't have native exome support (no --regions or --intervals option).
    // For exome data, Delly will naturally focus on regions with coverage (exome targets).
    // Optional: use -x to exclude problematic regions (centromeres, telomeres) or post-filter VCF.
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
    END_VERSIONS
    """
}