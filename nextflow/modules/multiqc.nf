process MULTIQC {
    label 'process_single'

    publishDir "${params.outdir}/01_QC/multiqc", mode: 'copy'

    container 'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'multiqc'
    
    """
    multiqc \\
        $args \\
        --filename ${prefix}_report.html \\
        --force \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^.*multiqc, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: 'multiqc'
    """
    mkdir ${prefix}_data
    touch ${prefix}_report.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^.*multiqc, version //')
    END_VERSIONS
    """
}