process FASTQC_TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/04_trimmomatic/fastqc", mode: 'copy'

    container 'staphb/fastqc:0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        fastqc \\
            $args \\
            --threads $task.cpus \\
            ${reads[0]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v\\([0-9\\.]*\\).*/\\1/' )
        END_VERSIONS
        """
    } else {
        """
        fastqc \\
            $args \\
            --threads $task.cpus \\
            ${reads[0]} \\
            ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v\\([0-9\\.]*\\).*/\\1/' )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fastqc.html
    touch ${prefix}_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v\\([0-9\\.]*\\).*/\\1/' )
    END_VERSIONS
    """
}