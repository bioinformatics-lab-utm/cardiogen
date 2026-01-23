process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/02_cutadapt", mode: 'copy'

    container 'kfdrc/cutadapt:latest'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_cutadapt_{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"),                     emit: log
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        cutadapt \\
            $args \\
            --cores $task.cpus \\
            -o ${prefix}_cutadapt.fastq.gz \\
            ${reads[0]} \\
            > ${prefix}_cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        cutadapt \\
            $args \\
            --cores $task.cpus \\
            -o ${prefix}_cutadapt_1.fastq.gz \\
            -p ${prefix}_cutadapt_2.fastq.gz \\
            ${reads[0]} \\
            ${reads[1]} \\
            > ${prefix}_cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}
