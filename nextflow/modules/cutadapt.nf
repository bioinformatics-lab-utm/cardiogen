process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/03_cutadapt", mode: 'copy'

    container 'quay.io/biocontainers/cutadapt:4.4--py39hf95cd2a_1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*_cutadapt.log"),     emit: log
    path "versions.yml",                         emit: versions

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
            --quality-cutoff 20 \\
            --minimum-length 36 \\
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
            --output ${prefix}_trimmed.fastq.gz \\
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
            --quality-cutoff 20 \\
            --minimum-length 36 \\
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \\
            --output ${prefix}_1_trimmed.fastq.gz \\
            --paired-output ${prefix}_2_trimmed.fastq.gz \\
            ${reads[0]} \\
            ${reads[1]} \\
            > ${prefix}_cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}_trimmed.fastq.gz
        touch ${prefix}_cutadapt.log
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_1_trimmed.fastq.gz
        touch ${prefix}_2_trimmed.fastq.gz
        touch ${prefix}_cutadapt.log
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}