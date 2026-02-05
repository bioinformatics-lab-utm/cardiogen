process CUTADAPT_CLEAN {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}/05_cutadapt_trimmomatic", mode: 'copy', pattern: "*_cutadapt.log"

    container 'quay.io/biocontainers/cutadapt:4.4--py39hf95cd2a_1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_cutadapt.fastq.gz"), emit: reads
    tuple val(meta), path("*_cutadapt.log"),      emit: log
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        cutadapt \\
            --cores $task.cpus \\
            --minimum-length 36 \\
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
            --output ${prefix}_cutadapt.fastq.gz \\
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
            --cores $task.cpus \\
            --minimum-length 36 \\
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \\
            --output ${prefix}_1_cutadapt.fastq.gz \\
            --paired-output ${prefix}_2_cutadapt.fastq.gz \\
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
        touch ${prefix}_cutadapt.fastq.gz
        touch ${prefix}_cutadapt.log
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_1_cutadapt.fastq.gz
        touch ${prefix}_2_cutadapt.fastq.gz
        touch ${prefix}_cutadapt.log
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}