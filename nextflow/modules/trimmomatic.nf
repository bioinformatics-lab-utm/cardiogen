process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/04_trimmomatic", mode: 'copy'

    container 'staphb/trimmomatic:0.39'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_paired.fastq.gz"),   emit: paired_reads
    tuple val(meta), path("*_unpaired.fastq.gz"), emit: unpaired_reads
    tuple val(meta), path("*_trimmomatic.log"),   emit: log
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        trimmomatic SE \\
            -threads $task.cpus \\
            $args \\
            ${reads[0]} \\
            ${prefix}_trimmed.fastq.gz \\
            ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \\
            LEADING:3 \\
            TRAILING:3 \\
            SLIDINGWINDOW:4:15 \\
            MINLEN:36 \\
            2> ${prefix}_trimmomatic.log

        # Rename for consistency with paired-end output
        mv ${prefix}_trimmed.fastq.gz ${prefix}_paired.fastq.gz
        touch ${prefix}_unpaired.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version 2>&1 | sed 's/^Trimmomatic //')
        END_VERSIONS
        """
    } else {
        """
        trimmomatic PE \\
            -threads $task.cpus \\
            $args \\
            ${reads[0]} \\
            ${reads[1]} \\
            ${prefix}_1_paired.fastq.gz \\
            ${prefix}_1_unpaired.fastq.gz \\
            ${prefix}_2_paired.fastq.gz \\
            ${prefix}_2_unpaired.fastq.gz \\
            ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \\
            LEADING:3 \\
            TRAILING:3 \\
            SLIDINGWINDOW:4:15 \\
            MINLEN:36 \\
            2> ${prefix}_trimmomatic.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version 2>&1 | sed 's/^Trimmomatic //')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}_paired.fastq.gz
        touch ${prefix}_unpaired.fastq.gz
        touch ${prefix}_trimmomatic.log
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version 2>&1 | sed 's/^Trimmomatic //')
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_1_paired.fastq.gz
        touch ${prefix}_1_unpaired.fastq.gz
        touch ${prefix}_2_paired.fastq.gz
        touch ${prefix}_2_unpaired.fastq.gz
        touch ${prefix}_trimmomatic.log
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version 2>&1 | sed 's/^Trimmomatic //')
        END_VERSIONS
        """
    }
}