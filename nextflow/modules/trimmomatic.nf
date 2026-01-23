process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/03_trimmomatic", mode: 'copy'

    container 'staphb/trimmomatic:0.39'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_*_paired.fastq.gz"),   emit: paired_reads
    tuple val(meta), path("*_*_unpaired.fastq.gz"), emit: unpaired_reads
    tuple val(meta), path("*.log"),                 emit: log
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_path = "/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
    
    if (meta.single_end) {
        """
        trimmomatic SE \\
            -threads $task.cpus \\
            ${reads[0]} \\
            ${prefix}_trimmed.fastq.gz \\
            ILLUMINACLIP:${adapter_path}:2:30:10 \\
            $args \\
            2> ${prefix}_trimmomatic.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version)
        END_VERSIONS
        """
    } else {
        """
        trimmomatic PE \\
            -threads $task.cpus \\
            ${reads[0]} \\
            ${reads[1]} \\
            ${prefix}_1_paired.fastq.gz \\
            ${prefix}_1_unpaired.fastq.gz \\
            ${prefix}_2_paired.fastq.gz \\
            ${prefix}_2_unpaired.fastq.gz \\
            ILLUMINACLIP:${adapter_path}:2:30:10 \\
            $args \\
            2> ${prefix}_trimmomatic.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version)
        END_VERSIONS
        """
    }
}
