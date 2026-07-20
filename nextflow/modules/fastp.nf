process FASTP {
    tag "$meta.run - $meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/02_fastp/${meta.run}", mode: 'copy'
    container 'staphb/fastp:0.23.4'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_fastp_{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),                 emit: json
    tuple val(meta), path("*.html"),                 emit: html
    path "versions.yml",                             emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}_fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}_fastp.json \\
            --html ${prefix}_fastp.html \\
            $args
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | head -n1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_fastp_1.fastq.gz \\
            --out2 ${prefix}_fastp_2.fastq.gz \\
            --fix_mgi_id \\
            --thread $task.cpus \\
            --json ${prefix}_fastp.json \\
            --html ${prefix}_fastp.html \\
            $args
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | head -n1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
}