process FASTP {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/02_fastp", mode: 'copy'

    container 'staphb/fastp:0.23.4'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastp.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),           emit: json
    tuple val(meta), path("*.html"),           emit: html
    path "versions.yml",                       emit: versions

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
            --json ${prefix}_fastp.json \\
            --html ${prefix}_fastp.html \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            --correction \\
            --cut_front \\
            --cut_tail \\
            --trim_poly_g \\
            --trim_poly_x \\
            --n_base_limit 5 \\
            --qualified_quality_phred 15 \\
            --unqualified_percent_limit 40 \\
            --length_required 36 \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^.*fastp //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_1_fastp.fastq.gz \\
            --out2 ${prefix}_2_fastp.fastq.gz \\
            --json ${prefix}_fastp.json \\
            --html ${prefix}_fastp.html \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            --correction \\
            --cut_front \\
            --cut_tail \\
            --trim_poly_g \\
            --trim_poly_x \\
            --n_base_limit 5 \\
            --qualified_quality_phred 15 \\
            --unqualified_percent_limit 40 \\
            --length_required 36 \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^.*fastp //; s/ .*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}_fastp.fastq.gz
        touch ${prefix}_fastp.json
        touch ${prefix}_fastp.html
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^.*fastp //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_1_fastp.fastq.gz
        touch ${prefix}_2_fastp.fastq.gz
        touch ${prefix}_fastp.json
        touch ${prefix}_fastp.html
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^.*fastp //; s/ .*\$//')
        END_VERSIONS
        """
    }
}