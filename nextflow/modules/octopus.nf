process OCTOPUS_HG37_BOWTIE2 {
    tag "$meta.id - $meta.qc_tool - hg37 - bowtie2"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/octopus/bowtie2/${meta.qc_tool}/hg37", mode: 'copy'

    container 'dancooke/octopus:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg37_bowtie2"
    def reference = "/reference/hg37/hg19.fa"
    """
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | grep -oP 'octopus \\K[0-9.]+')
    END_VERSIONS
    """
}

process OCTOPUS_HG38_BOWTIE2 {
    tag "$meta.id - $meta.qc_tool - hg38 - bowtie2"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/octopus/bowtie2/${meta.qc_tool}/hg38", mode: 'copy'

    container 'dancooke/octopus:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bowtie2"
    def reference = "/reference/hg38/hg38.fa"
    """
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | grep -oP 'octopus \\K[0-9.]+')
    END_VERSIONS
    """
}

process OCTOPUS_T2T_BOWTIE2 {
    tag "$meta.id - $meta.qc_tool - t2t - bowtie2"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/octopus/bowtie2/${meta.qc_tool}/t2t", mode: 'copy'

    container 'dancooke/octopus:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_t2t_bowtie2"
    def reference = "/reference/t2t/hs1.fa"
    """
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | grep -oP 'octopus \\K[0-9.]+')
    END_VERSIONS
    """
}

process OCTOPUS_HG37_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg37 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/octopus/bwamem/${meta.qc_tool}/hg37", mode: 'copy'

    container 'dancooke/octopus:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg37_bwamem"
    def reference = "/reference/hg37/hg19.fa"
    """
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | grep -oP 'octopus \\K[0-9.]+')
    END_VERSIONS
    """
}

process OCTOPUS_HG38_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/octopus/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'dancooke/octopus:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    def reference = "/reference/hg38/hg38.fa"
    """
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | grep -oP 'octopus \\K[0-9.]+')
    END_VERSIONS
    """
}

process OCTOPUS_T2T_BWAMEM {
    tag "$meta.id - $meta.qc_tool - t2t - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/octopus/bwamem/${meta.qc_tool}/t2t", mode: 'copy'

    container 'dancooke/octopus:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_t2t_bwamem"
    def reference = "/reference/t2t/hs1.fa"
    """
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | grep -oP 'octopus \\K[0-9.]+')
    END_VERSIONS
    """
}
