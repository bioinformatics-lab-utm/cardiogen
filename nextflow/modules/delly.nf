process DELLY_HG37_BOWTIE2 {
    tag "$meta.id - $meta.qc_tool - hg37 - bowtie2"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bowtie2/${meta.qc_tool}/hg37", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg37_bowtie2"
    def reference = "/reference/hg37/hg19.fa"
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    # Install bcftools and htslib for conversion (Alpine Linux)
    apk add --no-cache bcftools htslib > /dev/null 2>&1

    # Convert BCF to VCF.gz
    bcftools view ${prefix}.bcf | bgzip -c > ${prefix}.vcf.gz
    
    # Index VCF
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process DELLY_HG38_BOWTIE2 {
    tag "$meta.id - $meta.qc_tool - hg38 - bowtie2"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bowtie2/${meta.qc_tool}/hg38", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bowtie2"
    def reference = "/reference/hg38/hg38.fa"
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    # Install bcftools and htslib for conversion (Alpine Linux)
    apk add --no-cache bcftools htslib > /dev/null 2>&1

    # Convert BCF to VCF.gz
    bcftools view ${prefix}.bcf | bgzip -c > ${prefix}.vcf.gz
    
    # Index VCF
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process DELLY_T2T_BOWTIE2 {
    tag "$meta.id - $meta.qc_tool - t2t - bowtie2"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bowtie2/${meta.qc_tool}/t2t", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_t2t_bowtie2"
    def reference = "/reference/t2t/hs1.fa"
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    # Install bcftools and htslib for conversion (Alpine Linux)
    apk add --no-cache bcftools htslib > /dev/null 2>&1

    # Convert BCF to VCF.gz
    bcftools view ${prefix}.bcf | bgzip -c > ${prefix}.vcf.gz
    
    # Index VCF
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process DELLY_HG37_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg37 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bwamem/${meta.qc_tool}/hg37", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg37_bwamem"
    def reference = "/reference/hg37/hg19.fa"
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    # Install bcftools and htslib for conversion (Alpine Linux)
    apk add --no-cache bcftools htslib > /dev/null 2>&1

    # Convert BCF to VCF.gz
    bcftools view ${prefix}.bcf | bgzip -c > ${prefix}.vcf.gz
    
    # Index VCF
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process DELLY_HG38_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    def reference = "/reference/hg38/hg38.fa"
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    # Install bcftools and htslib for conversion (Alpine Linux)
    apk add --no-cache bcftools htslib > /dev/null 2>&1

    # Convert BCF to VCF.gz
    bcftools view ${prefix}.bcf | bgzip -c > ${prefix}.vcf.gz
    
    # Index VCF
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process DELLY_T2T_BWAMEM {
    tag "$meta.id - $meta.qc_tool - t2t - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bwamem/${meta.qc_tool}/t2t", mode: 'copy'

    container 'dellytools/delly:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_t2t_bwamem"
    def reference = "/reference/t2t/hs1.fa"
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    # Install bcftools and htslib for conversion (Alpine Linux)
    apk add --no-cache bcftools htslib > /dev/null 2>&1

    # Convert BCF to VCF.gz
    bcftools view ${prefix}.bcf | bgzip -c > ${prefix}.vcf.gz
    
    # Index VCF
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
