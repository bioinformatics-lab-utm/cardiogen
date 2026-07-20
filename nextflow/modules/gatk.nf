// COMMENTED OUT - BOWTIE2 not used for exome analysis
// process GATK_HAPLOTYPECALLER_HG37_BOWTIE2 {
//     tag "$meta.id - $meta.qc_tool - hg37 - bowtie2"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/gatk/bowtie2/${meta.qc_tool}/hg37", mode: 'copy'
//
//     container 'broadinstitute/gatk:4.5.0.0'
//
//     input:
//     tuple val(meta), path(bam), path(bai)
//
//     output:
//     tuple val(meta), path("*.vcf.gz"),     emit: vcf
//     tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
//     path "versions.yml",                   emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_hg37_bowtie2"
//     def reference = "/reference/hg37/hg19.fa"
//     def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
//     """
//     gatk ${avail_mem} HaplotypeCaller \\
//         -R ${reference} \\
//         -I ${bam} \\
//         -O ${prefix}.vcf.gz \\
//         --native-pair-hmm-threads $task.cpus \\
//         $args
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
//     END_VERSIONS
//     """
// }

// COMMENTED OUT - BOWTIE2 not used for exome analysis
// process GATK_HAPLOTYPECALLER_HG38_BOWTIE2 {
//     tag "$meta.id - $meta.qc_tool - hg38 - bowtie2"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/gatk/bowtie2/${meta.qc_tool}/hg38", mode: 'copy'
//
//     container 'broadinstitute/gatk:4.5.0.0'
//
//     input:
//     tuple val(meta), path(bam), path(bai)
//
//     output:
//     tuple val(meta), path("*.vcf.gz"),     emit: vcf
//     tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
//     path "versions.yml",                   emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_hg38_bowtie2"
//     def reference = "/reference/hg38/hg38.fa"
//     def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
//     """
//     gatk ${avail_mem} HaplotypeCaller \\
//         -R ${reference} \\
//         -I ${bam} \\
//         -O ${prefix}.vcf.gz \\
//         --native-pair-hmm-threads $task.cpus \\
//         $args
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
//     END_VERSIONS
//     """
// }

// COMMENTED OUT - BOWTIE2 not used for exome analysis
// process GATK_HAPLOTYPECALLER_T2T_BOWTIE2 {
//     tag "$meta.id - $meta.qc_tool - t2t - bowtie2"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/gatk/bowtie2/${meta.qc_tool}/t2t", mode: 'copy'
//
//     container 'broadinstitute/gatk:4.5.0.0'
//
//     input:
//     tuple val(meta), path(bam), path(bai)
//
//     output:
//     tuple val(meta), path("*.vcf.gz"),     emit: vcf
//     tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
//     path "versions.yml",                   emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_t2t_bowtie2"
//     def reference = "/reference/t2t/hs1.fa"
//     def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
//     """
//     gatk ${avail_mem} HaplotypeCaller \\
//         -R ${reference} \\
//         -I ${bam} \\
//         -O ${prefix}.vcf.gz \\
//         --native-pair-hmm-threads $task.cpus \\
//         $args
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
//     END_VERSIONS
//     """
// }

// COMMENTED OUT - HG37 not needed for exome analysis
// process GATK_HAPLOTYPECALLER_HG37_BWAMEM {
//     tag "$meta.id - $meta.qc_tool - hg37 - bwamem"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/gatk/bwamem/${meta.qc_tool}/hg37", mode: 'copy'
//
//     container 'broadinstitute/gatk:4.5.0.0'
//
//     input:
//     tuple val(meta), path(bam), path(bai)
//
//     output:
//     tuple val(meta), path("*.vcf.gz"),     emit: vcf
//     tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
//     path "versions.yml",                   emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_hg37_bwamem"
//     def reference = "/reference/hg37/hg19.fa"
//     def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
//     def intervals = params.exome_mode && params.exome_bed_hg37 ? "-L ${params.exome_bed_hg37} --interval-padding ${params.exome_padding}" : ''
//     """
//     gatk ${avail_mem} HaplotypeCaller \\
//         -R ${reference} \\
//         -I ${bam} \\
//         -O ${prefix}.vcf.gz \\
//         --native-pair-hmm-threads $task.cpus \\
//         ${intervals} \\
//         $args
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
//     END_VERSIONS
//     """
// }

process GATK_HAPLOTYPECALLER_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/05_variant_calling/${meta.run}/gatk/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path exome_bed
    path exome_bed_tbi

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
    def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
    def intervals = params.exome_mode && exome_bed.name != 'NO_FILE' ? "-L ${exome_bed} --interval-padding ${params.exome_padding}" : ''
    """
    gatk ${avail_mem} HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        -O ${prefix}.vcf.gz \\
        --native-pair-hmm-threads $task.cpus \\
        ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
    END_VERSIONS
    """
}

// COMMENTED OUT - T2T not needed for exome analysis
// process GATK_HAPLOTYPECALLER_T2T_BWAMEM {
//     tag "$meta.id - $meta.qc_tool - t2t - bwamem"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/gatk/bwamem/${meta.qc_tool}/t2t", mode: 'copy'
//
//     container 'broadinstitute/gatk:4.5.0.0'
//
//     input:
//     tuple val(meta), path(bam), path(bai)
//
//     output:
//     tuple val(meta), path("*.vcf.gz"),     emit: vcf
//     tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
//     path "versions.yml",                   emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_t2t_bwamem"
//     def reference = "/reference/t2t/hs1.fa"
//     def avail_mem = task.memory ? "--java-options \"-Xmx${task.memory.toGiga()}g\"" : ''
//     """
//     gatk ${avail_mem} HaplotypeCaller \\
//         -R ${reference} \\
//         -I ${bam} \\
//         -O ${prefix}.vcf.gz \\
//         --native-pair-hmm-threads $task.cpus \\
//         $args
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+')
//     END_VERSIONS
//     """
// }
