// COMMENTED OUT - BOWTIE2 not used for exome analysis
// process MANTA_HG37_BOWTIE2 {
//     tag "$meta.id - $meta.qc_tool - hg37 - bowtie2"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/manta/bowtie2/${meta.qc_tool}/hg37", mode: 'copy'
//
//     container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
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
//     """
//     # Configure Manta workflow
//     configManta.py \\
//         --bam ${bam} \\
//         --referenceFasta ${reference} \\
//         --runDir manta_run \\
//         $args
//
//     # Run Manta workflow
//     chmod +x manta_run/runWorkflow.py
//     python2 manta_run/runWorkflow.py \\
//         -m local \\
//         -j $task.cpus
//
//     # Rename and move output files
//     mv manta_run/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
//     mv manta_run/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
//     END_VERSIONS
//     """
// }

// COMMENTED OUT - BOWTIE2 not used for exome analysis
// process MANTA_HG38_BOWTIE2 {
//     tag "$meta.id - $meta.qc_tool - hg38 - bowtie2"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/manta/bowtie2/${meta.qc_tool}/hg38", mode: 'copy'
//
//     container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
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
//     """
//     # Configure Manta workflow
//     configManta.py \\
//         --bam ${bam} \\
//         --referenceFasta ${reference} \\
//         --runDir manta_run \\
//         $args
//
//     # Run Manta workflow
//     chmod +x manta_run/runWorkflow.py
//     python2 manta_run/runWorkflow.py \\
//         -m local \\
//         -j $task.cpus
//
//     # Rename and move output files
//     mv manta_run/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
//     mv manta_run/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
//     END_VERSIONS
//     """
// }

// COMMENTED OUT - BOWTIE2 not used for exome analysis
// process MANTA_T2T_BOWTIE2 {
//     tag "$meta.id - $meta.qc_tool - t2t - bowtie2"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/manta/bowtie2/${meta.qc_tool}/t2t", mode: 'copy'
//
//     container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
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
//     """
//     # Configure Manta workflow
//     configManta.py \\
//         --bam ${bam} \\
//         --referenceFasta ${reference} \\
//         --runDir manta_run \\
//         $args
//
//     # Run Manta workflow
//     chmod +x manta_run/runWorkflow.py
//     python2 manta_run/runWorkflow.py \\
//         -m local \\
//         -j $task.cpus
//
//     # Rename and move output files
//     mv manta_run/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
//     mv manta_run/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
//     END_VERSIONS
//     """
// }

// COMMENTED OUT - HG37 not needed for exome analysis
// process MANTA_HG37_BWAMEM {
//     tag "$meta.id - $meta.qc_tool - hg37 - bwamem"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/manta/bwamem/${meta.qc_tool}/hg37", mode: 'copy'
//
//     container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
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
//     """
//     # Configure Manta workflow
//     configManta.py \\
//         --bam ${bam} \\
//         --referenceFasta ${reference} \\
//         --runDir manta_run \\
//         $args
//
//     # Run Manta workflow
//     chmod +x manta_run/runWorkflow.py
//     python2 manta_run/runWorkflow.py \\
//         -m local \\
//         -j $task.cpus
//
//     # Rename and move output files
//     mv manta_run/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
//     mv manta_run/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
//     END_VERSIONS
//     """
// }

process MANTA_HG38_BWAMEM {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/05_variant_calling/${meta.run}/manta/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'

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
    def exome_flag = params.exome_mode ? "--exome" : ''
    def call_regions = params.exome_mode && exome_bed.name != 'NO_FILE' ? "--callRegions ${exome_bed}" : ''
    """
    # Configure Manta workflow
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${reference} \\
        --runDir manta_run \\
        ${exome_flag} \\
        ${call_regions} \\
        $args

    # Run Manta workflow
    chmod +x manta_run/runWorkflow.py
    python2 manta_run/runWorkflow.py \\
        -m local \\
        -j $task.cpus

    # Rename and move output files
    mv manta_run/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
    mv manta_run/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
    END_VERSIONS
    """
}

// COMMENTED OUT - T2T not needed for exome analysis
// process MANTA_T2T_BWAMEM {
//     tag "$meta.id - $meta.qc_tool - t2t - bwamem"
//     label 'process_high'
//
//     publishDir "${params.outdir}/04_variant_calling/manta/bwamem/${meta.qc_tool}/t2t", mode: 'copy'
//
//     container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
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
//     """
//     # Configure Manta workflow
//     configManta.py \\
//         --bam ${bam} \\
//         --referenceFasta ${reference} \\
//         --runDir manta_run \\
//         $args
//
//     # Run Manta workflow
//     chmod +x manta_run/runWorkflow.py
//     python2 manta_run/runWorkflow.py \\
//         -m local \\
//         -j $task.cpus
//
//     # Rename and move output files
//     mv manta_run/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
//     mv manta_run/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi
//
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
//     END_VERSIONS
//     """
// }
