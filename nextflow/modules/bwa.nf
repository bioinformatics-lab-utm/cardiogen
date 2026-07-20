// COMMENTED OUT - HG37 not needed for exome analysis
// process BWAMEM_HG37 {
//     tag "$meta.id - $meta.qc_tool - hg37"
//     label 'process_high'
//
//     publishDir "${params.outdir}/03_bwamem/${meta.qc_tool}/hg37", mode: 'copy'
//
//     container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'
//
//     input:
//     tuple val(meta), path(reads)
//
//     output:
//     tuple val(meta), path("*.bam"),       emit: bam
//     tuple val(meta), path("*.bam.bai"),   emit: bai
//     path "versions.yml",                  emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_hg37"
//     def reference = "/reference/hg37/hg19.fa.gz"
//     def sort_mem_mb = Math.max(1, (task.memory.toMega() * 0.5 / task.cpus).intValue())
//
//     if (meta.single_end) {
//         """
//         # Check if BWA index exists, if not create it
//         if [ ! -f "${reference}.bwt" ]; then
//             echo "BWA index not found. Creating index for hg37..."
//             bwa index ${reference}
//             echo "BWA index created successfully."
//         fi
//
//         bwa mem \\
//             -t $task.cpus \\
//             -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
//             $args \\
//             ${reference} \\
//             ${reads[0]} \\
//             | samtools view -@ $task.cpus -bS - \\
//             | samtools sort -@ $task.cpus -m ${sort_mem_mb}M -o ${prefix}.bam -
//
//         samtools index -@ $task.cpus ${prefix}.bam
//
//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
//             samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
//         END_VERSIONS
//         """
//     } else {
//         """
//         # Check if BWA index exists, if not create it
//         if [ ! -f "${reference}.bwt" ]; then
//             echo "BWA index not found. Creating index for hg37..."
//             bwa index ${reference}
//             echo "BWA index created successfully."
//         fi
//
//         bwa mem \\
//             -t $task.cpus \\
//             -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
//             $args \\
//             ${reference} \\
//             ${reads[0]} \\
//             ${reads[1]} \\
//             | samtools view -@ $task.cpus -bS - \\
//             | samtools sort -@ $task.cpus -m ${sort_mem_mb}M -o ${prefix}.bam -
//
//         samtools index -@ $task.cpus ${prefix}.bam
//
//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
//             samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
//         END_VERSIONS
//         """
//     }
// }

process BWAMEM_HG38 {
    tag "$meta.run - $meta.id - $meta.qc_tool - hg38"
    label 'process_high'

    publishDir "${params.outdir}/04_bwamem/${meta.run}/${meta.qc_tool}/hg38", mode: 'copy'

    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"),       emit: bam
    tuple val(meta), path("*.bam.bai"),   emit: bai
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38"
    def reference = "/reference/hg38/hg38.fa.gz"
    def sort_mem_mb = Math.max(1, (task.memory.toMega() * 0.5 / task.cpus).intValue())

    if (meta.single_end) {
        """
        # Check if BWA index exists, if not create it
        if [ ! -f "${reference}.bwt" ]; then
            echo "BWA index not found. Creating index for hg38..."
            bwa index ${reference}
            echo "BWA index created successfully."
        fi

        bwa mem \\
            -t $task.cpus \\
            -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
            $args \\
            ${reference} \\
            ${reads[0]} \\
            | samtools view -@ $task.cpus -bS - \\
            | samtools sort -@ $task.cpus -m ${sort_mem_mb}M -o ${prefix}.bam -

        samtools index -@ $task.cpus ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
            samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        # Check if BWA index exists, if not create it
        if [ ! -f "${reference}.bwt" ]; then
            echo "BWA index not found. Creating index for hg38..."
            bwa index ${reference}
            echo "BWA index created successfully."
        fi

        bwa mem \\
            -t $task.cpus \\
            -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
            $args \\
            ${reference} \\
            ${reads[0]} \\
            ${reads[1]} \\
            | samtools view -@ $task.cpus -bS - \\
            | samtools sort -@ $task.cpus -m ${sort_mem_mb}M -o ${prefix}.bam -

        samtools index -@ $task.cpus ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
            samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
        END_VERSIONS
        """
    }
}

// COMMENTED OUT - T2T not needed for exome analysis
// process BWAMEM_T2T {
//     tag "$meta.id - $meta.qc_tool - t2t"
//     label 'process_high'
//
//     publishDir "${params.outdir}/03_bwamem/${meta.qc_tool}/t2t", mode: 'copy'
//
//     container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'
//
//     input:
//     tuple val(meta), path(reads)
//
//     output:
//     tuple val(meta), path("*.bam"),       emit: bam
//     tuple val(meta), path("*.bam.bai"),   emit: bai
//     path "versions.yml",                  emit: versions
//
//     when:
//     task.ext.when == null || task.ext.when
//
//     script:
//     def args = task.ext.args ?: ''
//     def prefix = "${meta.id}_${meta.qc_tool}_t2t"
//     def reference = "/reference/t2t/hs1.fa.gz"
//     def sort_mem_mb = Math.max(1, (task.memory.toMega() * 0.5 / task.cpus).intValue())
//
//     if (meta.single_end) {
//         """
//         # Check if BWA index exists, if not create it
//         if [ ! -f "${reference}.bwt" ]; then
//             echo "BWA index not found. Creating index for T2T..."
//             bwa index ${reference}
//             echo "BWA index created successfully."
//         fi
//
//         bwa mem \\
//             -t $task.cpus \\
//             -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
//             $args \\
//             ${reference} \\
//             ${reads[0]} \\
//             | samtools view -@ $task.cpus -bS - \\
//             | samtools sort -@ $task.cpus -m ${sort_mem_mb}M -o ${prefix}.bam -
//
//         samtools index -@ $task.cpus ${prefix}.bam
//
//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
//             samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
//         END_VERSIONS
//         """
//     } else {
//         """
//         # Check if BWA index exists, if not create it
//         if [ ! -f "${reference}.bwt" ]; then
//             echo "BWA index not found. Creating index for T2T..."
//             bwa index ${reference}
//             echo "BWA index created successfully."
//         fi
//
//         bwa mem \\
//             -t $task.cpus \\
//             -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
//             $args \\
//             ${reference} \\
//             ${reads[0]} \\
//             ${reads[1]} \\
//             | samtools view -@ $task.cpus -bS - \\
//             | samtools sort -@ $task.cpus -m ${sort_mem_mb}M -o ${prefix}.bam -
//
//         samtools index -@ $task.cpus ${prefix}.bam
//
//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
//             samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
//         END_VERSIONS
//         """
//     }
// }
