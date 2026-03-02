process BOWTIE2_ALIGN_HG37 {
    tag "$meta.id - $meta.qc_tool - hg37"
    label 'process_high'

    publishDir "${params.outdir}/04_bowtie2/${meta.qc_tool}/hg37", mode: 'copy'

    container 'staphb/bowtie2:2.5.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"),       emit: bam
    tuple val(meta), path("*.bam.bai"),   emit: bai
    tuple val(meta), path("*.log"),       emit: log
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg37"
    def index_path = workflow.containerEngine ? "/reference/hg37/hg19" : "${params.hg37_index}/hg19"
    def reference = workflow.containerEngine ? "/reference/hg37/hg19.fa" : "${params.hg37_index}/hg19.fa"
    // Use fewer threads for samtools to give more memory per thread
    def samtools_threads = Math.min(32, task.cpus as int)
    def sort_mem_mb = Math.max(768, (task.memory.toMega() * 0.55 / samtools_threads).intValue())
    """
    # Check if Bowtie2 index exists, if not create it
    if [ ! -f "${index_path}.1.bt2" ]; then
        echo "Bowtie2 index not found. Creating index for hg37..."
        bowtie2-build --threads $task.cpus ${reference} ${index_path}
        echo "Bowtie2 index created successfully."
    fi

    bowtie2 \\
        -x ${index_path} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --rg-id ${meta.id} \\
        --rg SM:${meta.id} \\
        --rg PL:ILLUMINA \\
        --rg LB:${meta.id} \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}_bowtie2.log \\
        | samtools view -@ ${samtools_threads} -bS - \\
        | samtools sort -@ ${samtools_threads} -m ${sort_mem_mb}M -o ${prefix}.bam -

    samtools index -@ ${samtools_threads} ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
    END_VERSIONS
    """
}

process BOWTIE2_ALIGN_HG38 {
    tag "$meta.id - $meta.qc_tool - hg38"
    label 'process_high'

    publishDir "${params.outdir}/04_bowtie2/${meta.qc_tool}/hg38", mode: 'copy'

    container 'staphb/bowtie2:2.5.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"),       emit: bam
    tuple val(meta), path("*.bam.bai"),   emit: bai
    tuple val(meta), path("*.log"),       emit: log
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38"
    def index_path = workflow.containerEngine ? "/reference/hg38/hg38" : "${params.hg38_index}/hg38"
    def reference = workflow.containerEngine ? "/reference/hg38/hg38.fa" : "${params.hg38_index}/hg38.fa"
    // Use fewer threads for samtools to give more memory per thread
    def samtools_threads = Math.min(32, task.cpus as int)
    def sort_mem_mb = Math.max(768, (task.memory.toMega() * 0.55 / samtools_threads).intValue())
    """
    # Check if Bowtie2 index exists, if not create it
    if [ ! -f "${index_path}.1.bt2" ]; then
        echo "Bowtie2 index not found. Creating index for hg38..."
        bowtie2-build --threads $task.cpus ${reference} ${index_path}
        echo "Bowtie2 index created successfully."
    fi

    bowtie2 \\
        -x ${index_path} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --rg-id ${meta.id} \\
        --rg SM:${meta.id} \\
        --rg PL:ILLUMINA \\
        --rg LB:${meta.id} \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}_bowtie2.log \\
        | samtools view -@ ${samtools_threads} -bS - \\
        | samtools sort -@ ${samtools_threads} -m ${sort_mem_mb}M -o ${prefix}.bam -

    samtools index -@ ${samtools_threads} ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
    END_VERSIONS
    """
}

process BOWTIE2_ALIGN_T2T {
    tag "$meta.id - $meta.qc_tool - t2t"
    label 'process_high'

    publishDir "${params.outdir}/04_bowtie2/${meta.qc_tool}/t2t", mode: 'copy'

    container 'staphb/bowtie2:2.5.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"),       emit: bam
    tuple val(meta), path("*.bam.bai"),   emit: bai
    tuple val(meta), path("*.log"),       emit: log
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_t2t"
    def index_path = workflow.containerEngine ? "/reference/t2t/hs1" : "${params.t2t_index}/hs1"
    def reference = workflow.containerEngine ? "/reference/t2t/hs1.fa" : "${params.t2t_index}/hs1.fa"
    // Use fewer threads for samtools to give more memory per thread
    def samtools_threads = Math.min(32, task.cpus as int)
    def sort_mem_mb = Math.max(768, (task.memory.toMega() * 0.55 / samtools_threads).intValue())
    """
    # Check if Bowtie2 index exists, if not create it
    if [ ! -f "${index_path}.1.bt2" ]; then
        echo "Bowtie2 index not found. Creating index for T2T..."
        bowtie2-build --threads $task.cpus ${reference} ${index_path}
        echo "Bowtie2 index created successfully."
    fi

    bowtie2 \\
        -x ${index_path} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --rg-id ${meta.id} \\
        --rg SM:${meta.id} \\
        --rg PL:ILLUMINA \\
        --rg LB:${meta.id} \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}_bowtie2.log \\
        | samtools view -@ ${samtools_threads} -bS - \\
        | samtools sort -@ ${samtools_threads} -m ${sort_mem_mb}M -o ${prefix}.bam -

    samtools index -@ ${samtools_threads} ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/ .*\$//')
    END_VERSIONS
    """
}
