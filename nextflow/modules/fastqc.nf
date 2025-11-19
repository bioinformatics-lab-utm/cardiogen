process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/01_QC/fastqc", mode: 'copy', pattern: "*_fastqc.html"
    publishDir "${params.outdir}/01_QC/fastqc/plots", mode: 'copy', pattern: "*/Images/*.png"

    container 'staphb/fastqc:0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    path "*/Images/*.png", emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Ensure we handle both single and paired-end data
    def input_files = meta.single_end ? reads : reads.join(' ')
    
    """
    fastqc \\
        $args \\
        --threads $task.cpus \\
        $input_files

    # Extract plots from ZIP files
    for zip_file in *_fastqc.zip; do
        unzip -q "\$zip_file"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/^.*FastQC v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1_fastqc.html
    touch ${prefix}_1_fastqc.zip
    if [ ! "${meta.single_end}" = "true" ]; then
        touch ${prefix}_2_fastqc.html
        touch ${prefix}_2_fastqc.zip
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/^.*FastQC v//')
    END_VERSIONS
    """
}