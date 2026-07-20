process COUNT_READS {
    tag "$meta.run - $meta.id - $meta.qc_tool"
    label 'process_low'

    // Reuse the bwa/samtools mulled container so samtools is available.
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(bam), path(bai), env(READ_COUNT), emit: counted

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Count primary, properly-paired, mapped reads (-f 0x2 = properly paired,
    # -F 0x900 = exclude secondary + supplementary). Variant callers - Delly in
    # particular - need enough properly-paired reads to estimate insert-size library
    # parameters; empty/failed samples (e.g. negative controls) have ~0 such reads and
    # would otherwise abort and kill the whole run. The caller uses this count to skip
    # those samples gracefully.
    READ_COUNT=\$(samtools view -c -f 0x2 -F 0x900 ${bam})
    """
}
