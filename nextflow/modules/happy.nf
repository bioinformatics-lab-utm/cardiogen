#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         hap.py Benchmarking Module
========================================================================================
 Compare variant calling results against truth VCF using hap.py
 Used for: DeepVariant, GATK HaplotypeCaller, Octopus
----------------------------------------------------------------------------------------
*/

process HAPPY_COMPARE {
    tag "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    label 'process_medium'
    publishDir "${params.outdir}/05_happy_comparison/${caller}/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'jmcdani20/hap.py:v0.3.12'

    input:
    tuple val(sample_id), val(caller), val(aligner), val(qc), val(reference), path(query_vcf), path(query_vcf_idx)
    path truth_vcf
    path truth_vcf_idx
    path confident_bed
    path reference_fasta
    path reference_fai

    output:
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}.*", emit: results
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}.summary.csv", emit: summary
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}.extended.csv", emit: extended
    path "versions.yml", emit: versions

    script:
    """
    # Run hap.py comparison
    /opt/hap.py/bin/hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        -f ${confident_bed} \\
        -r ${reference_fasta} \\
        -o ${sample_id}_${caller}_${aligner}_${qc}_${reference} \\
        --threads ${task.cpus} \\
        --engine=vcfeval \\
        --pass-only

    # Generate version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: \$(/opt/hap.py/bin/hap.py --version 2>&1 | grep "hap.py" | sed 's/.*hap.py *//')
    END_VERSIONS
    """
}
