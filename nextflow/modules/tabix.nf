#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         Tabix Indexing Module
========================================================================================
 Index VCF files using tabix
----------------------------------------------------------------------------------------
*/

process TABIX_INDEX {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf), path("${vcf}.tbi"), emit: indexed_vcf
    path "versions.yml", emit: versions

    script:
    """
    tabix -p vcf ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

/*
========================================================================================
                         BCF to VCF.GZ Conversion Module
========================================================================================
 Convert BCF files to VCF.GZ format for hap.py compatibility
----------------------------------------------------------------------------------------
*/

process BCF_TO_VCFGZ {
    tag "\${sample_id}_\${caller}_\${aligner}_\${qc}_\${reference}"
    label 'process_medium'

    container 'quay.io/biocontainers/bcftools:1.17--haef29d1_0'

    input:
    tuple val(sample_id), val(caller), val(aligner), val(qc), val(reference), path(bcf), path(bcf_idx)

    output:
    tuple val(sample_id), val(caller), val(aligner), val(qc), val(reference), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def output_name = "${sample_id}_${caller}_${aligner}_${qc}_${reference}.vcf.gz"
    """
    bcftools view -O z -o ${output_name} ${bcf}
    bcftools index -t ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools *//')
    END_VERSIONS
    """
}
