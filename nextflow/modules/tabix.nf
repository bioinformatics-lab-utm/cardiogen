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
