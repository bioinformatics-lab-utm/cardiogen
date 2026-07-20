#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         CrossMap Liftover Module
========================================================================================
 Lift the laboratory reference VCFs (reference/vcf_hdd/<run>/<sample>.vcf) from
 GRCh37/hg19 to hg38 so they can be compared positionally against the hg38 pipeline
 calls (hap.py, bcftools isec).

 Uses the same chain that was already used to lift the TruSight Cardio BED panel
 (reference/bed_file/hg19ToHg38.over.chain.gz). Both the chain and the hg38 FASTA are
 read from the /reference mount provided by the docker profile.
----------------------------------------------------------------------------------------
*/

process CROSSMAP_LIFTOVER {
    tag "${run} - ${sample}"
    label 'process_low'
    publishDir "${params.outdir}/06_comparison/${run}/lifted_reference", mode: 'copy'

    container 'python:3.11-slim'

    input:
    tuple val(key), val(run), val(sample), path(ref_vcf_hg19)

    output:
    tuple val(key), path("${sample}.hg38.vcf.gz"), path("${sample}.hg38.vcf.gz.tbi"), emit: lifted
    path "versions.yml", emit: versions

    script:
    def chain   = "/reference/bed_file/hg19ToHg38.over.chain.gz"
    def hg38_fa = "/reference/hg38/hg38.fa"
    """
    # Install liftover + VCF tooling
    apt-get update -qq
    apt-get install -y -qq bcftools tabix > /dev/null 2>&1
    pip install --quiet CrossMap==0.7.0

    # Lift hg19 -> hg38. CrossMap reads the target FASTA to update REF alleles/contigs.
    # The lifted VCF can be out of coordinate order, so sort + bgzip + index afterwards.
    CrossMap vcf ${chain} ${ref_vcf_hg19} ${hg38_fa} ${sample}.lifted.vcf

    bcftools sort ${sample}.lifted.vcf -O z -o ${sample}.hg38.vcf.gz
    tabix -p vcf ${sample}.hg38.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CrossMap: \$(CrossMap --version 2>&1 | head -n1)
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools *//')
    END_VERSIONS
    """
}
