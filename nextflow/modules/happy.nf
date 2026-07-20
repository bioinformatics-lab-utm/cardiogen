#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         hap.py Comparison Module
========================================================================================
 Compare the pipeline SNP/indel calls (DeepVariant, GATK - hg38) against the laboratory
 reference VCFs on a per-sample basis. The reference VCF (originally GRCh37/hg19) is
 lifted to hg38 by CROSSMAP_LIFTOVER beforehand and passed in as the "truth" call set.

 The native hg38 TruSight Cardio target panel (reference/bed_hg38) is used as the confident
 region, and the hg38 FASTA + panel BED are read from the /reference mount (docker profile).
 Both call sets are hg38 and chr-prefixed, so no chromosome renaming/fixchr is needed.

 Comparison uses hap.py's default xcmp engine (not vcfeval): vcfeval aborts when the
 baseline (lifted reference) carries variants on a contig the query never references
 (e.g. chr21/chrY with no pipeline call, or liftover artifacts on *_alt contigs). xcmp
 instead counts those as FN/FP, and the panel BED keeps the analysis on the target.
----------------------------------------------------------------------------------------
*/

process HAPPY_COMPARE {
    tag "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    label 'process_medium'
    publishDir "${params.outdir}/06_comparison/happy/${caller}/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'jmcdani20/hap.py:v0.3.12'

    input:
    tuple val(sample_id), val(caller), val(aligner), val(qc), val(reference), path(query_vcf), path(query_vcf_idx), path(truth_vcf), path(truth_vcf_idx)

    output:
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}.*",             emit: results
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}.summary.csv",   emit: summary
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}.extended.csv",  emit: extended
    path "versions.yml",                                                       emit: versions

    script:
    def prefix  = "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    def hg38_fa = "/reference/hg38/hg38.fa"
    def bed     = "/reference/bed_hg38/TruSight_Cardio_TargetedRegions_v1.0.hg38.bed"
    """
    # Truth = lifted laboratory reference VCF, Query = pipeline call set (both hg38)
    /opt/hap.py/bin/hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        -f ${bed} \\
        -r ${hg38_fa} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        --pass-only

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: \$(/opt/hap.py/bin/hap.py --version 2>&1 | grep "hap.py" | sed 's/.*hap.py *//')
    END_VERSIONS
    """
}
