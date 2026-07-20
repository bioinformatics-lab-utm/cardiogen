#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         bcftools isec Comparison Module
========================================================================================
 Intersect the pipeline SNP/indel calls (DeepVariant, GATK - hg38) against the
 laboratory reference VCFs (lifted to hg38 by CROSSMAP_LIFTOVER) on a per-sample basis.

 Both call sets are normalized (split multiallelics + left-align) for a fair
 position+allele intersection. Output (in the per-comparison folder):
   0000.vcf -> records private to the pipeline call set   (pipeline-only / "FP-like")
   0001.vcf -> records private to the reference call set   (reference-only / "FN-like")
   0002.vcf -> shared records, as represented in pipeline  (concordant)
   0003.vcf -> shared records, as represented in reference
   summary.txt -> counts of the above
----------------------------------------------------------------------------------------
*/

process BCFTOOLS_ISEC {
    tag "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    label 'process_low'
    publishDir "${params.outdir}/06_comparison/bcftools_isec/${caller}/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(caller), val(aligner), val(qc), val(reference), path(query_vcf), path(query_vcf_idx), path(truth_vcf), path(truth_vcf_idx)

    output:
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}",             emit: results
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}/summary.txt", emit: summary
    path "versions.yml",                                                     emit: versions

    script:
    def outdir  = "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    def hg38_fa = "/reference/hg38/hg38.fa"
    def bed     = "/reference/bed_hg38/TruSight_Cardio_TargetedRegions_v1.0.hg38.bed"
    """
    # Normalize both call sets (split multiallelics + left-align) for a fair comparison.
    # -c w: downgrade REF/duplicate-allele mismatches from a fatal error to a warning
    #       (the lifted reference and panel VCFs can contain duplicate alleles after the
    #        -m-any split, which otherwise aborts bcftools with exit 255).
    # second pass -d exact: drop the exact-duplicate records that the split can create.
    bcftools norm -m-any -c w -f ${hg38_fa} ${query_vcf} -O u \\
        | bcftools norm -d exact -O z -o query.norm.vcf.gz
    bcftools index -t query.norm.vcf.gz
    bcftools norm -m-any -c w -f ${hg38_fa} ${truth_vcf} -O u \\
        | bcftools norm -d exact -O z -o reference.norm.vcf.gz
    bcftools index -t reference.norm.vcf.gz

    mkdir -p ${outdir}
    # Index 0 = pipeline (query), index 1 = reference (truth).
    # -R: restrict to the native hg38 TruSight Cardio panel so isec compares the same
    #     region as hap.py (the pipeline only called within the panel, so reference
    #     variants outside it are a coverage difference, not a calling disagreement).
    bcftools isec -R ${bed} -p ${outdir} query.norm.vcf.gz reference.norm.vcf.gz

    # Human-readable summary (|| true: grep -c exits 1 when there are zero records)
    {
        echo "Comparison      : ${sample_id} (${caller} / ${aligner} / ${qc} / ${reference})"
        echo "Pipeline (query): ${query_vcf}"
        echo "Reference (ref) : ${truth_vcf}"
        echo ""
        echo "Pipeline-only (FP-like) : \$(grep -vc '^#' ${outdir}/0000.vcf || true)"
        echo "Reference-only (FN-like): \$(grep -vc '^#' ${outdir}/0001.vcf || true)"
        echo "Shared (concordant)     : \$(grep -vc '^#' ${outdir}/0002.vcf || true)"
    } > ${outdir}/summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools *//')
    END_VERSIONS
    """
}
