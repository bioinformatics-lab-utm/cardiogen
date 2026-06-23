#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         VCF Merging Module
========================================================================================
 Merge VCF files from multiple variant callers to create consensus callsets
 Uses bcftools merge with appropriate strategies:
 - SNP/INDEL: GATK + DeepVariant
 - SV: DELLY + MANTA
----------------------------------------------------------------------------------------
*/

process MERGE_SNPINDEL_CALLERS {
    tag "${sample_id}_${aligner}_${qc}_${reference}_gatk_deepvariant"
    label 'process_medium'
    publishDir "${params.outdir}/05_variant_calling/merged/snp_indel/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(aligner), val(qc), val(reference), path(gatk_vcf, stageAs: "gatk.vcf.gz"), path(gatk_idx, stageAs: "gatk.vcf.gz.tbi"), path(dv_vcf, stageAs: "deepvariant.vcf.gz"), path(dv_idx, stageAs: "deepvariant.vcf.gz.tbi")

    output:
    tuple val(sample_id), val(aligner), val(qc), val(reference), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: merged_vcf
    path "*.stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def prefix = "${sample_id}_${aligner}_${qc}_${reference}_gatk_deepvariant"
    """
    # Merge GATK and DeepVariant VCFs using bcftools
    # Using -m all to merge all records at the same position
    # --force-samples to allow merging when sample names differ
    
    echo "Merging GATK and DeepVariant VCFs..."
    bcftools merge \\
        --merge all \\
        --force-samples \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        gatk.vcf.gz deepvariant.vcf.gz
    
    # Index merged VCF
    bcftools index -t ${prefix}.merged.vcf.gz
    
    # Generate stats
    echo "=== GATK VCF Stats ===" > ${prefix}.stats.txt
    bcftools stats gatk.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== DeepVariant VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats deepvariant.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== Merged VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats ${prefix}.merged.vcf.gz | grep "^SN" >> ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}

process MERGE_SV_CALLERS {
    tag "${sample_id}_${aligner}_${qc}_${reference}_delly_manta"
    label 'process_medium'
    publishDir "${params.outdir}/05_variant_calling/merged/sv/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(aligner), val(qc), val(reference), path(delly_vcf, stageAs: "delly.bcf"), path(delly_idx, stageAs: "delly.bcf.csi"), path(manta_vcf, stageAs: "manta.vcf.gz"), path(manta_idx, stageAs: "manta.vcf.gz.tbi")

    output:
    tuple val(sample_id), val(aligner), val(qc), val(reference), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: merged_vcf
    path "*.stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def prefix = "${sample_id}_${aligner}_${qc}_${reference}_delly_manta"
    """
    # Convert Delly BCF to VCF.GZ
    echo "Converting Delly BCF to VCF.GZ..."
    bcftools view -O z -o delly_converted.vcf.gz delly.bcf
    bcftools index -t delly_converted.vcf.gz
    
    # For SV, we concatenate instead of merge because SV callers
    # often call different variants at similar positions
    # This keeps all unique variants from both callers
    
    echo "Concatenating Delly and Manta VCFs..."
    bcftools concat \\
        --allow-overlaps \\
        --remove-duplicates \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        delly_converted.vcf.gz manta.vcf.gz
    
    # Index merged VCF
    bcftools index -t ${prefix}.merged.vcf.gz
    
    # Generate stats
    echo "=== Delly VCF Stats ===" > ${prefix}.stats.txt
    bcftools stats delly_converted.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== Manta VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats manta.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== Merged VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats ${prefix}.merged.vcf.gz | grep "^SN" >> ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}

process MERGE_SNPINDEL_CALLERS_COMBINED {
    tag "${sample_id}_${aligner}_combined_${reference}_gatk_deepvariant"
    label 'process_medium'
    publishDir "${params.outdir}/05_variant_calling/merged/snp_indel_combined/${aligner}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(aligner), val(reference), path(gatk_vcf, stageAs: "gatk.vcf.gz"), path(gatk_idx, stageAs: "gatk.vcf.gz.tbi"), path(dv_vcf, stageAs: "deepvariant.vcf.gz"), path(dv_idx, stageAs: "deepvariant.vcf.gz.tbi")

    output:
    tuple val(sample_id), val(aligner), val(reference), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: merged_vcf
    path "*.stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def prefix = "${sample_id}_${aligner}_combined_${reference}_gatk_deepvariant"
    """
    # Merge GATK and DeepVariant VCFs using bcftools
    echo "Merging GATK and DeepVariant VCFs (combined QC)..."
    bcftools merge \\
        --merge all \\
        --force-samples \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        gatk.vcf.gz deepvariant.vcf.gz
    
    # Index merged VCF
    bcftools index -t ${prefix}.merged.vcf.gz
    
    # Generate stats
    echo "=== GATK VCF Stats ===" > ${prefix}.stats.txt
    bcftools stats gatk.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== DeepVariant VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats deepvariant.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== Merged VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats ${prefix}.merged.vcf.gz | grep "^SN" >> ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}

process MERGE_SV_CALLERS_COMBINED {
    tag "${sample_id}_${aligner}_combined_${reference}_delly_manta"
    label 'process_medium'
    publishDir "${params.outdir}/05_variant_calling/merged/sv_combined/${aligner}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(aligner), val(reference), path(delly_vcf, stageAs: "delly.bcf"), path(delly_idx, stageAs: "delly.bcf.csi"), path(manta_vcf, stageAs: "manta.vcf.gz"), path(manta_idx, stageAs: "manta.vcf.gz.tbi")

    output:
    tuple val(sample_id), val(aligner), val(reference), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: merged_vcf
    path "*.stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def prefix = "${sample_id}_${aligner}_combined_${reference}_delly_manta"
    """
    # Convert Delly BCF to VCF.GZ
    echo "Converting Delly BCF to VCF.GZ..."
    bcftools view -O z -o delly_converted.vcf.gz delly.bcf
    bcftools index -t delly_converted.vcf.gz
    
    # For SV, we concatenate instead of merge because SV callers
    # often call different variants at similar positions
    # This keeps all unique variants from both callers
    
    echo "Concatenating Delly and Manta VCFs (combined QC)..."
    bcftools concat \\
        --allow-overlaps \\
        --remove-duplicates \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        delly_converted.vcf.gz manta.vcf.gz
    
    # Index merged VCF
    bcftools index -t ${prefix}.merged.vcf.gz
    
    # Generate stats
    echo "=== Delly VCF Stats ===" > ${prefix}.stats.txt
    bcftools stats delly_converted.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== Manta VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats manta.vcf.gz | grep "^SN" >> ${prefix}.stats.txt
    
    echo "" >> ${prefix}.stats.txt
    echo "=== Merged VCF Stats ===" >> ${prefix}.stats.txt
    bcftools stats ${prefix}.merged.vcf.gz | grep "^SN" >> ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
