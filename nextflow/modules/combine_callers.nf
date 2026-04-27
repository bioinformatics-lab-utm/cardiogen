#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                    Combine Variant Callers Module
========================================================================================
 Combine results from multiple variant callers to create consensus callsets
 
 Strategy:
 - SNP/INDEL: Combine GATK + DeepVariant using bcftools isec (intersection)
 - SV: Combine DELLY + MANTA using SURVIVOR merge
========================================================================================
*/

process COMBINE_SNP_INDEL {
    tag "${sample_id}_${aligner}_${qc}_${reference}"
    label 'process_medium'
    publishDir "${params.outdir}/04_variant_calling/combined/snp_indel/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(aligner), val(qc), val(reference), 
          path(gatk_vcf, stageAs: "gatk.vcf.gz"), path(gatk_idx, stageAs: "gatk.vcf.gz.tbi"), 
          path(dv_vcf, stageAs: "deepvariant.vcf.gz"), path(dv_idx, stageAs: "deepvariant.vcf.gz.tbi")

    output:
    tuple val(sample_id), val(aligner), val(qc), val(reference), 
          path("*_combined.vcf.gz"), path("*_combined.vcf.gz.tbi"), emit: combined_vcf
    path "*_stats.txt", emit: stats
    path "*_sites_union.txt", emit: union_sites, optional: true
    path "*_sites_intersection.txt", emit: intersection_sites, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = "${sample_id}_${aligner}_${qc}_${reference}"
    """
    # Create working directory for bcftools isec
    mkdir -p isec_output
    
    # Find intersection and union of variants from GATK and DeepVariant
    # -p: output directory
    # -n +2: output sites present in at least 2 files (both callers agree)
    # -w: list of sites
    
    echo "Finding common variants between GATK and DeepVariant..."
    bcftools isec \\
        -p isec_output \\
        -n +2 \\
        -w 1,2 \\
        gatk.vcf.gz \\
        deepvariant.vcf.gz
    
    # Merge the intersection variants (high confidence)
    # These are variants called by both GATK and DeepVariant
    echo "Creating high-confidence combined callset..."
    
    if [ -f "isec_output/0002.vcf" ] && [ -s "isec_output/0002.vcf" ]; then
        # Compress and index intersection
        bgzip -c isec_output/0002.vcf > temp_combined.vcf.gz
        bcftools index -t temp_combined.vcf.gz
        
        # Save intersection sites
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' temp_combined.vcf.gz > ${prefix}_sites_intersection.txt
    else
        # If no intersection, create union instead
        # Use --force-samples because GATK and DeepVariant have the same sample name
        echo "No common variants found, creating union instead..."
        bcftools merge -m all --force-samples gatk.vcf.gz deepvariant.vcf.gz -O z -o temp_combined.vcf.gz
        bcftools index -t temp_combined.vcf.gz
    fi
    
    # Clean the combined VCF to make it compatible with hap.py
    # Keep only the first sample (both samples are the same individual anyway)
    # This removes duplicate FORMAT fields that cause hap.py to crash
    echo "Cleaning combined VCF for hap.py compatibility..."
    FIRST_SAMPLE=\$(bcftools query -l temp_combined.vcf.gz | head -1)
    bcftools view -s \${FIRST_SAMPLE} -O z -o ${prefix}_combined.vcf.gz temp_combined.vcf.gz
    bcftools index -t ${prefix}_combined.vcf.gz
    
    # Remove temporary files to avoid conflicts with output pattern
    rm -f temp_combined.vcf.gz temp_combined.vcf.gz.tbi
    
    # Generate statistics
    echo "=== Combined VCF Statistics ===" > ${prefix}_stats.txt
    bcftools stats ${prefix}_combined.vcf.gz | grep "^SN" >> ${prefix}_stats.txt
    
    echo "" >> ${prefix}_stats.txt
    echo "=== GATK Input Stats ===" >> ${prefix}_stats.txt
    bcftools stats gatk.vcf.gz | grep "^SN" >> ${prefix}_stats.txt
    
    echo "" >> ${prefix}_stats.txt
    echo "=== DeepVariant Input Stats ===" >> ${prefix}_stats.txt
    bcftools stats deepvariant.vcf.gz | grep "^SN" >> ${prefix}_stats.txt
    
    # Count variants
    COMBINED_COUNT=\$(bcftools view -H ${prefix}_combined.vcf.gz | wc -l)
    GATK_COUNT=\$(bcftools view -H gatk.vcf.gz | wc -l)
    DV_COUNT=\$(bcftools view -H deepvariant.vcf.gz | wc -l)
    
    echo "" >> ${prefix}_stats.txt
    echo "=== Variant Counts ===" >> ${prefix}_stats.txt
    echo "GATK variants: \${GATK_COUNT}" >> ${prefix}_stats.txt
    echo "DeepVariant variants: \${DV_COUNT}" >> ${prefix}_stats.txt
    echo "Combined variants: \${COMBINED_COUNT}" >> ${prefix}_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}

process COMBINE_SV {
    tag "${sample_id}_${aligner}_${qc}_${reference}"
    label 'process_medium'
    publishDir "${params.outdir}/04_variant_calling/combined/sv/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(sample_id), val(aligner), val(qc), val(reference),
          path(delly_vcf, stageAs: "delly.vcf.gz"), path(delly_idx, stageAs: "delly.vcf.gz.tbi"),
          path(manta_vcf, stageAs: "manta.vcf.gz"), path(manta_idx, stageAs: "manta.vcf.gz.tbi")

    output:
    tuple val(sample_id), val(aligner), val(qc), val(reference),
          path("*_combined.vcf.gz"), path("*_combined.vcf.gz.tbi"), emit: combined_vcf
    path "*_stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def prefix = "${sample_id}_${aligner}_${qc}_${reference}_sv"
    """
    # Delly may output BCF format, normalize to VCF.GZ
    # bcftools view handles both BCF and VCF.GZ transparently
    echo "Normalizing Delly output format..."
    bcftools view -O z -o delly_normalized.vcf.gz delly.vcf.gz
    bcftools index -t delly_normalized.vcf.gz
    
    # For SV, concatenate all variants from both callers
    # SV callers often detect different/complementary variants
    echo "Combining SV calls from DELLY and MANTA..."
    
    bcftools concat \\
        --allow-overlaps \\
        --remove-duplicates \\
        --output-type z \\
        --output ${prefix}_combined.vcf.gz \\
        delly_normalized.vcf.gz manta.vcf.gz
    
    # Index combined VCF
    bcftools index -t ${prefix}_combined.vcf.gz
    
    # Generate statistics
    echo "=== Combined SV Statistics ===" > ${prefix}_stats.txt
    bcftools stats ${prefix}_combined.vcf.gz | grep "^SN" >> ${prefix}_stats.txt
    
    echo "" >> ${prefix}_stats.txt
    echo "=== DELLY Input Stats ===" >> ${prefix}_stats.txt
    bcftools stats delly_normalized.vcf.gz | grep "^SN" >> ${prefix}_stats.txt
    
    echo "" >> ${prefix}_stats.txt
    echo "=== MANTA Input Stats ===" >> ${prefix}_stats.txt
    bcftools stats manta.vcf.gz | grep "^SN" >> ${prefix}_stats.txt
    
    # Count SVs
    COMBINED_COUNT=\$(bcftools view -H ${prefix}_combined.vcf.gz | wc -l)
    DELLY_COUNT=\$(bcftools view -H delly_normalized.vcf.gz | wc -l)
    MANTA_COUNT=\$(bcftools view -H manta.vcf.gz | wc -l)
    
    echo "" >> ${prefix}_stats.txt
    echo "=== SV Counts ===" >> ${prefix}_stats.txt
    echo "DELLY SVs: \${DELLY_COUNT}" >> ${prefix}_stats.txt
    echo "MANTA SVs: \${MANTA_COUNT}" >> ${prefix}_stats.txt
    echo "Combined SVs: \${COMBINED_COUNT}" >> ${prefix}_stats.txt
    
    # Remove temporary files to avoid conflicts with output pattern
    rm -f delly_normalized.vcf.gz delly_normalized.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
