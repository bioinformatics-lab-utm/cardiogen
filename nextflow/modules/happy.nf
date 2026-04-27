#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         hap.py Benchmarking Module
========================================================================================
 Compare variant calling results against truth VCF using hap.py
 Used for: DeepVariant, GATK HaplotypeCaller, Delly (for indels/small variants)
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
    # For hg37: normalize chromosome names (remove "chr" prefix to match GIAB truth VCF format)
    # For hg38/t2t: use fixchr flag to handle chromosome naming automatically
    if [[ "${reference}" == "hg37" ]]; then
        echo "HG37 detected: Normalizing chromosome names to match GIAB truth format..."
        
        # Install bcftools and samtools for chromosome normalization
        echo "Installing bcftools and samtools..."
        apt-get update -qq > /dev/null 2>&1
        apt-get install -y -qq bcftools samtools > /dev/null 2>&1
        
        # Create chromosome mapping file (chr1->1, chr2->2, etc.)
        cat > chr_rename.txt <<'CHRMAP'
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr10	10
chr11	11
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr19	19
chr20	20
chr21	21
chr22	22
chrX	X
chrY	Y
chrM	MT
CHRMAP

        # Check if query VCF uses "chr" prefix
        FIRST_CHR=\$(bcftools view -H ${query_vcf} 2>/dev/null | head -n1 | awk '{print \$1}' || true)
        
        if [[ "\${FIRST_CHR}" == chr* ]]; then
            echo "Query VCF uses 'chr' prefix, normalizing to match truth VCF (numeric only)..."
            bcftools annotate --rename-chrs chr_rename.txt -O z -o query_normalized.vcf.gz ${query_vcf}
            bcftools index -t query_normalized.vcf.gz
            QUERY_VCF="query_normalized.vcf.gz"
            echo "Chromosome names normalized"
        else
            echo "Query VCF already uses numeric chromosome names"
            QUERY_VCF="${query_vcf}"
        fi
        
        # Create reference version with numeric chromosome names
        echo "Creating reference FASTA with numeric chromosome names..."
        cat ${reference_fasta} | sed 's/>chr/>/' > reference_nochr.fa
        samtools faidx reference_nochr.fa
        REF_FASTA="reference_nochr.fa"
        FIXCHR_FLAG=""
    else
        echo "Using --fixchr flag for chromosome name handling"
        QUERY_VCF="${query_vcf}"
        REF_FASTA="${reference_fasta}"
        FIXCHR_FLAG="--fixchr"
    fi

    # Run hap.py comparison
    /opt/hap.py/bin/hap.py \\
        ${truth_vcf} \\
        \${QUERY_VCF} \\
        -f ${confident_bed} \\
        -r \${REF_FASTA} \\
        -o ${sample_id}_${caller}_${aligner}_${qc}_${reference} \\
        --threads ${task.cpus} \\
        --engine=vcfeval \\
        --pass-only \\
        \${FIXCHR_FLAG}

    # Generate version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: \$(/opt/hap.py/bin/hap.py --version 2>&1 | grep "hap.py" | sed 's/.*hap.py *//')
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools *//')
    END_VERSIONS
    """
}
