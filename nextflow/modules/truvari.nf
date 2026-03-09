#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         Truvari Benchmarking Module
========================================================================================
 Compare structural variant calling results against truth VCF using Truvari
 Used for: Manta, Delly
----------------------------------------------------------------------------------------
*/

process TRUVARI_COMPARE {
    tag "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    label 'process_medium'
    publishDir "${params.outdir}/05_truvari_comparison/${caller}/${aligner}/${qc}/${reference}", mode: 'copy'

    container 'python:3.11-slim'

    input:
    tuple val(sample_id), val(caller), val(aligner), val(qc), val(reference), path(query_vcf), path(query_vcf_idx)
    path truth_vcf
    path truth_vcf_idx
    path confident_bed
    path reference_fasta
    path reference_fai

    output:
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}/*", emit: results
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}/summary.json", emit: summary
    path "${sample_id}_${caller}_${aligner}_${qc}_${reference}/summary.csv", emit: summary_csv
    path "versions.yml", emit: versions

    script:
    def output_dir = "${sample_id}_${caller}_${aligner}_${qc}_${reference}"
    """
    # Install dependencies
    apt-get update -qq
    apt-get install -y -qq bcftools build-essential gcc make zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev > /dev/null 2>&1
    pip install truvari==4.2.2 --quiet

    # Convert BCF to VCF.GZ if needed (for Delly output)
    if [[ ${query_vcf} == *.bcf ]]; then
        echo "Converting BCF to VCF.GZ..."
        bcftools view -O z -o query.vcf.gz ${query_vcf}
        bcftools index -t query.vcf.gz
        QUERY_VCF="query.vcf.gz"
    else
        QUERY_VCF="${query_vcf}"
    fi

    # Run Truvari benchmarking
    truvari bench \\
        -b ${truth_vcf} \\
        -c \${QUERY_VCF} \\
        -f ${reference_fasta} \\
        -o ${output_dir} \\
        --passonly \\
        -p 0.00 \\
        -P 0.70 \\
        -O 0.00 \\
        --pick multi \\
        --includebed ${confident_bed}

    # Convert summary.json to CSV for easier comparison with hap.py
    python3 << 'EOF'
import json
import csv

with open('${output_dir}/summary.json', 'r') as f:
    data = json.load(f)

# Write CSV with metrics
with open('${output_dir}/summary.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Metric', 'Value'])
    for key, value in data.items():
        if not isinstance(value, dict):
            writer.writerow([key, value])
EOF

    # Generate version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(truvari version 2>&1 | head -n1)
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools *//')
    END_VERSIONS
    """
}
