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
        bcftools view -O z -o query_temp.vcf.gz ${query_vcf}
        bcftools index -t query_temp.vcf.gz
        TEMP_VCF="query_temp.vcf.gz"
    else
        TEMP_VCF="${query_vcf}"
    fi

    # Normalize chromosome names (remove "chr" prefix to match GIAB truth VCF format)
    echo "Normalizing chromosome names..."
    
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
    # Ignore SIGPIPE error from head closing the pipe early
    FIRST_CHR=\$(bcftools view -H \${TEMP_VCF} 2>/dev/null | head -n1 | awk '{print \$1}' || true)
    
    if [[ "\${FIRST_CHR}" == chr* ]]; then
        echo "Query VCF uses 'chr' prefix, normalizing to match truth VCF..."
        bcftools annotate --rename-chrs chr_rename.txt -O z -o query.vcf.gz \${TEMP_VCF}
        bcftools index -t query.vcf.gz
        echo "Chromosome names normalized"
    else
        echo "Query VCF already uses numeric chromosome names, no normalization needed"
        cp \${TEMP_VCF} query.vcf.gz
        cp \${TEMP_VCF}.tbi query.vcf.gz.tbi
    fi
    
    QUERY_VCF="query.vcf.gz"

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
