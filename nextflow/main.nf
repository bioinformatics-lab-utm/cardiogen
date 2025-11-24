#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
                         Cardiogen ONT Data Processing Workflow
========================================================================================
 Nextflow workflow for processing Oxford Nanopore sequencing data
 Started with FastQC and MultiQC quality control steps
----------------------------------------------------------------------------------------
*/

// Import modules
include { FASTP } from './modules/fastp'
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'
include { CUTADAPT } from './modules/cutadapt'
include { TRIMMOMATIC } from './modules/trimmomatic'
include { CUTADAPT_CLEAN } from './modules/cutadapt_clean'
include { TRIMMOMATIC_FINAL } from './modules/trimmomatic_final'
include { FASTQC_FASTP } from './modules/fastqc_fastp'
include { FASTQC_CUTADAPT } from './modules/fastqc_cutadapt'
include { FASTQC_TRIMMOMATIC } from './modules/fastqc_trimmomatic'
include { FASTQC_COMBINED } from './modules/fastqc_combined'
include { MULTIQC_FASTP } from './modules/multiqc_fastp'
include { MULTIQC_CUTADAPT } from './modules/multiqc_cutadapt'
include { MULTIQC_TRIMMOMATIC } from './modules/multiqc_trimmomatic'
include { MULTIQC_COMBINED } from './modules/multiqc_combined'

// Parameters
params.input_dir = "${projectDir}/test_data/ont_data"
params.outdir = "${projectDir}/results"
params.pattern = "*_{1,2}.fastq.gz"

// Help message
def helpMessage() {
    log.info"""
    ================================================================
                Cardiogen ONT Data Processing Workflow
    ================================================================
    
    Usage:
    nextflow run main.nf [options]
    
    Options:
    --input_dir     Path to directory containing FASTQ files (default: ${params.input_dir})
    --outdir        Output directory for results (default: ${params.outdir})
    --pattern       File pattern to match FASTQ files (default: ${params.pattern})
    --help          Show this help message
    
    Example:
    nextflow run main.nf --input_dir /path/to/fastq --outdir /path/to/results
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Main workflow
workflow {
    // Print workflow information
    log.info """
    ================================================================
                Cardiogen ONT Data Processing Workflow
    ================================================================
    Input directory : ${params.input_dir}
    Output directory: ${params.outdir}
    File pattern    : ${params.pattern}
    ================================================================
    """.stripIndent()

    // Create input channel from FASTQ files
    fastq_ch = Channel
        .fromFilePairs("${params.input_dir}/${params.pattern}", checkIfExists: true)
        .map { sample_id, files ->
            def meta = [:]
            meta.id = sample_id
            meta.single_end = false  // Paired-end data
            [meta, files]
        }

    // Run FastQC on raw reads first
    FASTQC(fastq_ch)

    // Collect FastQC reports for MultiQC (extract only the zip files)
    fastqc_reports = FASTQC.out.zip.map { meta, files -> files }.collect()

    // Run MultiQC on FastQC results
    MULTIQC(fastqc_reports)

    // Run all preprocessing tools independently on raw data
    FASTP(fastq_ch)
    CUTADAPT(fastq_ch)
    TRIMMOMATIC(fastq_ch)
    
    // Combined Cutadapt + Trimmomatic process (two-step chain)
    CUTADAPT_CLEAN(fastq_ch)
    TRIMMOMATIC_FINAL(CUTADAPT_CLEAN.out.reads)

    // Run FastQC on processed results
    FASTQC_FASTP(FASTP.out.reads)
    FASTQC_CUTADAPT(CUTADAPT.out.reads)
    FASTQC_TRIMMOMATIC(TRIMMOMATIC.out.paired_reads)
    FASTQC_COMBINED(TRIMMOMATIC_FINAL.out.paired_reads)

    // Run separate MultiQC for each processing step
    MULTIQC_FASTP(FASTQC_FASTP.out.zip.map { meta, files -> files }.collect())
    MULTIQC_CUTADAPT(FASTQC_CUTADAPT.out.zip.map { meta, files -> files }.collect())
    MULTIQC_TRIMMOMATIC(FASTQC_TRIMMOMATIC.out.zip.map { meta, files -> files }.collect())
    MULTIQC_COMBINED(FASTQC_COMBINED.out.zip.map { meta, files -> files }.collect())

    // Print completion message
    workflow.onComplete {
        log.info """
        ================================================================
                        Pipeline completed successfully!
        ================================================================
        Results are saved in: results/
        
        FastQC reports (raw data): results/01_QC/fastqc/
        MultiQC report (raw data QC): results/01_QC/multiqc/multiqc_report.html
        
        Fastp results: results/02_fastp/
        Fastp FastQC reports: results/02_fastp/fastqc/
        Fastp MultiQC report: results/02_fastp/multiqc/multiqc_report.html
        
        Cutadapt results: results/03_cutadapt/
        Cutadapt FastQC reports: results/03_cutadapt/fastqc/
        Cutadapt MultiQC report: results/03_cutadapt/multiqc/multiqc_report.html
        
        Trimmomatic results: results/04_trimmomatic/
        Trimmomatic FastQC reports: results/04_trimmomatic/fastqc/
        Trimmomatic MultiQC report: results/04_trimmomatic/multiqc/multiqc_report.html
        
        Combined Cutadapt+Trimmomatic results: results/05_cutadapt_trimmomatic/
        Combined FastQC reports: results/05_cutadapt_trimmomatic/fastqc/
        Combined MultiQC report: results/05_cutadapt_trimmomatic/multiqc/multiqc_report.html
        ================================================================
        """.stripIndent()
    }
}

workflow.onError {
    log.error "Pipeline execution failed: ${workflow.errorMessage}"
}
