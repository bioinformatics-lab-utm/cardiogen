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
include { CUTADAPT } from './modules/cutadapt'
include { TRIMMOMATIC } from './modules/trimmomatic'
include { BOWTIE2_ALIGN_HG37 } from './modules/bowtie2'
include { BOWTIE2_ALIGN_HG38 } from './modules/bowtie2'
include { BOWTIE2_ALIGN_T2T } from './modules/bowtie2'
include { BWAMEM_HG37 } from './modules/bwa'
include { BWAMEM_HG38 } from './modules/bwa'
include { BWAMEM_T2T } from './modules/bwa'
include { OCTOPUS_HG37_BOWTIE2 } from './modules/octopus'
include { OCTOPUS_HG38_BOWTIE2 } from './modules/octopus'
include { OCTOPUS_T2T_BOWTIE2 } from './modules/octopus'
include { OCTOPUS_HG37_BWAMEM } from './modules/octopus'
include { OCTOPUS_HG38_BWAMEM } from './modules/octopus'
include { OCTOPUS_T2T_BWAMEM } from './modules/octopus'
include { MANTA_HG37_BOWTIE2 } from './modules/manta'
include { MANTA_HG38_BOWTIE2 } from './modules/manta'
include { MANTA_T2T_BOWTIE2 } from './modules/manta'
include { MANTA_HG37_BWAMEM } from './modules/manta'
include { MANTA_HG38_BWAMEM } from './modules/manta'
include { MANTA_T2T_BWAMEM } from './modules/manta'

// Parameters
params.input_dir = "${projectDir}/test_data/ont_data"
params.outdir = "${projectDir}/results"
params.pattern = "*_{1,2}.fastq.gz"

// Reference genome parameters
params.hg37_index = "${projectDir}/reference/hg37"
params.hg38_index = "${projectDir}/reference/hg38"
params.t2t_index = "${projectDir}/reference/t2t"

// Help message
def helpMessage() {
    log.info"""
    ================================================================
                Cardiogen ONT Data Processing Workflow
    ================================================================
    
    Usage:
    nextflow run main.nf [options]
    
    Options:
    --input_dir       Path to directory containing FASTQ files (default: ${params.input_dir})
    --outdir          Output directory for results (default: ${params.outdir})
    --pattern         File pattern to match FASTQ files (default: ${params.pattern})
    --hg37_index      Path to hg37 bowtie2 index directory (default: ${params.hg37_index})
    --hg38_index      Path to hg38 bowtie2 index directory (default: ${params.hg38_index})
    --t2t_index       Path to T2T bowtie2 index directory (default: ${params.t2t_index})
    --help            Show this help message
    
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

    // Run preprocessing tools on raw data
    FASTP(fastq_ch)
    CUTADAPT(fastq_ch)
    TRIMMOMATIC(fastq_ch)

    // Prepare metadata for alignments with each genome
    
    // FASTP outputs
    fastp_hg37 = FASTP.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "fastp"
            [meta, tuple[1]]
        }
    
    fastp_hg38 = FASTP.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "fastp"
            [meta, tuple[1]]
        }
    
    fastp_t2t = FASTP.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "fastp"
            [meta, tuple[1]]
        }

    // CUTADAPT outputs
    cutadapt_hg37 = CUTADAPT.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "cutadapt"
            [meta, tuple[1]]
        }
    
    cutadapt_hg38 = CUTADAPT.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "cutadapt"
            [meta, tuple[1]]
        }
    
    cutadapt_t2t = CUTADAPT.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "cutadapt"
            [meta, tuple[1]]
        }

    // TRIMMOMATIC outputs
    trimmomatic_hg37 = TRIMMOMATIC.out.paired_reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "trimmomatic"
            [meta, tuple[1]]
        }
    
    trimmomatic_hg38 = TRIMMOMATIC.out.paired_reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "trimmomatic"
            [meta, tuple[1]]
        }
    
    trimmomatic_t2t = TRIMMOMATIC.out.paired_reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "trimmomatic"
            [meta, tuple[1]]
        }

    // Combine all inputs for each genome and run alignments
    hg37_inputs = fastp_hg37.mix(cutadapt_hg37, trimmomatic_hg37)
    hg38_inputs = fastp_hg38.mix(cutadapt_hg38, trimmomatic_hg38)
    t2t_inputs = fastp_t2t.mix(cutadapt_t2t, trimmomatic_t2t)

    // Run Bowtie2 alignments
    BOWTIE2_ALIGN_HG37(hg37_inputs)
    BOWTIE2_ALIGN_HG38(hg38_inputs)
    BOWTIE2_ALIGN_T2T(t2t_inputs)

    // Run BWA-MEM alignments
    BWAMEM_HG37(hg37_inputs)
    BWAMEM_HG38(hg38_inputs)
    BWAMEM_T2T(t2t_inputs)

    // ========================================
    // VARIANT CALLING WITH OCTOPUS
    // ========================================

    // Add aligner metadata to BAM files from Bowtie2
    bowtie2_hg37_bams = BOWTIE2_ALIGN_HG37.out.bam
        .join(BOWTIE2_ALIGN_HG37.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bowtie2"
            [new_meta, bam, bai]
        }

    bowtie2_hg38_bams = BOWTIE2_ALIGN_HG38.out.bam
        .join(BOWTIE2_ALIGN_HG38.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bowtie2"
            [new_meta, bam, bai]
        }

    bowtie2_t2t_bams = BOWTIE2_ALIGN_T2T.out.bam
        .join(BOWTIE2_ALIGN_T2T.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bowtie2"
            [new_meta, bam, bai]
        }

    // Add aligner metadata to BAM files from BWA-MEM
    bwamem_hg37_bams = BWAMEM_HG37.out.bam
        .join(BWAMEM_HG37.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bwamem"
            [new_meta, bam, bai]
        }

    bwamem_hg38_bams = BWAMEM_HG38.out.bam
        .join(BWAMEM_HG38.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bwamem"
            [new_meta, bam, bai]
        }

    bwamem_t2t_bams = BWAMEM_T2T.out.bam
        .join(BWAMEM_T2T.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bwamem"
            [new_meta, bam, bai]
        }

    // Run Octopus variant calling on Bowtie2 alignments
    OCTOPUS_HG37_BOWTIE2(bowtie2_hg37_bams)
    OCTOPUS_HG38_BOWTIE2(bowtie2_hg38_bams)
    OCTOPUS_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run Octopus variant calling on BWA-MEM alignments
    OCTOPUS_HG37_BWAMEM(bwamem_hg37_bams)
    OCTOPUS_HG38_BWAMEM(bwamem_hg38_bams)
    OCTOPUS_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH MANTA
    // ========================================

    // Run Manta variant calling on Bowtie2 alignments
    MANTA_HG37_BOWTIE2(bowtie2_hg37_bams)
    MANTA_HG38_BOWTIE2(bowtie2_hg38_bams)
    MANTA_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run Manta variant calling on BWA-MEM alignments
    MANTA_HG37_BWAMEM(bwamem_hg37_bams)
    MANTA_HG38_BWAMEM(bwamem_hg38_bams)
    MANTA_T2T_BWAMEM(bwamem_t2t_bams)
}

// Print completion message
workflow.onComplete {
    log.info """
    ================================================================
                    Pipeline completed successfully!
    ================================================================
    Results are saved in: ${params.outdir}/

    PREPROCESSING (3 tools × 1 sample = 3 outputs):
      - Fastp: ${params.outdir}/01_fastp/
      - Cutadapt: ${params.outdir}/02_cutadapt/
      - Trimmomatic: ${params.outdir}/03_trimmomatic/

    ALIGNMENT - BOWTIE2 (3 QC × 3 genomes = 9 combinations):
      
      FASTP alignments:
        - hg37: ${params.outdir}/04_bowtie2/fastp/hg37/
        - hg38: ${params.outdir}/04_bowtie2/fastp/hg38/
        - T2T:  ${params.outdir}/04_bowtie2/fastp/t2t/
      
      CUTADAPT alignments:
        - hg37: ${params.outdir}/04_bowtie2/cutadapt/hg37/
        - hg38: ${params.outdir}/04_bowtie2/cutadapt/hg38/
        - T2T:  ${params.outdir}/04_bowtie2/cutadapt/t2t/
      
      TRIMMOMATIC alignments:
        - hg37: ${params.outdir}/04_bowtie2/trimmomatic/hg37/
        - hg38: ${params.outdir}/04_bowtie2/trimmomatic/hg38/
        - T2T:  ${params.outdir}/04_bowtie2/trimmomatic/t2t/

    ALIGNMENT - BWA-MEM (3 QC × 3 genomes = 9 combinations):
      
      FASTP alignments:
        - hg37: ${params.outdir}/05_bwamem/fastp/hg37/
        - hg38: ${params.outdir}/05_bwamem/fastp/hg38/
        - T2T:  ${params.outdir}/05_bwamem/fastp/t2t/
      
      CUTADAPT alignments:
        - hg37: ${params.outdir}/05_bwamem/cutadapt/hg37/
        - hg38: ${params.outdir}/05_bwamem/cutadapt/hg38/
        - T2T:  ${params.outdir}/05_bwamem/cutadapt/t2t/
      
      TRIMMOMATIC alignments:
        - hg37: ${params.outdir}/05_bwamem/trimmomatic/hg37/
        - hg38: ${params.outdir}/05_bwamem/trimmomatic/hg38/
        - T2T:  ${params.outdir}/05_bwamem/trimmomatic/t2t/

    VARIANT CALLING - OCTOPUS (2 aligners × 3 QC × 3 genomes = 18 combinations):

      BOWTIE2 variant calls:
        FASTP:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bowtie2/fastp/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bowtie2/fastp/hg38/
          - T2T:  ${params.outdir}/06_variant_calling/octopus/bowtie2/fastp/t2t/
        
        CUTADAPT:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bowtie2/cutadapt/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bowtie2/cutadapt/hg38/
          - T2T:  ${params.outdir}/06_variant_calling/octopus/bowtie2/cutadapt/t2t/
        
        TRIMMOMATIC:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bowtie2/trimmomatic/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bowtie2/trimmomatic/hg38/
          - T2T:  ${params.outdir}/06_variant_calling/octopus/bowtie2/trimmomatic/t2t/

      BWA-MEM variant calls:
        FASTP:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bwamem/fastp/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bwamem/fastp/hg38/
          - T2T:  ${params.outdir}/06_variant_calling/octopus/bwamem/fastp/t2t/
        
        CUTADAPT:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bwamem/cutadapt/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bwamem/cutadapt/hg38/
          - T2T:  ${params.outdir}/06_variant_calling/octopus/bwamem/cutadapt/t2t/
        
        TRIMMOMATIC:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bwamem/trimmomatic/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bwamem/trimmomatic/hg38/
          - T2T:  ${params.outdir}/06_variant_calling/octopus/bwamem/trimmomatic/t2t/

    TOTAL OUTPUT FILES: 30 (3 preprocessing + 9 bowtie2 + 9 bwamem + 9 octopus/bowtie2 + 9 octopus/bwamem)
    ================================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution failed: ${workflow.errorMessage}"
}
