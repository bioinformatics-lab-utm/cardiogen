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
// include { CUTADAPT } from './modules/cutadapt'
// include { TRIMMOMATIC } from './modules/trimmomatic'
include { FASTQC } from './modules/fastqc'
include { FASTQC_FASTP } from './modules/fastqc_fastp'
// include { FASTQC_CUTADAPT } from './modules/fastqc_cutadapt'
// include { FASTQC_TRIMMOMATIC } from './modules/fastqc_trimmomatic'
// include { BOWTIE2_ALIGN_HG37 } from './modules/bowtie2'
include { BOWTIE2_ALIGN_HG38 } from './modules/bowtie2'
// include { BOWTIE2_ALIGN_T2T } from './modules/bowtie2'
// include { BWAMEM_HG37 } from './modules/bwa'
include { BWAMEM_HG38 } from './modules/bwa'
// include { BWAMEM_T2T } from './modules/bwa'
// include { OCTOPUS_HG37_BOWTIE2 } from './modules/octopus'
// include { OCTOPUS_HG38_BOWTIE2 } from './modules/octopus'
// include { OCTOPUS_HG37_BWAMEM } from './modules/octopus'
// include { OCTOPUS_HG38_BWAMEM } from './modules/octopus'
// include { MANTA_HG37_BOWTIE2 } from './modules/manta'
include { MANTA_HG38_BOWTIE2 } from './modules/manta'
// include { MANTA_T2T_BOWTIE2 } from './modules/manta'
// include { MANTA_HG37_BWAMEM } from './modules/manta'
include { MANTA_HG38_BWAMEM } from './modules/manta'
// include { MANTA_T2T_BWAMEM } from './modules/manta'
// include { DELLY_HG37_BOWTIE2 } from './modules/delly'
include { DELLY_HG38_BOWTIE2 } from './modules/delly'
// include { DELLY_T2T_BOWTIE2 } from './modules/delly'
// include { DELLY_HG37_BWAMEM } from './modules/delly'
include { DELLY_HG38_BWAMEM } from './modules/delly'
// include { DELLY_T2T_BWAMEM } from './modules/delly'
// include { DEEPVARIANT_HG37_BOWTIE2 } from './modules/deepvariant'
include { DEEPVARIANT_HG38_BOWTIE2 } from './modules/deepvariant'
// include { DEEPVARIANT_T2T_BOWTIE2 } from './modules/deepvariant'
// include { DEEPVARIANT_HG37_BWAMEM } from './modules/deepvariant'
include { DEEPVARIANT_HG38_BWAMEM } from './modules/deepvariant'
// include { DEEPVARIANT_T2T_BWAMEM } from './modules/deepvariant'
// include { GATK_HAPLOTYPECALLER_HG37_BOWTIE2 } from './modules/gatk'
include { GATK_HAPLOTYPECALLER_HG38_BOWTIE2 } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_T2T_BOWTIE2 } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_HG37_BWAMEM } from './modules/gatk'
include { GATK_HAPLOTYPECALLER_HG38_BWAMEM } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_T2T_BWAMEM } from './modules/gatk'
include { HAPPY_COMPARE } from './modules/happy'
include { HAPPY_COMPARE as HAPPY_COMPARE_COMBINED } from './modules/happy'
// include { TABIX_INDEX as TABIX_OCTOPUS_HG38_BOWTIE2 } from './modules/tabix'
// include { TABIX_INDEX as TABIX_OCTOPUS_HG38_BWAMEM } from './modules/tabix'

// Parameters
params.input_dir = "${projectDir}/test_data/ont_data"
params.outdir = "${projectDir}/results"
params.pattern = "HG002_R{1,2}_complete.fastq.gz"

// Reference genome parameters
// params.hg37_index = "${projectDir}/reference/hg37"
params.hg38_index = "${projectDir}/reference/hg38"
// params.t2t_index = "${projectDir}/reference/t2t"

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
    // --hg37_index      Path to hg37 bowtie2 index directory (default: ${params.hg37_index})
    --hg38_index      Path to hg38 bowtie2 index directory (default: ${params.hg38_index})
    // --t2t_index       Path to T2T bowtie2 index directory (default: ${params.t2t_index})
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
    // CUTADAPT(fastq_ch)
    // TRIMMOMATIC(fastq_ch)

    // Run FastQC on raw data
    FASTQC(fastq_ch)

    // Run FastQC on processed data
    FASTQC_FASTP(FASTP.out.reads)
    // FASTQC_CUTADAPT(CUTADAPT.out.reads)
    // FASTQC_TRIMMOMATIC(TRIMMOMATIC.out.paired_reads)

    // Prepare metadata for alignments with each genome
    
    // FASTP outputs
    // fastp_hg37 = FASTP.out.reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "fastp"
    //         [meta, tuple[1]]
    //     }
    
    fastp_hg38 = FASTP.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "fastp"
            [meta, tuple[1]]
        }
    
    // fastp_t2t = FASTP.out.reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "fastp"
    //         [meta, tuple[1]]
    //     }

    // CUTADAPT outputs
    // cutadapt_hg37 = CUTADAPT.out.reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "cutadapt"
    //         [meta, tuple[1]]
    //     }
    
    // cutadapt_hg38 = CUTADAPT.out.reads
    //     .map { tuple ->
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "cutadapt"
    //         [meta, tuple[1]]
    //     }

    // cutadapt_t2t = CUTADAPT.out.reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "cutadapt"
    //         [meta, tuple[1]]
    //     }

    // TRIMMOMATIC outputs
    // trimmomatic_hg37 = TRIMMOMATIC.out.paired_reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "trimmomatic"
    //         [meta, tuple[1]]
    //     }
    
    // trimmomatic_hg38 = TRIMMOMATIC.out.paired_reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "trimmomatic"
    //         [meta, tuple[1]]
    //     }
    
    // trimmomatic_t2t = TRIMMOMATIC.out.paired_reads
    //     .map { tuple -> 
    //         def meta = tuple[0].clone()
    //         meta.qc_tool = "trimmomatic"
    //         [meta, tuple[1]]
    //     }

    // Combine all inputs for each genome and run alignments
    // hg37_inputs = fastp_hg37 // .mix(cutadapt_hg37, trimmomatic_hg37)
    hg38_inputs = fastp_hg38 // .mix(cutadapt_hg38) // .mix(trimmomatic_hg38)
    // t2t_inputs = fastp_t2t // .mix(cutadapt_t2t, trimmomatic_t2t)

    // Run Bowtie2 alignments
    // BOWTIE2_ALIGN_HG37(hg37_inputs)
    BOWTIE2_ALIGN_HG38(hg38_inputs)
    // BOWTIE2_ALIGN_T2T(t2t_inputs)

    // Run BWA-MEM alignments
    // BWAMEM_HG37(hg37_inputs)
    BWAMEM_HG38(hg38_inputs)
    // BWAMEM_T2T(t2t_inputs)

    // ========================================
    // VARIANT CALLING WITH OCTOPUS
    // ========================================

    // Add aligner metadata to BAM files from Bowtie2
    // bowtie2_hg37_bams = BOWTIE2_ALIGN_HG37.out.bam
    //     .join(BOWTIE2_ALIGN_HG37.out.bai)
    //     .map { meta, bam, bai ->
    //         def new_meta = meta.clone()
    //         new_meta.aligner = "bowtie2"
    //         [new_meta, bam, bai]
    //     }

    bowtie2_hg38_bams = BOWTIE2_ALIGN_HG38.out.bam
        .join(BOWTIE2_ALIGN_HG38.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bowtie2"
            [new_meta, bam, bai]
        }

    // bowtie2_t2t_bams = BOWTIE2_ALIGN_T2T.out.bam
    //     .join(BOWTIE2_ALIGN_T2T.out.bai)
    //     .map { meta, bam, bai ->
    //         def new_meta = meta.clone()
    //         new_meta.aligner = "bowtie2"
    //         [new_meta, bam, bai]
    //     }

    // Add aligner metadata to BAM files from BWA-MEM
    // bwamem_hg37_bams = BWAMEM_HG37.out.bam
    //     .join(BWAMEM_HG37.out.bai)
    //     .map { meta, bam, bai ->
    //         def new_meta = meta.clone()
    //         new_meta.aligner = "bwamem"
    //         [new_meta, bam, bai]
    //     }

    bwamem_hg38_bams = BWAMEM_HG38.out.bam
        .join(BWAMEM_HG38.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bwamem"
            [new_meta, bam, bai]
        }

    // bwamem_t2t_bams = BWAMEM_T2T.out.bam
    //     .join(BWAMEM_T2T.out.bai)
    //     .map { meta, bam, bai ->
    //         def new_meta = meta.clone()
    //         new_meta.aligner = "bwamem"
    //         [new_meta, bam, bai]
    //     }

    // Run Octopus variant calling on Bowtie2 alignments
    // OCTOPUS_HG37_BOWTIE2(bowtie2_hg37_bams)
    // OCTOPUS_HG38_BOWTIE2(bowtie2_hg38_bams)

    // Run Octopus variant calling on BWA-MEM alignments
    // OCTOPUS_HG37_BWAMEM(bwamem_hg37_bams)
    // OCTOPUS_HG38_BWAMEM(bwamem_hg38_bams)

    // ========================================
    // VARIANT CALLING WITH MANTA
    // ========================================

    // Run Manta variant calling on Bowtie2 alignments
    // MANTA_HG37_BOWTIE2(bowtie2_hg37_bams)
    MANTA_HG38_BOWTIE2(bowtie2_hg38_bams)
    // MANTA_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run Manta variant calling on BWA-MEM alignments
    // MANTA_HG37_BWAMEM(bwamem_hg37_bams)
    MANTA_HG38_BWAMEM(bwamem_hg38_bams)
    // MANTA_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH DELLY
    // ========================================

    // Run Delly variant calling on Bowtie2 alignments
    // DELLY_HG37_BOWTIE2(bowtie2_hg37_bams)
    DELLY_HG38_BOWTIE2(bowtie2_hg38_bams)
    // DELLY_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run Delly variant calling on BWA-MEM alignments
    // DELLY_HG37_BWAMEM(bwamem_hg37_bams)
    DELLY_HG38_BWAMEM(bwamem_hg38_bams)
    // DELLY_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH DEEPVARIANT
    // ========================================

    // Run DeepVariant variant calling on Bowtie2 alignments
    // DEEPVARIANT_HG37_BOWTIE2(bowtie2_hg37_bams)
    DEEPVARIANT_HG38_BOWTIE2(bowtie2_hg38_bams)
    // DEEPVARIANT_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run DeepVariant variant calling on BWA-MEM alignments
    // DEEPVARIANT_HG37_BWAMEM(bwamem_hg37_bams)
    DEEPVARIANT_HG38_BWAMEM(bwamem_hg38_bams)
    // DEEPVARIANT_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH GATK HAPLOTYPECALLER
    // ========================================

    // Run GATK HaplotypeCaller variant calling on Bowtie2 alignments
    // GATK_HAPLOTYPECALLER_HG37_BOWTIE2(bowtie2_hg37_bams)
    GATK_HAPLOTYPECALLER_HG38_BOWTIE2(bowtie2_hg38_bams)
    // GATK_HAPLOTYPECALLER_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run GATK HaplotypeCaller variant calling on BWA-MEM alignments
    // GATK_HAPLOTYPECALLER_HG37_BWAMEM(bwamem_hg37_bams)
    GATK_HAPLOTYPECALLER_HG38_BWAMEM(bwamem_hg38_bams)
    // GATK_HAPLOTYPECALLER_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // BENCHMARKING WITH hap.py (HG38 only)
    // ========================================

    // Load truth VCF file and reference for HG38
    def truth_hg38_vcf = file(params.truth_hg38_vcf)
    def truth_hg38_vcf_idx = file("${params.truth_hg38_vcf}.tbi")
    def truth_hg38_bed = file(params.truth_hg38_bed)
    def hg38_ref = file("${params.hg38_index}/hg38.fa")
    def hg38_ref_fai = file("${params.hg38_index}/hg38.fa.fai")

    // Index Octopus VCFs (they don't output .tbi)
    // TABIX_OCTOPUS_HG38_BOWTIE2(OCTOPUS_HG38_BOWTIE2.out.vcf)
    // TABIX_OCTOPUS_HG38_BWAMEM(OCTOPUS_HG38_BWAMEM.out.vcf)

    // Prepare ALL variant channels with QC information (18 combinations: 3 callers × 2 aligners × 3 QC)
    // DeepVariant HG38
    deepvariant_hg38_bwamem_fastp = DEEPVARIANT_HG38_BWAMEM.out.vcf
        .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
        .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
        .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "fastp", "hg38", vcf, idx] }
    
    deepvariant_hg38_bowtie2_fastp = DEEPVARIANT_HG38_BOWTIE2.out.vcf
        .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
        .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
        .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "fastp", "hg38", vcf, idx] }
    
    // deepvariant_hg38_bowtie2_cutadapt = DEEPVARIANT_HG38_BOWTIE2.out.vcf
    //     .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'cutadapt' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "cutadapt", "hg38", vcf, idx] }
    
    // deepvariant_hg38_bowtie2_trimmomatic = DEEPVARIANT_HG38_BOWTIE2.out.vcf
    //     .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'trimmomatic' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "trimmomatic", "hg38", vcf, idx] }
    
    // deepvariant_hg38_bwamem_cutadapt = DEEPVARIANT_HG38_BWAMEM.out.vcf
    //     .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'cutadapt' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "cutadapt", "hg38", vcf, idx] }
    
    // deepvariant_hg38_bwamem_trimmomatic = DEEPVARIANT_HG38_BWAMEM.out.vcf
    //     .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'trimmomatic' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "trimmomatic", "hg38", vcf, idx] }

    // GATK HG38
    gatk_hg38_bwamem_fastp = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
        .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
        .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
        .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "fastp", "hg38", vcf, idx] }
    
    gatk_hg38_bowtie2_fastp = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
        .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
        .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
        .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "fastp", "hg38", vcf, idx] }
    
    // gatk_hg38_bowtie2_cutadapt = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'cutadapt' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "cutadapt", "hg38", vcf, idx] }
    
    // gatk_hg38_bowtie2_trimmomatic = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'trimmomatic' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "trimmomatic", "hg38", vcf, idx] }
    
    // gatk_hg38_bwamem_cutadapt = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'cutadapt' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "cutadapt", "hg38", vcf, idx] }
    
    // gatk_hg38_bwamem_trimmomatic = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'trimmomatic' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "trimmomatic", "hg38", vcf, idx] }

    // Octopus HG38 (using indexed VCFs)
    // octopus_hg38_bowtie2_fastp = TABIX_OCTOPUS_HG38_BOWTIE2.out.indexed_vcf
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bowtie2", "fastp", "hg38", vcf, idx] }
    
    // octopus_hg38_bowtie2_cutadapt = TABIX_OCTOPUS_HG38_BOWTIE2.out.indexed_vcf
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'cutadapt' }
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bowtie2", "cutadapt", "hg38", vcf, idx] }
    
    // octopus_hg38_bowtie2_trimmomatic = TABIX_OCTOPUS_HG38_BOWTIE2.out.indexed_vcf
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'trimmomatic' }
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bowtie2", "trimmomatic", "hg38", vcf, idx] }
    
    // octopus_hg38_bwamem_fastp = TABIX_OCTOPUS_HG38_BWAMEM.out.indexed_vcf
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bwamem", "fastp", "hg38", vcf, idx] }
    
    // octopus_hg38_bwamem_cutadapt = TABIX_OCTOPUS_HG38_BWAMEM.out.indexed_vcf
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'cutadapt' }
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bwamem", "cutadapt", "hg38", vcf, idx] }
    
    // octopus_hg38_bwamem_trimmomatic = TABIX_OCTOPUS_HG38_BWAMEM.out.indexed_vcf
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'trimmomatic' }
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bwamem", "trimmomatic", "hg38", vcf, idx] }

    // Combine all HG38 variants for separate comparisons
    all_hg38_variants = deepvariant_hg38_bwamem_fastp.mix(
        deepvariant_hg38_bowtie2_fastp,
        // deepvariant_hg38_bowtie2_cutadapt,
        // deepvariant_hg38_bowtie2_trimmomatic,
        // deepvariant_hg38_bwamem_cutadapt,
        // deepvariant_hg38_bwamem_trimmomatic,
        gatk_hg38_bowtie2_fastp,
        // gatk_hg38_bowtie2_cutadapt,
        // gatk_hg38_bowtie2_trimmomatic,
        // gatk_hg38_bwamem_cutadapt,
        // gatk_hg38_bwamem_trimmomatic,
        // octopus_hg38_bowtie2_fastp,
        // octopus_hg38_bowtie2_cutadapt,
        // octopus_hg38_bowtie2_trimmomatic,
        // octopus_hg38_bwamem_fastp,
        // octopus_hg38_bwamem_cutadapt,
        // octopus_hg38_bwamem_trimmomatic
        gatk_hg38_bwamem_fastp
    )

    // Run hap.py comparison for each HG38 variant (18 separate comparisons)
    HAPPY_COMPARE(
        all_hg38_variants,
        truth_hg38_vcf,
        truth_hg38_vcf_idx,
        truth_hg38_bed,
        hg38_ref,
        hg38_ref_fai
    )

    // ========================================
    // COMBINED COMPARISON (all QC tools mixed)
    // ========================================

    // Prepare combined channels (without QC distinction) - 6 comparisons
    deepvariant_hg38_bwamem_combined = DEEPVARIANT_HG38_BWAMEM.out.vcf
        .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
        .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "combined", "hg38", vcf, idx] }

    gatk_hg38_bwamem_combined = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
        .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
        .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "combined", "hg38", vcf, idx] }

    deepvariant_hg38_bowtie2_combined = DEEPVARIANT_HG38_BOWTIE2.out.vcf
        .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
        .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "combined", "hg38", vcf, idx] }

    gatk_hg38_bowtie2_combined = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
        .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
        .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "combined", "hg38", vcf, idx] }

    // octopus_hg38_bowtie2_combined = TABIX_OCTOPUS_HG38_BOWTIE2.out.indexed_vcf
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bowtie2", "combined", "hg38", vcf, idx] }
    
    // octopus_hg38_bwamem_combined = TABIX_OCTOPUS_HG38_BWAMEM.out.indexed_vcf
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bwamem", "combined", "hg38", vcf, idx] }

    // Mix all for combined comparison
    all_hg38_variants_combined = deepvariant_hg38_bwamem_combined.mix(
        deepvariant_hg38_bowtie2_combined,
        gatk_hg38_bowtie2_combined,
        // octopus_hg38_bowtie2_combined,
        // octopus_hg38_bwamem_combined
        gatk_hg38_bwamem_combined
    )

    // Run hap.py for combined variants (6 comparisons: 3 callers × 2 aligners)
    HAPPY_COMPARE_COMBINED(
        all_hg38_variants_combined,
        truth_hg38_vcf,
        truth_hg38_vcf_idx,
        truth_hg38_bed,
        hg38_ref,
        hg38_ref_fai
    )
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
        
        CUTADAPT:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bowtie2/cutadapt/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bowtie2/cutadapt/hg38/
        
        TRIMMOMATIC:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bowtie2/trimmomatic/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bowtie2/trimmomatic/hg38/

      BWA-MEM variant calls:
        FASTP:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bwamem/fastp/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bwamem/fastp/hg38/
        
        CUTADAPT:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bwamem/cutadapt/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bwamem/cutadapt/hg38/
        
        TRIMMOMATIC:
          - hg37: ${params.outdir}/06_variant_calling/octopus/bwamem/trimmomatic/hg37/
          - hg38: ${params.outdir}/06_variant_calling/octopus/bwamem/trimmomatic/hg38/

    VARIANT CALLING - MANTA (2 aligners × 3 QC × 3 genomes = 18 combinations):
      Output directories similar to Octopus under:
        ${params.outdir}/06_variant_calling/manta/

    VARIANT CALLING - DELLY (2 aligners × 3 QC × 3 genomes = 18 combinations):
      Output directories similar to Octopus under:
        ${params.outdir}/06_variant_calling/delly/

    VARIANT CALLING - DEEPVARIANT (2 aligners × 3 QC × 3 genomes = 18 combinations):
      Output directories similar to Octopus under:
        ${params.outdir}/06_variant_calling/deepvariant/

    VARIANT CALLING - GATK HAPLOTYPECALLER (2 aligners × 3 QC × 3 genomes = 18 combinations):
      Output directories similar to Octopus under:
        ${params.outdir}/06_variant_calling/gatk/
    ================================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution failed: ${workflow.errorMessage}"
}
