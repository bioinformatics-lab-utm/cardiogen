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
// include { BOWTIE2_ALIGN_HG38 } from './modules/bowtie2'
// include { BOWTIE2_ALIGN_T2T } from './modules/bowtie2'
// include { BWAMEM_HG37 } from './modules/bwa'
include { BWAMEM_HG38 } from './modules/bwa'
// include { BWAMEM_T2T } from './modules/bwa'
include { COUNT_READS } from './modules/count_reads'
// include { OCTOPUS_HG37_BOWTIE2 } from './modules/octopus'
// include { OCTOPUS_HG38_BOWTIE2 } from './modules/octopus'
// include { OCTOPUS_HG37_BWAMEM } from './modules/octopus'
// include { OCTOPUS_HG38_BWAMEM } from './modules/octopus'
// include { MANTA_HG37_BOWTIE2 } from './modules/manta'
// include { MANTA_HG38_BOWTIE2 } from './modules/manta'
// include { MANTA_T2T_BOWTIE2 } from './modules/manta'
// include { MANTA_HG37_BWAMEM } from './modules/manta'
include { MANTA_HG38_BWAMEM } from './modules/manta'
// include { MANTA_T2T_BWAMEM } from './modules/manta'
// include { DELLY_HG37_BOWTIE2 } from './modules/delly'
// include { DELLY_HG38_BOWTIE2 } from './modules/delly'
// include { DELLY_T2T_BOWTIE2 } from './modules/delly'
// include { DELLY_HG37_BWAMEM } from './modules/delly'
include { DELLY_HG38_BWAMEM } from './modules/delly'
// include { DELLY_T2T_BWAMEM } from './modules/delly'
// include { DEEPVARIANT_HG37_BOWTIE2 } from './modules/deepvariant'
// include { DEEPVARIANT_HG38_BOWTIE2 } from './modules/deepvariant'
// include { DEEPVARIANT_T2T_BOWTIE2 } from './modules/deepvariant'
// include { DEEPVARIANT_HG37_BWAMEM } from './modules/deepvariant'
include { DEEPVARIANT_HG38_BWAMEM } from './modules/deepvariant'
// include { DEEPVARIANT_T2T_BWAMEM } from './modules/deepvariant'
// include { GATK_HAPLOTYPECALLER_HG37_BOWTIE2 } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_HG38_BOWTIE2 } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_T2T_BOWTIE2 } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_HG37_BWAMEM } from './modules/gatk'
include { GATK_HAPLOTYPECALLER_HG38_BWAMEM } from './modules/gatk'
// include { GATK_HAPLOTYPECALLER_T2T_BWAMEM } from './modules/gatk'
include { CROSSMAP_LIFTOVER } from './modules/liftover'
include { HAPPY_COMPARE } from './modules/happy'
include { BCFTOOLS_ISEC } from './modules/bcftools_isec'
// Truvari is intentionally NOT wired in: the reference VCFs (reference/vcf_hdd) contain
// only SNPs/indels (no structural variants), so there is no SV truth for Manta/Delly to
// be benchmarked against.
// include { TRUVARI_COMPARE } from './modules/truvari'
// include { BCF_TO_VCFGZ } from './modules/tabix'  // Not needed - Delly not benchmarked with Happy
// include { COMBINE_SNP_INDEL } from './modules/combine_callers'
// include { COMBINE_SV } from './modules/combine_callers'
// include { TABIX_INDEX as TABIX_OCTOPUS_HG38_BOWTIE2 } from './modules/tabix'
// include { TABIX_INDEX as TABIX_OCTOPUS_HG38_BWAMEM } from './modules/tabix'

// Parameters
// `input_dir` holds one sub-folder per sequencing run; each run folder contains many
// paired-end fastq.gz files. Results are mirrored per-run using meta.run (see below).
params.input_dir = "${projectDir}/data"
params.outdir = "${projectDir}/results"
params.pattern = "*_R{1,2}_001.fastq.gz"

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
    --hg38_index      Path to hg38 bwa index directory (default: ${params.hg38_index})
    // --t2t_index       Path to T2T bowtie2 index directory (default: ${params.t2t_index})
    
    Exome Sequencing Options:
    --exome_mode      Enable exome sequencing mode (default: ${params.exome_mode})
    --exome_bed_hg37  Path to hg37 exome target regions BED file (default: ${params.exome_bed_hg37})
    --exome_bed_hg38  Path to hg38 exome target regions BED file (default: ${params.exome_bed_hg38})
    --exome_padding   Padding (bp) around exome regions (default: ${params.exome_padding})
    --seq_platform    Sequencing platform for DeepVariant model (default: ${params.seq_platform})
                      Options: WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA
    
    --help            Show this help message
    
    Example (WGS):
    nextflow run main.nf --input_dir /path/to/fastq --outdir /path/to/results
    
    Example (Exome):
    nextflow run main.nf --input_dir /path/to/fastq --outdir /path/to/results \\
        --exome_mode true --seq_platform WES \\
        --exome_bed_hg38 /path/to/exome_regions.bed
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
    Exome mode      : ${params.exome_mode}
    Seq platform    : ${params.seq_platform}
    ================================================================
    """.stripIndent()

    // Create input channel from FASTQ files.
    // Files live in `${input_dir}/<run>/<sample>_R{1,2}_001.fastq.gz`. The grouping
    // key embeds the run folder so identically-named samples in different runs never
    // get paired together, and meta.run is used to mirror every result per run folder.
    fastq_ch = Channel
        .fromFilePairs("${params.input_dir}/*/${params.pattern}", checkIfExists: true) { file ->
            def run = file.parent.name
            def sample = file.name.replaceAll(/_R[12]_001\.fastq\.gz$/, '')
            "${run}__${sample}"
        }
        .map { key, files ->
            def parts = key.split('__', 2)
            def meta = [:]
            meta.run = parts[0]
            meta.id = parts[1]
            meta.single_end = false  // Paired-end data
            [meta, files]
        }

    // Step 1: Run FastQC on raw data (quality control first)
    FASTQC(fastq_ch)

    // Step 2: Run preprocessing tools on raw data
    FASTP(fastq_ch)
    // CUTADAPT(fastq_ch)
    // TRIMMOMATIC(fastq_ch)

    // Step 3: Run FastQC on processed data
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
    // BOWTIE2_ALIGN_HG38(hg38_inputs)
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

    // bowtie2_hg38_bams = BOWTIE2_ALIGN_HG38.out.bam
    //     .join(BOWTIE2_ALIGN_HG38.out.bai)
    //     .map { meta, bam, bai ->
    //         def new_meta = meta.clone()
    //         new_meta.aligner = "bowtie2"
    //         [new_meta, bam, bai]
    //     }

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

    // Create channel with BAM files and broadcast to multiple variant callers
    bwamem_hg38_bams = BWAMEM_HG38.out.bam
        .join(BWAMEM_HG38.out.bai)
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.aligner = "bwamem"
            [new_meta, bam, bai]
        }

    // Gate: drop empty/failed samples before variant calling. Callers such as Delly
    // abort ("Sample has not enough data to estimate library parameters!") on near-empty
    // BAMs, and because the global errorStrategy is 'finish', a single junk sample (e.g. a
    // negative control or a sample that failed sequencing - only a handful of reads) would
    // otherwise kill the entire run. COUNT_READS reports the number of properly-paired
    // mapped reads; samples below params.min_mapped_reads are logged and excluded here,
    // while genuine failures on real samples still surface normally.
    bwamem_hg38_bams_gated = COUNT_READS(bwamem_hg38_bams).counted
        .branch { meta, bam, bai, count ->
            pass: count.toInteger() >= params.min_mapped_reads
            skip: true
        }

    bwamem_hg38_bams_gated.skip.subscribe { meta, bam, bai, count ->
        log.warn "Skipping variant calling for ${meta.run} / ${meta.id} (${meta.qc_tool}): " +
                 "only ${count} properly-paired mapped reads (< params.min_mapped_reads=${params.min_mapped_reads})"
    }

    // Using multiMap to create separate copies for each variant caller
    bwamem_hg38_bams_all = bwamem_hg38_bams_gated.pass
        .map { meta, bam, bai, count -> [meta, bam, bai] }
        .multiMap { meta, bam, bai ->
            manta: [meta, bam, bai]
            delly: [meta, bam, bai]
            deepvariant: [meta, bam, bai]
            gatk: [meta, bam, bai]
        }
    
    // Assign each channel copy
    bwamem_hg38_bams_manta = bwamem_hg38_bams_all.manta
    bwamem_hg38_bams_delly = bwamem_hg38_bams_all.delly
    bwamem_hg38_bams_deepvariant = bwamem_hg38_bams_all.deepvariant
    bwamem_hg38_bams_gatk = bwamem_hg38_bams_all.gatk

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
    // MANTA_HG38_BOWTIE2(bowtie2_hg38_bams)
    // MANTA_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run Manta variant calling on BWA-MEM alignments
    // MANTA_HG37_BWAMEM(bwamem_hg37_bams)
    exome_bed_ch = params.exome_mode && params.exome_bed_hg38 ? Channel.value(file(params.exome_bed_hg38)) : Channel.value(file('NO_FILE'))
    exome_bed_tbi_ch = params.exome_mode && params.exome_bed_hg38 ? Channel.value(file(params.exome_bed_hg38 + '.tbi')) : Channel.value(file('NO_FILE'))
    MANTA_HG38_BWAMEM(bwamem_hg38_bams_manta, exome_bed_ch, exome_bed_tbi_ch)
    // MANTA_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH DELLY
    // ========================================

    // Run Delly variant calling on Bowtie2 alignments
    // DELLY_HG37_BOWTIE2(bowtie2_hg37_bams)
    // DELLY_HG38_BOWTIE2(bowtie2_hg38_bams)
    // DELLY_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run Delly variant calling on BWA-MEM alignments
    // DELLY_HG37_BWAMEM(bwamem_hg37_bams)
    DELLY_HG38_BWAMEM(bwamem_hg38_bams_delly)
    // DELLY_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH DEEPVARIANT
    // ========================================

    // Run DeepVariant variant calling on Bowtie2 alignments
    // DEEPVARIANT_HG37_BOWTIE2(bowtie2_hg37_bams)
    // DEEPVARIANT_HG38_BOWTIE2(bowtie2_hg38_bams)
    // DEEPVARIANT_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run DeepVariant variant calling on BWA-MEM alignments
    // DEEPVARIANT_HG37_BWAMEM(bwamem_hg37_bams)
    DEEPVARIANT_HG38_BWAMEM(bwamem_hg38_bams_deepvariant, exome_bed_ch, exome_bed_tbi_ch)
    // DEEPVARIANT_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // VARIANT CALLING WITH GATK HAPLOTYPECALLER
    // ========================================

    // Run GATK HaplotypeCaller variant calling on Bowtie2 alignments
    // GATK_HAPLOTYPECALLER_HG37_BOWTIE2(bowtie2_hg37_bams)
    // GATK_HAPLOTYPECALLER_HG38_BOWTIE2(bowtie2_hg38_bams)
    // GATK_HAPLOTYPECALLER_T2T_BOWTIE2(bowtie2_t2t_bams)

    // Run GATK HaplotypeCaller variant calling on BWA-MEM alignments
    // GATK_HAPLOTYPECALLER_HG37_BWAMEM(bwamem_hg37_bams)
    GATK_HAPLOTYPECALLER_HG38_BWAMEM(bwamem_hg38_bams_gatk, exome_bed_ch, exome_bed_tbi_ch)
    // GATK_HAPLOTYPECALLER_T2T_BWAMEM(bwamem_t2t_bams)

    // ========================================
    // COMPARISON: PIPELINE CALLS vs LABORATORY REFERENCE VCFs (SNP/INDEL)
    // ========================================
    // The laboratory reference VCFs live in reference/vcf_hdd/<run>/<sample>.vcf. They are
    // GRCh37/hg19 and contain SNP/indel calls only, so they are:
    //   1. matched to the pipeline calls per-sample (by run + sample, stripping the fastq
    //      lane suffix: e.g. meta.id "001TC_S1_L001" -> reference "001TC_S1.vcf"),
    //   2. lifted hg19 -> hg38 with CrossMap (CROSSMAP_LIFTOVER), then
    //   3. compared against the hg38 pipeline SNP/indel calls (DeepVariant, GATK) with
    //      both hap.py (HAPPY_COMPARE) and bcftools isec (BCFTOOLS_ISEC).
    // Samples without a matching reference VCF (e.g. empty run folders) are dropped by the
    // .exists() filter and simply not compared.

    // Strip the Illumina lane suffix (_L001, ...) so meta.id matches the reference basename
    def strip_lane = { id -> id.replaceAll(/_L0*\d+$/, '') }

    // SNP/indel pipeline calls, keyed by "<run>__<sample>" for matching to the reference
    dv_snpindel = DEEPVARIANT_HG38_BWAMEM.out.vcf
        .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
        .map { meta, vcf, idx -> ["${meta.run}__${strip_lane(meta.id)}", meta, 'deepvariant', vcf, idx] }

    gatk_snpindel = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
        .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
        .map { meta, vcf, idx -> ["${meta.run}__${strip_lane(meta.id)}", meta, 'gatk', vcf, idx] }

    snpindel_calls = dv_snpindel.mix(gatk_snpindel)

    // One reference VCF per run+sample (deduped; DeepVariant and GATK share the same one),
    // keeping only those that actually exist on disk.
    reference_vcfs = snpindel_calls
        .map { key, meta, caller, vcf, idx ->
            [key, meta.run, strip_lane(meta.id),
             file("${projectDir}/reference/vcf_hdd/${meta.run}/${strip_lane(meta.id)}.vcf")]
        }
        .filter { key, run, sample, ref_vcf -> ref_vcf.exists() }
        .unique { it[0] }

    // Lift the reference VCFs hg19 -> hg38
    CROSSMAP_LIFTOVER(reference_vcfs)

    // Pair each pipeline call with its lifted reference (by run+sample key)
    comparison_inputs = snpindel_calls
        .combine(CROSSMAP_LIFTOVER.out.lifted, by: 0)
        .map { key, meta, caller, vcf, idx, ref_vcf, ref_tbi ->
            [meta.id, caller, meta.aligner, meta.qc_tool, 'hg38', vcf, idx, ref_vcf, ref_tbi]
        }

    // Run both comparison methods
    HAPPY_COMPARE(comparison_inputs)
    BCFTOOLS_ISEC(comparison_inputs)

    // ========================================
    // COMBINE VARIANT CALLERS - COMMENTED OUT
    // ========================================

    // Combine GATK + DeepVariant for SNP/INDEL (HG38, FASTP, BWAMEM) - COMMENTED OUT
    // combined_snpindel_hg38_bwamem_fastp = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .join(
    //         DEEPVARIANT_HG38_BWAMEM.out.vcf
    //             .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    //             .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' },
    //         by: 0
    //     )
    //     .map { meta, gatk_vcf, gatk_idx, dv_vcf, dv_idx ->
    //         [meta.id, meta.aligner, meta.qc_tool, 'hg38', gatk_vcf, gatk_idx, dv_vcf, dv_idx]
    //     }

    // Combine GATK + DeepVariant for SNP/INDEL (HG37, FASTP, BOWTIE2) - COMMENTED OUT
    // combined_snpindel_hg37_bowtie2_fastp = GATK_HAPLOTYPECALLER_HG37_BOWTIE2.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG37_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .join(
    //         DEEPVARIANT_HG37_BOWTIE2.out.vcf
    //             .join(DEEPVARIANT_HG37_BOWTIE2.out.tbi)
    //             .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' },
    //         by: 0
    //     )
    //     .map { meta, gatk_vcf, gatk_idx, dv_vcf, dv_idx ->
    //         [meta.id, meta.aligner, meta.qc_tool, 'hg37', gatk_vcf, gatk_idx, dv_vcf, dv_idx]
    //     }

    // Run COMBINE_SNP_INDEL process - COMMENTED OUT
    // COMBINE_SNP_INDEL(
    //     combined_snpindel_hg38_bwamem_fastp // .mix(combined_snpindel_hg37_bowtie2_fastp)
    // )

    // Combine DELLY + MANTA for SV (HG38, FASTP, BWAMEM) - COMMENTED OUT
    // combined_sv_hg38_bwamem_fastp = DELLY_HG38_BWAMEM.out.vcf
    //     .join(DELLY_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .join(
    //         MANTA_HG38_BWAMEM.out.vcf
    //             .join(MANTA_HG38_BWAMEM.out.tbi)
    //             .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' },
    //         by: 0
    //     )
    //     .map { meta, delly_vcf, delly_idx, manta_vcf, manta_idx ->
    //         [meta.id, meta.aligner, meta.qc_tool, 'hg38', delly_vcf, delly_idx, manta_vcf, manta_idx]
    //     }

    // Combine DELLY + MANTA for SV (HG37, FASTP, BOWTIE2) - COMMENTED OUT
    // combined_sv_hg37_bowtie2_fastp = DELLY_HG37_BOWTIE2.out.vcf
    //     .join(DELLY_HG37_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .join(
    //         MANTA_HG37_BOWTIE2.out.vcf
    //             .join(MANTA_HG37_BOWTIE2.out.tbi)
    //             .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' },
    //         by: 0
    //     )
    //     .map { meta, delly_vcf, delly_idx, manta_vcf, manta_idx ->
    //         [meta.id, meta.aligner, meta.qc_tool, 'hg37', delly_vcf, delly_idx, manta_vcf, manta_idx]
    //     }

    // Run COMBINE_SV process - COMMENTED OUT
    // COMBINE_SV(
    //     combined_sv_hg38_bwamem_fastp // .mix(combined_sv_hg37_bowtie2_fastp)
    // )

    // ========================================
    // BENCHMARKING WITH hap.py (HG38) - COMMENTED OUT
    // ========================================

    // Load truth VCF file and reference for HG37 - COMMENTED OUT
    // def truth_hg37_vcf = file(params.truth_hg37_vcf)
    // def truth_hg37_vcf_idx = file("${params.truth_hg37_vcf}.tbi")
    // def truth_hg37_bed = file(params.truth_hg37_bed)
    // def hg37_ref = file("${params.hg37_index}/hg19.fa")
    // def hg37_ref_fai = file("${params.hg37_index}/hg19.fa.fai")

    // Load truth VCF file and reference for HG38 - COMMENTED OUT
    // def truth_hg38_vcf = file(params.truth_hg38_vcf)
    // def truth_hg38_vcf_idx = file("${params.truth_hg38_vcf}.tbi")
    // def truth_hg38_bed = file(params.truth_hg38_bed)
    // def hg38_ref = file("${params.hg38_index}/hg38.fa")
    // def hg38_ref_fai = file("${params.hg38_index}/hg38.fa.fai")

    // Index Octopus VCFs (they don't output .tbi) - COMMENTED OUT
    // TABIX_OCTOPUS_HG37_BOWTIE2(OCTOPUS_HG37_BOWTIE2.out.vcf)
    // TABIX_OCTOPUS_HG37_BWAMEM(OCTOPUS_HG37_BWAMEM.out.vcf)

    // Prepare variant channels for HG38 with QC information - COMMENTED OUT
    // DeepVariant HG38
    // deepvariant_hg38_bwamem_fastp = DEEPVARIANT_HG38_BWAMEM.out.vcf
    //     .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "fastp", "hg38", vcf, idx] }
    
    // deepvariant_hg38_bowtie2_fastp = DEEPVARIANT_HG38_BOWTIE2.out.vcf
    //     .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "fastp", "hg38", vcf, idx] }

    // GATK HG38
    // gatk_hg38_bwamem_fastp = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "fastp", "hg38", vcf, idx] }
    
    // gatk_hg38_bowtie2_fastp = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "fastp", "hg38", vcf, idx] }

    // Prepare combined SNP/INDEL variants for Happy benchmarking (1 comparison: 1 aligner) - COMMENTED OUT
    // Using combined GATK+DeepVariant callsets
    // combined_snpindel_for_happy = COMBINE_SNP_INDEL.out.combined_vcf
    //     .map { sample_id, aligner, qc, reference, vcf, idx ->
    //         [sample_id, "combined_gatk_deepvariant", aligner, qc, reference, vcf, idx]
    //     }

    // Combine individual callers + combined callset for Happy benchmarking - COMMENTED OUT
    // Total: 3 comparisons (2 individual + 1 combined)
    // all_hg38_variants = deepvariant_hg38_bwamem_fastp.mix(
    //     // deepvariant_hg38_bowtie2_fastp,
    //     // gatk_hg38_bowtie2_fastp,
    //     gatk_hg38_bwamem_fastp,
    //     combined_snpindel_for_happy
    // )

    // Run hap.py comparison on individual AND combined SNP/INDEL callsets - COMMENTED OUT
    // 3 comparisons: DeepVariant (1), GATK (1), Combined GATK+DeepVariant (1)
    // HAPPY_COMPARE(
    //     all_hg38_variants,
    //     truth_hg38_vcf,
    //     truth_hg38_vcf_idx,
    //     truth_hg38_bed,
    //     hg38_ref,
    //     hg38_ref_fai
    // )

    // Index Octopus VCFs (they don't output .tbi) - HG38 COMMENTED OUT
    // TABIX_OCTOPUS_HG38_BOWTIE2(OCTOPUS_HG38_BOWTIE2.out.vcf)
    // TABIX_OCTOPUS_HG38_BWAMEM(OCTOPUS_HG38_BWAMEM.out.vcf)

    // Prepare ALL variant channels with QC information (18 combinations: 3 callers × 2 aligners × 3 QC) - HG38 COMMENTED OUT
    // DeepVariant HG38
    // deepvariant_hg38_bwamem_fastp = DEEPVARIANT_HG38_BWAMEM.out.vcf
    //     .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "fastp", "hg38", vcf, idx] }
    
    // deepvariant_hg38_bowtie2_fastp = DEEPVARIANT_HG38_BOWTIE2.out.vcf
    //     .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "fastp", "hg38", vcf, idx] }
    
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
    // gatk_hg38_bwamem_fastp = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "fastp", "hg38", vcf, idx] }
    
    // gatk_hg38_bowtie2_fastp = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "fastp", "hg38", vcf, idx] }
    
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
    // all_hg38_variants = deepvariant_hg38_bwamem_fastp.mix(
    //     deepvariant_hg38_bowtie2_fastp,
    //     // deepvariant_hg38_bowtie2_cutadapt,
    //     // deepvariant_hg38_bowtie2_trimmomatic,
    //     // deepvariant_hg38_bwamem_cutadapt,
    //     // deepvariant_hg38_bwamem_trimmomatic,
    //     gatk_hg38_bowtie2_fastp,
    //     // gatk_hg38_bowtie2_cutadapt,
    //     // gatk_hg38_bowtie2_trimmomatic,
    //     // gatk_hg38_bwamem_cutadapt,
    //     // gatk_hg38_bwamem_trimmomatic,
    //     // octopus_hg38_bowtie2_fastp,
    //     // octopus_hg38_bowtie2_cutadapt,
    //     // octopus_hg38_bowtie2_trimmomatic,
    //     // octopus_hg38_bwamem_fastp,
    //     // octopus_hg38_bwamem_cutadapt,
    //     // octopus_hg38_bwamem_trimmomatic
    //     gatk_hg38_bwamem_fastp
    // )

    // Run hap.py comparison for each HG38 variant (18 separate comparisons)
    // HAPPY_COMPARE(
    //     all_hg38_variants,
    //     truth_hg38_vcf,
    //     truth_hg38_vcf_idx,
    //     truth_hg38_bed,
    //     hg38_ref,
    //     hg38_ref_fai
    // )

    // ========================================
    // COMBINED COMPARISON (all QC tools mixed)
    // ========================================

    // Prepare combined channels (without QC distinction) - 6 comparisons
    // deepvariant_hg38_bwamem_combined = DEEPVARIANT_HG38_BWAMEM.out.vcf
    //     .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bwamem", "combined", "hg38", vcf, idx] }

    // gatk_hg38_bwamem_combined = GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bwamem", "combined", "hg38", vcf, idx] }

    // deepvariant_hg38_bowtie2_combined = DEEPVARIANT_HG38_BOWTIE2.out.vcf
    //     .join(DEEPVARIANT_HG38_BOWTIE2.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "deepvariant", "bowtie2", "combined", "hg38", vcf, idx] }

    // gatk_hg38_bowtie2_combined = GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.vcf
    //     .join(GATK_HAPLOTYPECALLER_HG38_BOWTIE2.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "gatk", "bowtie2", "combined", "hg38", vcf, idx] }

    // octopus_hg38_bowtie2_combined = TABIX_OCTOPUS_HG38_BOWTIE2.out.indexed_vcf
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bowtie2", "combined", "hg38", vcf, idx] }
    
    // octopus_hg38_bwamem_combined = TABIX_OCTOPUS_HG38_BWAMEM.out.indexed_vcf
    //     .map { meta, vcf, idx -> [meta.id, "octopus", "bwamem", "combined", "hg38", vcf, idx] }

    // Mix all for combined comparison
    // all_hg38_variants_combined = deepvariant_hg38_bwamem_combined.mix(
    //     deepvariant_hg38_bowtie2_combined,
    //     gatk_hg38_bowtie2_combined,
    //     // octopus_hg38_bowtie2_combined,
    //     // octopus_hg38_bwamem_combined
    //     gatk_hg38_bwamem_combined
    // )

    // Run hap.py for combined variants (6 comparisons: 3 callers × 2 aligners)
    // HAPPY_COMPARE_COMBINED(
    //     all_hg38_variants_combined,
    //     truth_hg38_vcf,
    //     truth_hg38_vcf_idx,
    //     truth_hg38_bed,
    //     hg38_ref,
    //     hg38_ref_fai
    // )

    // ========================================
    // BENCHMARKING WITH Truvari (Structural Variants - HG38) - COMMENTED OUT
    // ========================================

    // Load truth SV VCF file for HG37 (GIAB CMRG v1.00 structural variants) - COMMENTED OUT
    // def truth_hg37_sv_vcf = file(params.truth_hg37_sv_vcf)
    // def truth_hg37_sv_vcf_idx = file("${params.truth_hg37_sv_vcf}.tbi")
    // def truth_hg37_sv_bed = file(params.truth_hg37_sv_bed)

    // Load truth SV VCF file for HG38 (GIAB CMRG v1.00 structural variants) - COMMENTED OUT
    // def truth_hg38_sv_vcf = file(params.truth_hg38_sv_vcf)
    // def truth_hg38_sv_vcf_idx = file("${params.truth_hg38_sv_vcf}.tbi")
    // def truth_hg38_sv_bed = file(params.truth_hg38_sv_bed)

    // Prepare Manta HG38 variants (separate QC comparisons) - COMMENTED OUT
    // manta_hg38_bwamem_fastp = MANTA_HG38_BWAMEM.out.vcf
    //     .join(MANTA_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bwamem", "fastp", "hg38", vcf, idx] }

    // manta_hg38_bowtie2_fastp = MANTA_HG38_BOWTIE2.out.vcf
    //     .join(MANTA_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bowtie2", "fastp", "hg38", vcf, idx] }

    // Prepare Delly HG38 variants (separate QC comparisons) - COMMENTED OUT
    // delly_hg38_bwamem_fastp = DELLY_HG38_BWAMEM.out.vcf
    //     .join(DELLY_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bwamem", "fastp", "hg38", vcf, idx] }

    // delly_hg38_bowtie2_fastp = DELLY_HG38_BOWTIE2.out.vcf
    //     .join(DELLY_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bowtie2", "fastp", "hg38", vcf, idx] }

    // Prepare combined SV variants for Truvari benchmarking (1 comparison: 1 aligner) - COMMENTED OUT
    // Using combined DELLY+MANTA callsets
    // combined_sv_for_truvari = COMBINE_SV.out.combined_vcf
    //     .map { sample_id, aligner, qc, reference, vcf, idx ->
    //         [sample_id, "combined_delly_manta", aligner, qc, reference, vcf, idx]
    //     }

    // Combine individual callers + combined callset for Truvari benchmarking - COMMENTED OUT
    // Total: 3 comparisons (2 individual + 1 combined)
    // all_hg38_sv_variants = manta_hg38_bwamem_fastp.mix(
    //     // manta_hg38_bowtie2_fastp,
    //     delly_hg38_bwamem_fastp,
    //     // delly_hg38_bowtie2_fastp,
    //     combined_sv_for_truvari
    // )

    // Run Truvari comparison on individual AND combined SV callsets - COMMENTED OUT
    // 3 comparisons: MANTA (1), DELLY (1), Combined DELLY+MANTA (1)
    // TRUVARI_COMPARE(
    //     all_hg38_sv_variants,
    //     truth_hg38_sv_vcf,
    //     truth_hg38_sv_vcf_idx,
    //     truth_hg38_sv_bed,
    //     hg38_ref,
    //     hg38_ref_fai
    // )

    // ========================================
    // BENCHMARKING WITH Truvari (Structural Variants - HG38) - COMMENTED OUT
    // ========================================

    // Load truth SV VCF file for HG38 (GIAB CMRG v1.00 structural variants) - COMMENTED OUT
    // def truth_hg38_sv_vcf = file(params.truth_hg38_sv_vcf)
    // def truth_hg38_sv_vcf_idx = file("${params.truth_hg38_sv_vcf}.tbi")
    // def truth_hg38_sv_bed = file(params.truth_hg38_sv_bed)

    // Prepare Manta HG38 variants (separate QC comparisons) - COMMENTED OUT
    // manta_hg38_bwamem_fastp = MANTA_HG38_BWAMEM.out.vcf
    //     .join(MANTA_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bwamem", "fastp", "hg38", vcf, idx] }

    // manta_hg38_bowtie2_fastp = MANTA_HG38_BOWTIE2.out.vcf
    //     .join(MANTA_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bowtie2", "fastp", "hg38", vcf, idx] }

    // Prepare Delly HG38 variants (separate QC comparisons)
    // delly_hg38_bwamem_fastp = DELLY_HG38_BWAMEM.out.vcf
    //     .join(DELLY_HG38_BWAMEM.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bwamem", "fastp", "hg38", vcf, idx] }

    // delly_hg38_bowtie2_fastp = DELLY_HG38_BOWTIE2.out.vcf
    //     .join(DELLY_HG38_BOWTIE2.out.tbi)
    //     .filter { meta, vcf, idx -> meta.qc_tool == 'fastp' }
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bowtie2", "fastp", "hg38", vcf, idx] }

    // Combine all HG38 structural variants for separate comparisons
    // all_hg38_sv_variants = manta_hg38_bwamem_fastp.mix(
    //     manta_hg38_bowtie2_fastp,
    //     delly_hg38_bwamem_fastp,
    //     delly_hg38_bowtie2_fastp
    // )

    // Run Truvari comparison for each HG38 structural variant (4 comparisons: 2 callers × 2 aligners)
    // TRUVARI_COMPARE(
    //     all_hg38_sv_variants,
    //     truth_hg38_sv_vcf,
    //     truth_hg38_sv_vcf_idx,
    //     truth_hg38_sv_bed,
    //     hg38_ref,
    //     hg38_ref_fai
    // )

    // ========================================
    // COMBINED COMPARISON FOR SV (all QC tools mixed) - HG38 - COMMENTED OUT
    // ========================================

    // Prepare combined channels for structural variants - HG37 - COMMENTED OUT
    // manta_hg37_bwamem_combined = MANTA_HG37_BWAMEM.out.vcf
    //     .join(MANTA_HG37_BWAMEM.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bwamem", "combined", "hg37", vcf, idx] }

    // manta_hg37_bowtie2_combined = MANTA_HG37_BOWTIE2.out.vcf
    //     .join(MANTA_HG37_BOWTIE2.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bowtie2", "combined", "hg37", vcf, idx] }

    // delly_hg37_bwamem_combined = DELLY_HG37_BWAMEM.out.vcf
    //     .join(DELLY_HG37_BWAMEM.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bwamem", "combined", "hg37", vcf, idx] }

    // delly_hg37_bowtie2_combined = DELLY_HG37_BOWTIE2.out.vcf
    //     .join(DELLY_HG37_BOWTIE2.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bowtie2", "combined", "hg37", vcf, idx] }

    // Mix all for combined SV comparison - HG37
    // all_hg37_sv_variants_combined = manta_hg37_bwamem_combined.mix(
    //     manta_hg37_bowtie2_combined,
    //     delly_hg37_bwamem_combined,
    //     delly_hg37_bowtie2_combined
    // )

    // Run Truvari for combined SV variants (4 comparisons: 2 SV callers × 2 aligners) - HG37
    // TRUVARI_COMPARE_COMBINED(
    //     all_hg37_sv_variants_combined,
    //     truth_hg37_sv_vcf,
    //     truth_hg37_sv_vcf_idx,
    //     truth_hg37_sv_bed,
    //     hg37_ref,
    //     hg37_ref_fai
    // )

    // ========================================
    // COMBINED COMPARISON FOR SV (all QC tools mixed) - HG38 - COMMENTED OUT
    // ========================================

    // Prepare combined channels for structural variants - HG38 - COMMENTED OUT
    // manta_hg38_bwamem_combined = MANTA_HG38_BWAMEM.out.vcf
    //     .join(MANTA_HG38_BWAMEM.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bwamem", "combined", "hg38", vcf, idx] }

    // manta_hg38_bowtie2_combined = MANTA_HG38_BOWTIE2.out.vcf
    //     .join(MANTA_HG38_BOWTIE2.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "manta", "bowtie2", "combined", "hg38", vcf, idx] }

    // delly_hg38_bwamem_combined = DELLY_HG38_BWAMEM.out.vcf
    //     .join(DELLY_HG38_BWAMEM.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bwamem", "combined", "hg38", vcf, idx] }

    // delly_hg38_bowtie2_combined = DELLY_HG38_BOWTIE2.out.vcf
    //     .join(DELLY_HG38_BOWTIE2.out.tbi)
    //     .map { meta, vcf, idx -> [meta.id, "delly", "bowtie2", "combined", "hg38", vcf, idx] }

    // Mix all for combined SV comparison - COMMENTED OUT
    // all_hg38_sv_variants_combined = manta_hg38_bwamem_combined.mix(
    //     // manta_hg38_bowtie2_combined,
    //     delly_hg38_bwamem_combined
    //     // delly_hg38_bowtie2_combined
    // )

    // Run Truvari for combined SV variants (2 comparisons: 2 SV callers × 1 aligner) - COMMENTED OUT
    // TRUVARI_COMPARE_COMBINED(
    //     all_hg38_sv_variants_combined,
    //     truth_hg38_sv_vcf,
    //     truth_hg38_sv_vcf_idx,
    //     truth_hg38_sv_bed,
    //     hg38_ref,
    //     hg38_ref_fai
    // )
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

    COMBINED VARIANT CALLSETS (BENCHMARKED):
      SNP/INDEL (GATK + DeepVariant intersection):
        - Combined callsets: ${params.outdir}/04_variant_calling/combined/snp_indel/
        - Statistics: Check *_stats.txt files for comparison metrics
      
      SV (DELLY + MANTA union):
        - Combined callsets: ${params.outdir}/04_variant_calling/combined/sv/
        - Statistics: Check *_stats.txt files for comparison metrics

    BENCHMARKING RESULTS (individual + combined callsets):
      hap.py comparisons (SNP/INDEL vs GIAB truth):
        - Individual callers: DeepVariant (2 aligners), GATK (2 aligners)
        - Combined callset: GATK+DeepVariant (2 aligners)
        - Total: 6 comparisons
        - Results: ${params.outdir}/05_happy_comparison/
      
      Truvari comparisons (SV vs GIAB truth):
        - Individual callers: MANTA (2 aligners), DELLY (2 aligners)
        - Combined callset: DELLY+MANTA (2 aligners)
        - Total: 6 comparisons
        - Results: ${params.outdir}/05_truvari_comparison/
    ================================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution failed: ${workflow.errorMessage}"
}
