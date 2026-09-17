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
include { FASTQC_FASTP } from './modules/fastqc_fastp'
include { BWAMEM_HG38 } from './modules/bwa'
include { COUNT_READS } from './modules/count_reads'
include { MANTA_HG38_BWAMEM } from './modules/manta'
include { DELLY_HG38_BWAMEM } from './modules/delly'
include { DEEPVARIANT_HG38_BWAMEM } from './modules/deepvariant'
include { DV_PASS_FILTER_HG38_BWAMEM } from './modules/deepvariant'
include { GATK_HAPLOTYPECALLER_HG38_BWAMEM } from './modules/gatk'
include { GATK_HARDFILTER_HG38_BWAMEM } from './modules/gatk'
params.input_dir = "${projectDir}/data"
params.outdir = "${projectDir}/results"
params.pattern = "*_R{1,2}_001.fastq.gz"

// Reference genome parameters
params.hg38_index = "${projectDir}/reference/hg38"

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
    // Files live either in `${input_dir}/<run>/<sample>_R{1,2}_001.fastq.gz` or directly in
    // `${input_dir}/<sample>_R{1,2}_001.fastq.gz`. The grouping key embeds the run so
    // identically-named samples in different runs never get paired together, and meta.run
    // is used to mirror every result per run folder. For files directly in input_dir the run
    // is the reference/vcf_hdd/<run>/ folder holding that sample's laboratory VCF
    // (<sample>.vcf, i.e. meta.id without the lane suffix); fallback: the input_dir name.
    def input_dir = file(params.input_dir).toAbsolutePath().normalize()
    def lab_vcfs_by_sample = files("${projectDir}/reference/vcf_hdd/*/*.vcf").groupBy { it.baseName }

    fastq_ch = Channel
        .fromFilePairs(["${input_dir}/${params.pattern}", "${input_dir}/*/${params.pattern}"]) { file ->
            def sample = file.name.replaceAll(/_R[12]_001\.fastq\.gz$/, '')
            def run = file.parent.name
            if (file.parent.toAbsolutePath().normalize() == input_dir) {
                def lab_runs = lab_vcfs_by_sample.get(sample.replaceAll(/_L0*\d+$/, ''), [])*.parent*.name
                run = lab_runs.size() == 1 ? lab_runs[0] : input_dir.name
            }
            "${run}__${sample}"
        }
        .ifEmpty { error "No files match pattern `${params.pattern}` in ${input_dir} or its run sub-folders" }
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

    // Step 3: Run FastQC on processed data
    FASTQC_FASTP(FASTP.out.reads)

    // Prepare metadata for alignments with each genome
    
    
    fastp_hg38 = FASTP.out.reads
        .map { tuple -> 
            def meta = tuple[0].clone()
            meta.qc_tool = "fastp"
            [meta, tuple[1]]
        }

    // Combine all inputs for each genome and run alignments
    hg38_inputs = fastp_hg38 // .mix(cutadapt_hg38) // .mix(trimmomatic_hg38)


    // Run BWA-MEM alignments
    BWAMEM_HG38(hg38_inputs)

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

    // ========================================
    // VARIANT CALLING WITH MANTA
    // ========================================

    exome_bed_ch = params.exome_mode && params.exome_bed_hg38 ? Channel.value(file(params.exome_bed_hg38)) : Channel.value(file('NO_FILE'))
    exome_bed_tbi_ch = params.exome_mode && params.exome_bed_hg38 ? Channel.value(file(params.exome_bed_hg38 + '.tbi')) : Channel.value(file('NO_FILE'))
    MANTA_HG38_BWAMEM(bwamem_hg38_bams_manta, exome_bed_ch, exome_bed_tbi_ch)

    // ========================================
    // VARIANT CALLING WITH DELLY
    // ========================================
    DELLY_HG38_BWAMEM(bwamem_hg38_bams_delly)

    // ========================================
    // VARIANT CALLING WITH DEEPVARIANT
    // ========================================

    DEEPVARIANT_HG38_BWAMEM(bwamem_hg38_bams_deepvariant, exome_bed_ch, exome_bed_tbi_ch)

    // Keep only PASS DeepVariant calls (drop RefCall, LowQual, NoCall)
    DV_PASS_FILTER_HG38_BWAMEM(
        DEEPVARIANT_HG38_BWAMEM.out.vcf
            .join(DEEPVARIANT_HG38_BWAMEM.out.tbi)
    )

    // ========================================
    // VARIANT CALLING WITH GATK HAPLOTYPECALLER
    // ========================================

    GATK_HAPLOTYPECALLER_HG38_BWAMEM(bwamem_hg38_bams_gatk, exome_bed_ch, exome_bed_tbi_ch)

    GATK_HARDFILTER_HG38_BWAMEM(
        GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.vcf
            .join(GATK_HAPLOTYPECALLER_HG38_BWAMEM.out.tbi)
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
