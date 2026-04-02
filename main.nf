#!/usr/bin/env nextflow

/*
========================================================================================
    QTL2_NF: Nextflow pipeline for multiparental mouse QTL analysis
========================================================================================
    Github : https://github.com/Dpbudke/QTL2_NF
    Author : Dpbudke
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Print help message
def helpMessage() {
    log.info"""
    ========================================
     QTL2_NF Pipeline v${workflow.manifest.version}
    ========================================
    
    Usage:
        nextflow run main.nf --phenotype_file <file> --study_prefix <prefix> [options]
    
    Required Arguments:
        --phenotype_file     Path to master phenotype CSV (from your pre-pipeline RMD)
        --study_prefix       Study identifier prefix for output files

    Optional Arguments:
        --analysis_type     Analysis preset that auto-configures outdir, sample_filter, exclude_covariates,
                            and interactive_covar. Valid values:
                              - Diet_additive    : Full cohort, additive model (default settings)
                              - Diet_interactive : Full cohort, GxE interaction with Diet covariate
                              - AIN76_only       : AIN76 diet subset, excludes Diet covariate
                              - HC_only          : HC diet subset, excludes Diet covariate
        --finalreport_files Path to GeneSeek FinalReport file(s) - accepts glob patterns (e.g., 'Data/FinalReport*.txt')
        --outdir            Output directory [default: results, or Results_final/{analysis_type} if --analysis_type set]
        --study_type        Optional prefix for output directory (e.g. 'LCMS' or 'eQTL').
                            When set, outdir becomes Results_final/{study_type}/{analysis_type}
        --auto_prefix_samples    Automatically prefix sample IDs [default: false]
        --test_mode         Positive control: chr 2 only, coat_color as phenotype [default: false]
        --lod_threshold     LOD threshold for filtering QTLs before permutation testing [default: 7.0]
        --sample_filter     JSON filter for sample subsetting by covariates (e.g., '{"Sex":["male"],"Diet":["hc"]}') [default: null]
        --exclude_covariates  Comma-separated covariates to exclude from model (e.g., 'Diet,batch') [default: null]
        --interactive_covar Covariate to test as interactive term for GxE analysis (e.g., 'Diet') [default: null]
    
    Pipeline Modules:
        Module 1: Phenotype processing and validation
        Module 2: Genotype processing (requires --finalreport_files)
        Module 3: Control file generation
        Module 4: Cross2 object creation
        Module 5: Genome scan preparation (with increased resources)
        Module 6: HPC array-based genome scanning
        Module 7: Chunked permutation testing (DO_Pipe approach - 1000 permutations)
        Module 8: Significant QTL identification
        Module 9: QTL visualizations (generates plot_coefCC plots for 99% QTLs)
    
    Examples:
        # Using analysis_type preset (recommended):
        nextflow run main.nf --analysis_type Diet_additive --phenotype_file Data/pheno.csv --study_prefix DOChln
        nextflow run main.nf --analysis_type HC_only --phenotype_file Data/pheno.csv --study_prefix DOChln

        # Manual configuration (backwards compatible):
        nextflow run main.nf \\
            --phenotype_file Data/formatted_phenotypes_for_nextflow.csv \\
            --finalreport_files 'Data/FinalReport*.txt' \\
            --study_prefix DOChln \\
            --outdir results
    
    Profiles:
        -profile standard   Local execution with Docker (default)
        -profile test       Run with test data (chr 2 only)
        -profile docker     Explicit Docker profile
    
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

// Check required parameters
if (!params.phenotype_file) {
    log.error "Error: --phenotype_file is required"
    helpMessage()
    exit 1
}

if (!params.study_prefix) {
    log.error "Error: --study_prefix is required" 
    helpMessage()
    exit 1
}

// Validate input files exist
phenotype_input = file(params.phenotype_file)
if (!phenotype_input.exists()) {
    log.error "Error: Phenotype file not found: ${params.phenotype_file}"
    exit 1
}

/*
========================================================================================
    ANALYSIS TYPE PRESET CONFIGURATION
========================================================================================
*/

// Valid analysis types
def validAnalysisTypes = ['Diet_additive', 'Diet_interactive', 'AIN76_only', 'HC_only']

// Validate analysis_type if provided
if (params.analysis_type && !validAnalysisTypes.contains(params.analysis_type)) {
    log.error """
    Invalid --analysis_type: '${params.analysis_type}'
    Valid options: ${validAnalysisTypes.join(', ')}
    """.stripIndent()
    exit 1
}

// Compute effective parameters based on analysis_type
// outdir logic: use explicit --outdir if provided, else auto-compute from analysis_type (+study_type), else default to 'results'
def effective_outdir = params.outdir ?:
    (params.analysis_type ?
        (params.study_type ? "Results_final/${params.study_type}/${params.analysis_type}" : "Results_final/${params.analysis_type}")
        : 'results')

def effective_sample_filter = params.sample_filter ?: (
    params.analysis_type == 'AIN76_only' ? '{"Diet": ["ain76a"]}' :
    params.analysis_type == 'HC_only' ? '{"Diet": ["hc"]}' : null
)

def effective_exclude_covariates = params.exclude_covariates ?: (
    params.analysis_type in ['AIN76_only', 'HC_only'] ? 'Diet' : null
)

def effective_interactive_covar = params.interactive_covar ?: (
    params.analysis_type == 'Diet_interactive' ? 'Diet' : null
)

// Override params.outdir so modules use the computed value
params.outdir = effective_outdir

// Log analysis type configuration
if (params.analysis_type) {
    def analysisDescriptions = [
        'Diet_additive': 'full cohort, additive model',
        'Diet_interactive': 'full cohort, GxE interaction with Diet',
        'AIN76_only': 'AIN76 diet subset',
        'HC_only': 'HC diet subset'
    ]
    log.info "Analysis type: ${params.analysis_type} (${analysisDescriptions[params.analysis_type]})"
    log.info """
    Auto-configured settings for ${params.analysis_type}:
        outdir             : ${effective_outdir}
        sample_filter      : ${effective_sample_filter ?: 'null (full cohort)'}
        exclude_covariates : ${effective_exclude_covariates ?: 'null'}
        interactive_covar  : ${effective_interactive_covar ?: 'null (additive model)'}
    """.stripIndent()
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { PHENOTYPE_PROCESS } from './modules/01_phenotype_process.nf'
include { GENOTYPE_PROCESS  } from './modules/02_genotype_process.nf'
include { GENERATE_CONTROL_FILE } from './modules/03_control_file_generation.nf'
include { CREATE_CROSS2_OBJECT } from './modules/04_cross2_creation.nf'
include { PREPARE_GENOME_SCAN_SETUP } from './modules/05_prepare_genome_scan.nf'
include { GENOME_SCAN_SETUP; GENOME_SCAN_BATCH; COMBINE_BATCH_RESULTS } from './modules/06_qtl_analysis.nf'
include { FILTER_PEAKS_BY_REGION } from './modules/06b_filter_peaks_by_region.nf'
include { CHUNKED_PERMUTATION_TESTING } from './modules/07_permutation_testing.nf'
include { IDENTIFY_SIGNIFICANT_QTLS } from './modules/08_identify_significant_qtls.nf'
include { VISUALIZE_QTLS } from './modules/09_visualize.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    log.info """
    ========================================
     QTL2_NF Pipeline Started
    ========================================
    phenotype_file : ${params.phenotype_file}
    study_prefix   : ${params.study_prefix}
    outdir         : ${effective_outdir}
    test_mode      : ${params.test_mode ? 'true (chr 2 only)' : 'false'}
    lod_threshold  : ${params.lod_threshold}
    ========================================
    """.stripIndent()

    // Create input channels
    ch_phenotype_file = Channel.fromPath(params.phenotype_file, checkIfExists: true)
    ch_study_prefix   = Channel.value(params.study_prefix)
    ch_gtf_file       = Channel.fromPath(params.gtf_file, checkIfExists: true)

    // MODULE 1: Phenotype Processing
    PHENOTYPE_PROCESS(
        ch_phenotype_file,
        ch_study_prefix,
        Channel.value(effective_sample_filter ?: "null"),
        Channel.value(effective_exclude_covariates ?: "null")
    )
    
    // MODULE 2: Genotype Processing (if FinalReport files provided)
    if (params.finalreport_files) {
        ch_finalreport = Channel.fromPath(params.finalreport_files, checkIfExists: true)
                                .collect()  // Collect all files into a single emission
        ch_allele_codes = Channel.fromPath("${projectDir}/Data/GM_allelecodes.csv", checkIfExists: true)
        ch_founder_genos = Channel.fromPath("${projectDir}/Data/GM_foundergeno*.csv", checkIfExists: true).collect()
        ch_gmap_files = Channel.fromPath("${projectDir}/Data/GM_gmap*.csv", checkIfExists: true).collect()
        ch_pmap_files = Channel.fromPath("${projectDir}/Data/GM_pmap*.csv", checkIfExists: true).collect()

        GENOTYPE_PROCESS(
            ch_finalreport,
            PHENOTYPE_PROCESS.out.valid_samples,
            ch_allele_codes,
            ch_study_prefix
        )

        // MODULE 3: Control File and Cross2 Object Creation
        GENERATE_CONTROL_FILE(
            PHENOTYPE_PROCESS.out.pheno,
            PHENOTYPE_PROCESS.out.covar,
            GENOTYPE_PROCESS.out.geno_files.collect(),
            GENOTYPE_PROCESS.out.allele_codes,
            ch_founder_genos,
            ch_gmap_files,
            ch_pmap_files,
            ch_study_prefix
        )
        
        CREATE_CROSS2_OBJECT(
            GENERATE_CONTROL_FILE.out.control_file,
            ch_founder_genos.mix(
                ch_gmap_files,
                ch_pmap_files,
                PHENOTYPE_PROCESS.out.pheno,
                PHENOTYPE_PROCESS.out.covar,
                GENOTYPE_PROCESS.out.geno_files
            ).collect(),
            ch_study_prefix
        )
        
        // Display results for Module 2
        GENOTYPE_PROCESS.out.geno_files.view { "Genotype files created: $it" }
        GENOTYPE_PROCESS.out.validation_report.view { "Genotype validation report: $it" }
        GENOTYPE_PROCESS.out.allele_codes.view { "Reference allele codes downloaded: $it" }
        
        // Display results for Module 3
        GENERATE_CONTROL_FILE.out.control_file.view { "Control file created: $it" }
        CREATE_CROSS2_OBJECT.out.cross2_object.view { "Cross2 object created: $it" }
        CREATE_CROSS2_OBJECT.out.validation_report.view { "Cross2 validation report: $it" }
        
        // MODULE 5: Genome Scan Preparation (simplified single-process)
        PREPARE_GENOME_SCAN_SETUP(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            ch_study_prefix
        )

        // MODULE 6: Full Genome Scan (Batch Processing)
        GENOME_SCAN_SETUP(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            ch_study_prefix,
            Channel.value(params.lod_threshold)
        )

        GENOME_SCAN_SETUP.out.summary.view { "Chunking summary: $it" }

        // Use batch_ids from GENOME_SCAN_SETUP stdout
        GENOME_SCAN_SETUP.out.batch_ids
            .splitText()
            .map { it.trim() }
            .filter { it != "" }
            .set { ch_batch_ids }

        // Process batches in parallel
        // Use .first() to convert queue channels to value channels for reuse with each
        GENOME_SCAN_BATCH(
            ch_batch_ids,
            CREATE_CROSS2_OBJECT.out.cross2_object.first(),
            PREPARE_GENOME_SCAN_SETUP.out.genoprob.first(),
            PREPARE_GENOME_SCAN_SETUP.out.alleleprob.first(),
            PREPARE_GENOME_SCAN_SETUP.out.kinship_loco.first(),
            GENOME_SCAN_SETUP.out.chunk_file.first(),
            GENOME_SCAN_SETUP.out.batch_file.first(),
            ch_study_prefix,
            Channel.value(params.lod_threshold),
            Channel.value(effective_interactive_covar ?: "null")
        )

        // Combine all batch results
        COMBINE_BATCH_RESULTS(
            GENOME_SCAN_BATCH.out.batch_results.collect(),
            GENOME_SCAN_BATCH.out.batch_peaks.collect(),
            GENOME_SCAN_BATCH.out.batch_filtered_phenos.collect(),
            GENOME_SCAN_BATCH.out.batch_log.collect(),
            ch_study_prefix,
            Channel.value(params.lod_threshold)
        )

        COMBINE_BATCH_RESULTS.out.batch_summary.view { "Batch processing summary: $it" }
        COMBINE_BATCH_RESULTS.out.peaks.view { "Peaks found: $it" }

        // MODULE 6b: Optional Regional QTL Filtering (before permutation testing)
        def filtered_phenos_ch
        def regional_filter_file = "${effective_outdir}/06b_regional_filtering/${params.study_prefix}_regional_filtered_phenotypes.txt"

        if (params.qtl_region != null) {
            log.info "Filtering QTL peaks to genomic region: ${params.qtl_region}"

            FILTER_PEAKS_BY_REGION(
                COMBINE_BATCH_RESULTS.out.peaks,
                PREPARE_GENOME_SCAN_SETUP.out.genetic_map,
                Channel.fromPath(params.gtf_file),
                ch_study_prefix,
                Channel.value(params.qtl_region)
            )

            FILTER_PEAKS_BY_REGION.out.filter_report.view { "Regional filtering report: $it" }
            filtered_phenos_ch = FILTER_PEAKS_BY_REGION.out.filtered_phenotypes
        } else if (file(regional_filter_file).exists()) {
            log.info "Using existing regional filtered phenotypes from previous run: ${regional_filter_file}"
            filtered_phenos_ch = Channel.fromPath(regional_filter_file)
        } else {
            log.info "No regional filtering applied - using all LOD-filtered phenotypes"
            filtered_phenos_ch = COMBINE_BATCH_RESULTS.out.filtered_phenotypes
        }

        // MODULE 7: Chunked Permutation Testing (DO_Pipe Approach with Smart Pre-filtering)
        if (params.run_perm_benchmark) {
            log.info "Running permutation testing WITH performance benchmarking (adds 30-40 min)"
        } else {
            log.info "Running permutation testing with fast setup (skipping benchmark)"
        }

        CHUNKED_PERMUTATION_TESTING(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            PREPARE_GENOME_SCAN_SETUP.out.genoprob,
            PREPARE_GENOME_SCAN_SETUP.out.kinship_loco,
            filtered_phenos_ch,
            ch_study_prefix,
            Channel.value(params.lod_threshold),
            Channel.value(params.run_perm_benchmark),
            Channel.value(params.perm_per_chunk),
            Channel.value(params.chunks_per_batch)
        )

        // MODULE 8: Significant QTL Identification
        IDENTIFY_SIGNIFICANT_QTLS(
            CHUNKED_PERMUTATION_TESTING.out.filtered_cross2,
            COMBINE_BATCH_RESULTS.out.scan_results,
            CHUNKED_PERMUTATION_TESTING.out.permutation_matrix,
            CHUNKED_PERMUTATION_TESTING.out.permutation_thresholds,
            ch_study_prefix
        )

        // MODULE 9: QTL Visualizations
        // NOTE: This module generates individual PNG plots for significant QTLs
        log.info "Running QTL visualization (plot_coefCC for 99% significant QTLs)"

        VISUALIZE_QTLS(
            PREPARE_GENOME_SCAN_SETUP.out.alleleprob,
            COMBINE_BATCH_RESULTS.out.scan_results,
            CHUNKED_PERMUTATION_TESTING.out.filtered_cross2,
            IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls,
            ch_study_prefix
        )

        // Display results for Module 9
        VISUALIZE_QTLS.out.validation_report.view { "QTL visualization report: $it" }

        // Display results for Modules 5-8
        PREPARE_GENOME_SCAN_SETUP.out.setup_report.view { "Genome scan prep report: $it" }
        CHUNKED_PERMUTATION_TESTING.out.setup_log.view { "Permutation setup log: $it" }
        CHUNKED_PERMUTATION_TESTING.out.aggregation_log.view { "Permutation aggregation log: $it" }
        CHUNKED_PERMUTATION_TESTING.out.permutation_thresholds.view { "Significance thresholds: $it" }
        IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls.view { "Significant QTLs identified: $it" }
        IDENTIFY_SIGNIFICANT_QTLS.out.qtl_summary.view { "QTL summary report: $it" }
        
    } else {
        log.info "Skipping genotype processing - no FinalReport files specified"
        log.info "To include genotype processing, add: --finalreport_files 'Data/FinalReport*.txt'"
    }
    
    // Display results for Module 1
    PHENOTYPE_PROCESS.out.covar.view { "Covariate file created: $it" }
    PHENOTYPE_PROCESS.out.pheno.view { "Phenotype file created: $it" }
    PHENOTYPE_PROCESS.out.validation_report.view { "Validation report: $it" }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================
     QTL2_NF Pipeline Completed
    ========================================
    Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Time     : ${workflow.duration}
    Workdir  : ${workflow.workDir}
    Results  : ${effective_outdir}
    Test Mode: ${params.test_mode ? 'Enabled (chr 2 only)' : 'Disabled'}
    ========================================
    """.stripIndent()
    
    // Clean up work directory if successful (optional)
    if (workflow.success && params.cleanup) {
        log.info "Cleaning up work directory: ${workflow.workDir}"
        workflow.workDir.deleteDir()
    }
}

workflow.onError {
    log.error """
    ========================================
     Pipeline execution stopped with error
    ========================================
    Error message: ${workflow.errorMessage}
    Error report : ${workflow.errorReport}
    ========================================
    """.stripIndent()
}