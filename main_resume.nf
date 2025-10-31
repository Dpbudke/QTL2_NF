#!/usr/bin/env nextflow

/*
========================================================================================
    QTL2_NF: Enhanced pipeline with resume-from-step capability
========================================================================================
    Resume pipeline execution from any step using existing outputs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Default parameters
params.phenotype_file = null
params.finalreport_files = null
params.study_prefix = null
params.outdir = 'results'
params.auto_prefix_samples = false
params.test_mode = false
params.lod_threshold = 7.0
params.sample_filter = null
params.help = false

// Resume parameters
params.resume_from = null  // Options: phenotype, genotype, control, cross2, prepare_scan, genome_scan, permutation, perm_aggregate, significant_qtls, qtlviewer
params.input_dir = null    // Directory to load cached results from (defaults to outdir if not specified)

// QTL Viewer parameters (Module 9 - optional, supports HPC and local deployment)
params.run_qtlviewer = false  // Set to true to run QTL Viewer setup (supports HPC and local deployment)

// Print help message
def helpMessage() {
    log.info"""
    ========================================
     QTL2_NF Pipeline with Resume Feature
    ========================================

    Usage:
        nextflow run main_resume.nf --phenotype_file <file> --study_prefix <prefix> [options]

    Resume from specific step:
        nextflow run main_resume.nf --resume_from <step> --study_prefix <prefix> [options]

    Required Arguments (for full run):
        --phenotype_file     Path to master phenotype CSV
        --study_prefix       Study identifier prefix for output files

    Resume Options:
        --resume_from        Step to resume from:
                            phenotype       - Start from phenotype processing
                            genotype        - Start from genotype processing
                            control         - Start from control file generation
                            cross2          - Start from cross2 object creation
                            prepare_scan    - Start from genome scan preparation
                            genome_scan     - Start from genome scanning
                            regional_filter - Start from regional QTL filtering (requires --qtl_region)
                            permutation     - Start from permutation testing
                            significant_qtls - Start from QTL identification
                            qtlviewer       - Start from QTL Viewer setup

    Other Arguments:
        --finalreport_files  Path to GeneSeek FinalReport file(s)
        --outdir            Output directory for new results [default: results]
        --input_dir         Directory to load cached results from (defaults to outdir if not specified)
        --lod_threshold     LOD threshold for filtering [default: 7.0]
        --sample_filter     JSON filter for sample subsetting
        --run_qtlviewer     Enable QTL Viewer setup (Module 9 - supports HPC and local deployment) [default: false]

    Examples:
        # Full pipeline
        nextflow run main_resume.nf --phenotype_file Data/pheno.csv --study_prefix MyStudy

        # Resume from genome scan (skips all earlier steps)
        nextflow run main_resume.nf --resume_from genome_scan --study_prefix MyStudy

        # Resume from permutation testing
        nextflow run main_resume.nf --resume_from permutation --study_prefix MyStudy --lod_threshold 5.0

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

// Validate resume_from parameter
def valid_resume_steps = ['phenotype', 'genotype', 'control', 'cross2', 'prepare_scan',
                         'genome_scan', 'regional_filter', 'permutation', 'perm_aggregate', 'significant_qtls', 'qtlviewer']

if (params.resume_from && !(params.resume_from in valid_resume_steps)) {
    log.error "Invalid --resume_from value: ${params.resume_from}"
    log.error "Valid options: ${valid_resume_steps.join(', ')}"
    exit 1
}

// Check required parameters for full run
if (!params.resume_from) {
    if (!params.phenotype_file) {
        log.error "Error: --phenotype_file is required for full pipeline run"
        helpMessage()
        exit 1
    }
    if (!params.study_prefix) {
        log.error "Error: --study_prefix is required"
        helpMessage()
        exit 1
    }
}

if (params.resume_from && !params.study_prefix) {
    log.error "Error: --study_prefix is required when using --resume_from"
    exit 1
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
include { CHUNKED_PERMUTATION_TESTING; PERM_AGGREGATE } from './modules/07_permutation_testing.nf'
include { IDENTIFY_SIGNIFICANT_QTLS } from './modules/08_identify_significant_qtls.nf'
include { PREPARE_QTLVIEWER_DATA } from './modules/09_qtl_viewer.nf'

/*
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

def checkFileExists(filepath, description) {
    def file = file(filepath)
    if (!file.exists()) {
        log.error "Required file not found for resume: ${description}"
        log.error "Path: ${filepath}"
        log.error "Make sure previous pipeline steps completed successfully."
        exit 1
    }
    return file
}

def shouldRunStep(currentStep, resumeFrom) {
    if (!resumeFrom) return true

    def stepOrder = ['phenotype', 'genotype', 'control', 'cross2', 'prepare_scan',
                    'genome_scan', 'regional_filter', 'permutation', 'perm_aggregate', 'significant_qtls', 'qtlviewer']

    def currentIndex = stepOrder.indexOf(currentStep)
    def resumeIndex = stepOrder.indexOf(resumeFrom)

    return currentIndex >= resumeIndex
}

// Function to check if we should load files (step was skipped) or use process outputs (step ran)
def shouldLoadFromFiles(currentStep, resumeFrom) {
    if (!resumeFrom) return false

    def stepOrder = ['phenotype', 'genotype', 'control', 'cross2', 'prepare_scan',
                    'genome_scan', 'permutation', 'significant_qtls', 'qtlviewer']

    def currentIndex = stepOrder.indexOf(currentStep)
    def resumeIndex = stepOrder.indexOf(resumeFrom)

    return currentIndex < resumeIndex
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Set input_dir to outdir if not specified
    def input_dir = params.input_dir ?: params.outdir

    log.info """
    ========================================
     QTL2_NF Pipeline Started
    ========================================
    study_prefix   : ${params.study_prefix}
    input_dir      : ${input_dir}
    outdir         : ${params.outdir}
    resume_from    : ${params.resume_from ?: 'beginning (full run)'}
    lod_threshold  : ${params.lod_threshold}
    ========================================
    """.stripIndent()

    // Initialize channels
    ch_study_prefix = Channel.value(params.study_prefix)

    // MODULE 1: Phenotype Processing
    if (shouldRunStep('phenotype', params.resume_from)) {
        log.info "Running PHENOTYPE_PROCESS"
        ch_phenotype_file = Channel.fromPath(params.phenotype_file, checkIfExists: true)

        PHENOTYPE_PROCESS(
            ch_phenotype_file,
            ch_study_prefix,
            Channel.value(params.sample_filter ?: "null")
        )

        ch_pheno = PHENOTYPE_PROCESS.out.pheno
        ch_covar = PHENOTYPE_PROCESS.out.covar
        ch_valid_samples = PHENOTYPE_PROCESS.out.valid_samples

        PHENOTYPE_PROCESS.out.validation_report.view { "Validation report: $it" }
    } else {
        log.info "Skipping PHENOTYPE_PROCESS - loading from existing files"
        ch_pheno = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_pheno.csv", "phenotype file"))
        ch_covar = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_covar.csv", "covariate file"))
        ch_valid_samples = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_valid_samples.txt", "valid samples file"))
    }

    // MODULE 2: Genotype Processing
    if (shouldRunStep('genotype', params.resume_from) && params.finalreport_files) {
        log.info "Running GENOTYPE_PROCESS"
        ch_finalreport = Channel.fromPath(params.finalreport_files, checkIfExists: true).collect()

        GENOTYPE_PROCESS(
            ch_finalreport,
            ch_valid_samples,
            ch_study_prefix
        )

        ch_geno_files = GENOTYPE_PROCESS.out.geno_files
        ch_allele_codes = GENOTYPE_PROCESS.out.allele_codes

        GENOTYPE_PROCESS.out.validation_report.view { "Genotype validation report: $it" }
    } else if (shouldRunStep('genotype', params.resume_from)) {
        log.info "Skipping GENOTYPE_PROCESS - loading from existing files"
        ch_geno_files = Channel.fromPath("${input_dir}/02_genotype_processing/${params.study_prefix}_geno*.csv").collect()
        ch_allele_codes = Channel.fromPath(checkFileExists("${input_dir}/02_genotype_processing/GM_allelecodes.csv", "allele codes file"))
    }

    // Continue only if we have genotype data
    if (params.finalreport_files || params.resume_from) {

        // MODULE 3: Control File Generation
        if (shouldRunStep('control', params.resume_from)) {
            log.info "Running GENERATE_CONTROL_FILE"

            GENERATE_CONTROL_FILE(
                ch_pheno,
                ch_covar,
                ch_geno_files,
                ch_allele_codes,
                ch_study_prefix
            )

            ch_control_file = GENERATE_CONTROL_FILE.out.control_file
            ch_founder_genos = GENERATE_CONTROL_FILE.out.founder_genos
            ch_genetic_maps = GENERATE_CONTROL_FILE.out.genetic_maps
            ch_physical_maps = GENERATE_CONTROL_FILE.out.physical_maps

            GENERATE_CONTROL_FILE.out.control_file.view { "Control file created: $it" }
        } else {
            log.info "Skipping GENERATE_CONTROL_FILE - loading from existing files"
            ch_control_file = Channel.fromPath(checkFileExists("${input_dir}/03_control_file_generation/${params.study_prefix}_control.json", "control file"))
            ch_founder_genos = Channel.fromPath("${input_dir}/03_control_file_generation/GM_foundergeno*.csv").collect()
            ch_genetic_maps = Channel.fromPath("${input_dir}/03_control_file_generation/GM_gmap*.csv").collect()
            ch_physical_maps = Channel.fromPath("${input_dir}/03_control_file_generation/GM_pmap*.csv").collect()
        }

        // MODULE 4: Cross2 Object Creation
        if (shouldRunStep('cross2', params.resume_from)) {
            log.info "Running CREATE_CROSS2_OBJECT"

            CREATE_CROSS2_OBJECT(
                ch_control_file,
                ch_founder_genos.mix(ch_genetic_maps, ch_physical_maps, ch_pheno, ch_covar, ch_geno_files).collect(),
                ch_study_prefix
            )

            ch_cross2_object = CREATE_CROSS2_OBJECT.out.cross2_object
            CREATE_CROSS2_OBJECT.out.validation_report.view { "Cross2 validation report: $it" }
        } else {
            log.info "Skipping CREATE_CROSS2_OBJECT - loading from existing files"
            ch_cross2_object = Channel.fromPath(checkFileExists("${input_dir}/04_cross2_creation/${params.study_prefix}_cross2.rds", "cross2 object"))
        }

        // MODULE 5: Genome Scan Preparation - set up channels based on resume logic
        if (shouldRunStep('prepare_scan', params.resume_from)) {
            PREPARE_GENOME_SCAN_SETUP(
                ch_cross2_object,
                ch_study_prefix
            )
            // Note: Removed .view operation to prevent channel consumption issues
            // PREPARE_GENOME_SCAN_SETUP.out.setup_report.view { "Genome scan preparation report: $it" }

            // Use outputs directly from Module 5 process
            ch_genoprob = PREPARE_GENOME_SCAN_SETUP.out.genoprob
            ch_alleleprob = PREPARE_GENOME_SCAN_SETUP.out.alleleprob
            ch_kinship_loco = PREPARE_GENOME_SCAN_SETUP.out.kinship_loco
            ch_genetic_map = PREPARE_GENOME_SCAN_SETUP.out.genetic_map

            // Make GENOME_SCAN_SETUP depend on Module 5 completion using a different output
            ch_cross2_for_setup = PREPARE_GENOME_SCAN_SETUP.out.setup_report
                .combine(ch_cross2_object)
                .map { report, cross2 -> cross2 }
        } else {
            // Load from existing files when Module 5 skipped
            ch_genoprob = Channel.fromPath(checkFileExists("${input_dir}/05_genome_scan_preparation/${params.study_prefix}_genoprob.rds", "genotype probabilities"))
            ch_alleleprob = Channel.fromPath(checkFileExists("${input_dir}/05_genome_scan_preparation/${params.study_prefix}_alleleprob.rds", "allele probabilities"))
            ch_kinship_loco = Channel.fromPath(checkFileExists("${input_dir}/05_genome_scan_preparation/${params.study_prefix}_kinship_loco.rds", "kinship matrices"))
            ch_genetic_map = Channel.fromPath(checkFileExists("${input_dir}/05_genome_scan_preparation/${params.study_prefix}_genetic_map.rds", "genetic map"))

            // No dependency needed when Module 5 skipped
            ch_cross2_for_setup = ch_cross2_object
        }

        // Check if genome scan should run or use existing results
        if (shouldRunStep('genome_scan', params.resume_from)) {
            // MODULE 6: Genome Scan Setup - now properly depends on Module 5 completion
            GENOME_SCAN_SETUP(
                ch_cross2_for_setup,
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            GENOME_SCAN_SETUP.out.summary.view { "Chunking summary: $it" }

            log.info "Running GENOME_SCAN_BATCH processes"

            // Create batch IDs from batch file
            GENOME_SCAN_SETUP.out.batch_file
                .splitText()
                .map { it.trim() }
                .filter { !it.startsWith('batch_id') }
                .map { it.split('\t')[0] }
                .unique()
                .set { ch_batch_ids }

            // Process batches in parallel
            GENOME_SCAN_BATCH(
                ch_batch_ids,
                ch_cross2_object,
                ch_genoprob,
                ch_alleleprob,
                ch_kinship_loco,
                GENOME_SCAN_SETUP.out.chunk_file,
                GENOME_SCAN_SETUP.out.batch_file,
                ch_study_prefix,
                Channel.value(params.lod_threshold),
                Channel.value(params.interactive_covar)
            )

            // Check if we need to run COMBINE or if batch results already exist
            def batch_results_dir = "${params.outdir}/06_qtl_analysis/batches"
            def combined_results_file = "${params.outdir}/06_qtl_analysis/${params.study_prefix}_scan_results.rds"

            if (file(combined_results_file).exists()) {
                log.info "Combined results already exist - skipping COMBINE_BATCH_RESULTS"
                ch_scan_results = Channel.fromPath(combined_results_file)
                ch_peaks = Channel.fromPath("${params.outdir}/06_qtl_analysis/${params.study_prefix}_all_peaks.csv")
                ch_filtered_phenotypes = Channel.fromPath("${params.outdir}/06_qtl_analysis/${params.study_prefix}_filtered_phenotypes.txt")
            } else {
                log.info "Combining batch results from published files"

                // Instead of collecting from channels, read from published batch directory
                // This avoids the .collect() hang issue when batch jobs are cached
                Channel.fromPath("${batch_results_dir}/*_batch*_results.rds").collect().set { ch_batch_results_files }
                Channel.fromPath("${batch_results_dir}/*_batch*_peaks.csv").collect().set { ch_batch_peaks_files }
                Channel.fromPath("${batch_results_dir}/*_batch*_filtered_phenos.txt").collect().set { ch_batch_phenos_files }
                Channel.fromPath("${batch_results_dir}/*_batch*_log.txt").collect().set { ch_batch_log_files }

                COMBINE_BATCH_RESULTS(
                    ch_batch_results_files,
                    ch_batch_peaks_files,
                    ch_batch_phenos_files,
                    ch_batch_log_files,
                    ch_study_prefix,
                    Channel.value(params.lod_threshold)
                )

                ch_scan_results = COMBINE_BATCH_RESULTS.out.scan_results
                ch_peaks = COMBINE_BATCH_RESULTS.out.peaks
                ch_filtered_phenotypes = COMBINE_BATCH_RESULTS.out.filtered_phenotypes

                COMBINE_BATCH_RESULTS.out.batch_summary.view { "Batch processing summary: $it" }
                COMBINE_BATCH_RESULTS.out.peaks.view { "Peaks found: $it" }
            }

        } else {
            log.info "Skipping GENOME_SCAN_BATCH - loading from existing files"

            // Load existing combined results
            ch_scan_results = Channel.fromPath(checkFileExists("${input_dir}/06_qtl_analysis/${params.study_prefix}_scan_results.rds", "scan results"))
            ch_peaks = Channel.fromPath(checkFileExists("${input_dir}/06_qtl_analysis/${params.study_prefix}_all_peaks.csv", "peaks file"))
            ch_filtered_phenotypes = Channel.fromPath(checkFileExists("${input_dir}/06_qtl_analysis/${params.study_prefix}_filtered_phenotypes.txt", "filtered phenotypes"))
        }

        // MODULE 6b: Optional Regional QTL Filtering (before permutation testing)
        def filtered_phenos_for_perm
        def regional_filter_file = "${input_dir}/06b_regional_filtering/${params.study_prefix}_regional_filtered_phenotypes.txt"

        if (params.resume_from == 'regional_filter') {
            // Regional filter step explicitly requested
            if (params.qtl_region == null) {
                log.error "ERROR: --resume_from regional_filter requires --qtl_region parameter"
                exit 1
            }

            log.info "Running regional filtering on loaded peaks: ${params.qtl_region}"

            FILTER_PEAKS_BY_REGION(
                ch_peaks,
                ch_genetic_map,
                Channel.fromPath(params.gtf_file),
                ch_study_prefix,
                Channel.value(params.qtl_region)
            )

            FILTER_PEAKS_BY_REGION.out.filter_report.view { "Regional filtering report: $it" }
            filtered_phenos_for_perm = FILTER_PEAKS_BY_REGION.out.filtered_phenotypes

        } else if (file(regional_filter_file).exists()) {
            // Existing regional filter found - use it automatically
            log.info "Using existing regional filtered phenotypes from: ${regional_filter_file}"
            filtered_phenos_for_perm = Channel.fromPath(regional_filter_file)

        } else if (params.qtl_region != null) {
            // Regional filtering enabled via parameter but no existing filter
            if (shouldRunStep('genome_scan', params.resume_from)) {
                // Genome scan just ran - filter using its output
                log.info "Filtering QTL peaks to genomic region: ${params.qtl_region}"

                FILTER_PEAKS_BY_REGION(
                    COMBINE_BATCH_RESULTS.out.peaks,
                    ch_genetic_map,
                    Channel.fromPath(params.gtf_file),
                    ch_study_prefix,
                    Channel.value(params.qtl_region)
                )

                FILTER_PEAKS_BY_REGION.out.filter_report.view { "Regional filtering report: $it" }
                filtered_phenos_for_perm = FILTER_PEAKS_BY_REGION.out.filtered_phenotypes

            } else {
                log.info "No existing regional filter - using all LOD-filtered phenotypes"
                filtered_phenos_for_perm = ch_filtered_phenotypes
            }
        } else {
            log.info "No regional filtering applied - using all LOD-filtered phenotypes"
            filtered_phenos_for_perm = ch_filtered_phenotypes
        }

        // MODULE 7: Chunked Permutation Testing (DO_Pipe Approach)
        if (shouldRunStep('permutation', params.resume_from)) {
            log.info "Running CHUNKED_PERMUTATION_TESTING (DO_Pipe approach)"

            // Make permutation depend on genome scan completion if genome scan ran
            if (shouldRunStep('genome_scan', params.resume_from)) {
                // Use genome scan output as dependency trigger
                ch_cross2_for_perm = COMBINE_BATCH_RESULTS.out.batch_summary
                    .combine(ch_cross2_object)
                    .map { summary, cross2 -> cross2 }
            } else {
                ch_cross2_for_perm = ch_cross2_object
            }

            CHUNKED_PERMUTATION_TESTING(
                ch_cross2_for_perm,
                ch_genoprob,
                ch_kinship_loco,
                filtered_phenos_for_perm,
                ch_study_prefix,
                Channel.value(params.lod_threshold),
                Channel.value(params.run_perm_benchmark),
                Channel.value(params.perm_per_chunk),
                Channel.value(params.chunks_per_batch)
            )

            ch_perm_results = CHUNKED_PERMUTATION_TESTING.out.permutation_matrix
            ch_thresholds = CHUNKED_PERMUTATION_TESTING.out.permutation_thresholds
            ch_filtered_cross2 = CHUNKED_PERMUTATION_TESTING.out.filtered_cross2

        } else if (params.resume_from == 'perm_aggregate') {
            log.info "Running only PERM_AGGREGATE - using existing batch results"

            // Collect batch result files from the published batches directory
            ch_batch_results = Channel.fromPath("${input_dir}/07_permutation_testing/batches/${params.study_prefix}_*_batch_*.rds")
                .collect()

            PERM_AGGREGATE(
                ch_batch_results,
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            ch_perm_results = PERM_AGGREGATE.out.full_perm_matrix
            ch_thresholds = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_cross2 = Channel.fromPath(checkFileExists("${input_dir}/07_permutation_testing/${params.study_prefix}_filtered_cross2.rds", "filtered cross2"))

        } else {
            log.info "Skipping CHUNKED_PERMUTATION_TESTING - loading from existing files"

            // Load existing permutation results (updated file names for chunked approach)
            ch_perm_results = Channel.fromPath(checkFileExists("${input_dir}/07_permutation_testing/${params.study_prefix}_fullPerm.rds", "permutation matrix"))
            ch_thresholds = Channel.fromPath(checkFileExists("${input_dir}/07_permutation_testing/${params.study_prefix}_permThresh.rds", "permutation thresholds"))
            ch_filtered_cross2 = Channel.fromPath(checkFileExists("${input_dir}/07_permutation_testing/${params.study_prefix}_filtered_cross2.rds", "filtered cross2"))
        }

        // Display results (conditional based on whether permutation was run)
        if (shouldRunStep('permutation', params.resume_from)) {
            CHUNKED_PERMUTATION_TESTING.out.setup_log.view { "Permutation setup log: $it" }
            CHUNKED_PERMUTATION_TESTING.out.aggregation_log.view { "Permutation aggregation log: $it" }
            CHUNKED_PERMUTATION_TESTING.out.permutation_thresholds.view { "Significance thresholds: $it" }
        }

        // MODULE 8: Significant QTL Identification
        if (shouldRunStep('significant_qtls', params.resume_from)) {
            IDENTIFY_SIGNIFICANT_QTLS(
                ch_filtered_cross2,
                ch_scan_results,
                ch_perm_results,
                ch_thresholds,
                ch_study_prefix
            )

            IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls.view { "Significant QTLs identified: $it" }
            IDENTIFY_SIGNIFICANT_QTLS.out.qtl_summary.view { "QTL summary report: $it" }

            ch_significant_qtls = IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls
        } else {
            log.info "Skipping IDENTIFY_SIGNIFICANT_QTLS - loading from existing files"
            ch_significant_qtls = Channel.fromPath(checkFileExists("${input_dir}/08_significant_qtls/${params.study_prefix}_significant_qtls.csv", "significant QTLs"))
        }

        // MODULE 9: QTL Viewer Integration (Optional - supports HPC and local deployment)
        // NOTE: Module 9 can be deployed on HPC or locally via Docker/Singularity
        // Use --run_qtlviewer to enable this module
        if (params.run_qtlviewer || shouldRunStep('qtlviewer', params.resume_from)) {
            log.info "Running QTL Viewer setup (supports HPC and local deployment)"
            log.info "Using genoprobs (not alleleprobs) for QTL Viewer compatibility"

            PREPARE_QTLVIEWER_DATA(
                ch_genoprob,         // Changed from ch_alleleprob
                ch_kinship_loco,
                ch_genetic_map,
                ch_cross2_object,    // Added for pmap extraction
                ch_scan_results,
                ch_significant_qtls,
                Channel.fromPath(params.gtf_file, checkIfExists: true),
                ch_study_prefix
            )

            PREPARE_QTLVIEWER_DATA.out.validation_report.view { "QTL Viewer conversion report: $it" }
            PREPARE_QTLVIEWER_DATA.out.instructions.view { "QTL Viewer instructions: $it" }
        } else {
            log.info "Skipping QTL Viewer setup - use --run_qtlviewer to enable (supports HPC and local deployment)"
        }

    } else {
        log.info "Skipping genotype-dependent modules - no FinalReport files specified"
    }
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
    Status     : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Time       : ${workflow.duration}
    Workdir    : ${workflow.workDir}
    Results    : ${params.outdir}
    Resume from: ${params.resume_from ?: 'beginning'}
    ========================================
    """.stripIndent()
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