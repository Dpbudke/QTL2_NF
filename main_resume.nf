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
params.outdir = null          // Resolved from analysis_type/study_type (see config auto-compute block)
params.auto_prefix_samples = false
params.test_mode = false
params.lod_threshold = 7.0
params.sample_filter = null
params.interactive_covar = null
params.exclude_covariates = null
params.help = false

// Analysis type preset (mirrors main.nf — auto-configures outdir, interactive_covar, etc.)
params.analysis_type = null   // 'Diet_additive', 'Diet_interactive', 'AIN76_only', 'HC_only'
params.study_type = null      // Set via config default ('eQTL') or override (e.g. 'LCMS')

// Resume parameters
params.resume_from = null  // Options: phenotype, genotype, control, cross2, prepare_scan, genome_scan, permutation, perm_aggregate, significant_qtls, visualize, timbr
params.input_dir = null    // Directory to load cached results from (defaults to outdir if not specified)

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
                            visualize       - Start from QTL visualization
                            timbr           - Start from TIMBR allelic series analysis

    Other Arguments:
        --finalreport_files  Path to GeneSeek FinalReport file(s)
        --outdir            Output directory for new results [default: results]
        --input_dir         Directory to load cached results from (defaults to outdir if not specified)
        --lod_threshold     LOD threshold for filtering [default: 7.0]
        --sample_filter     JSON filter for sample subsetting

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
                         'genome_scan', 'regional_filter', 'permutation', 'perm_aggregate', 'significant_qtls', 'visualize', 'timbr']

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
    ANALYSIS TYPE PRESET CONFIGURATION (mirrors main.nf)
========================================================================================
*/

def validAnalysisTypes = ['Diet_additive', 'Diet_interactive', 'AIN76_only', 'HC_only']

if (params.analysis_type && !validAnalysisTypes.contains(params.analysis_type)) {
    log.error """
    Invalid --analysis_type: '${params.analysis_type}'
    Valid options: ${validAnalysisTypes.join(', ')}
    """.stripIndent()
    exit 1
}

// Note: params.outdir is already resolved by the config auto-compute block when analysis_type is set.
// The effective_* variables below are needed for interactive_covar / sample_filter / exclude_covariates
// which the config block does not set (those are only used at runtime, not in publishDir paths).

def effective_outdir = params.outdir ?: 'results'

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
        resume_from        : ${params.resume_from ?: 'beginning (full run)'}
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
include { CHUNKED_PERMUTATION_TESTING; PERM_AGGREGATE } from './modules/07_permutation_testing.nf'
include { FIND_PEAKS_CHR; GATHER_PEAKS } from './modules/08_identify_significant_qtls.nf'
include { VISUALIZE_QTLS } from './modules/09_visualize.nf'
include { TIMBR_ANALYSIS } from './modules/10_timbr.nf'
include { CLASSIFY_CIS_TRANS_EQTLS } from './modules/analyses/classify_cis_trans_eqtls.nf'

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
                    'genome_scan', 'regional_filter', 'permutation', 'perm_aggregate', 'significant_qtls', 'visualize', 'timbr']

    def currentIndex = stepOrder.indexOf(currentStep)
    def resumeIndex = stepOrder.indexOf(resumeFrom)

    return currentIndex >= resumeIndex
}

// Function to check if we should load files (step was skipped) or use process outputs (step ran)
def shouldLoadFromFiles(currentStep, resumeFrom) {
    if (!resumeFrom) return false

    def stepOrder = ['phenotype', 'genotype', 'control', 'cross2', 'prepare_scan',
                    'genome_scan', 'permutation', 'significant_qtls', 'visualize']

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
            Channel.value(effective_sample_filter ?: "null"),
            Channel.value(effective_exclude_covariates ?: "null"),
            Channel.value(params.single_sex ?: "null")
        )

        ch_pheno = PHENOTYPE_PROCESS.out.pheno
        ch_covar = PHENOTYPE_PROCESS.out.covar
        ch_valid_samples = PHENOTYPE_PROCESS.out.valid_samples
        ch_valid_samples_original = PHENOTYPE_PROCESS.out.valid_samples_original

        PHENOTYPE_PROCESS.out.validation_report.view { "Validation report: $it" }
    } else {
        log.info "Skipping PHENOTYPE_PROCESS - loading from existing files"
        ch_pheno = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_pheno.csv", "phenotype file"))
        ch_covar = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_covar.csv", "covariate file"))
        ch_valid_samples = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_valid_samples.txt", "valid samples file"))
        ch_valid_samples_original = Channel.fromPath(checkFileExists("${input_dir}/01_phenotype_processing/${params.study_prefix}_valid_samples_original.txt", "original valid samples file"))
    }

    // MODULE 2: Genotype Processing
    if (shouldRunStep('genotype', params.resume_from) && params.finalreport_files) {
        log.info "Running GENOTYPE_PROCESS"
        ch_finalreport = Channel.fromPath(params.finalreport_files, checkIfExists: true).collect()
        ch_allele_codes = Channel.fromPath("${projectDir}/Data/GM/GM_allelecodes.csv", checkIfExists: true)

        GENOTYPE_PROCESS(
            ch_finalreport,
            ch_valid_samples,
            ch_valid_samples_original,
            ch_allele_codes,
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
                Channel.value(effective_interactive_covar ?: "null")
            )

            // Connect COMBINE_BATCH_RESULTS directly to GENOME_SCAN_BATCH channel outputs
            // (mirrors main.nf — reactive channel graph ensures combine runs after all batches finish)
            COMBINE_BATCH_RESULTS(
                GENOME_SCAN_BATCH.out.batch_results.collect(),
                GENOME_SCAN_BATCH.out.batch_peaks.collect(),
                GENOME_SCAN_BATCH.out.batch_filtered_phenos.collect(),
                GENOME_SCAN_BATCH.out.batch_log.collect(),
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            ch_scan_results = COMBINE_BATCH_RESULTS.out.scan_results
            ch_peaks = COMBINE_BATCH_RESULTS.out.peaks
            ch_filtered_phenotypes = COMBINE_BATCH_RESULTS.out.filtered_phenotypes

            COMBINE_BATCH_RESULTS.out.batch_summary.view { "Batch processing summary: $it" }
            COMBINE_BATCH_RESULTS.out.peaks.view { "Peaks found: $it" }

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
            log.info "Running only PERM_AGGREGATE - all batches assumed complete"

            // Pass the phenotype list as the trigger (confirms PERM_SETUP ran and
            // batches exist). PERM_AGGREGATE now reads batch files directly from
            // the published batch directory rather than staging them as inputs.
            ch_agg_trigger = Channel.fromPath(
                checkFileExists("${input_dir}/07_permutation_testing/${params.study_prefix}_phenotype_list.txt",
                                "phenotype list"))

            PERM_AGGREGATE(
                ch_agg_trigger,
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

        // MODULE 8: Significant QTL Identification (scatter-gather)
        // Channel.of() drives 20 parallel SLURM jobs (one per DO chromosome).
        // Each job loads the full scan matrix and subsets to its chromosome internally.
        // This avoids the glob/flatten/join pattern that was unreliable in resume mode.
        if (shouldRunStep('significant_qtls', params.resume_from)) {
            ch_chrs = Channel.of("1","2","3","4","5","6","7","8","9","10",
                                  "11","12","13","14","15","16","17","18","19","X")

            FIND_PEAKS_CHR(ch_chrs, ch_scan_results.first(), ch_filtered_cross2.first(), ch_thresholds.first())

            GATHER_PEAKS(
                FIND_PEAKS_CHR.out.chr_peaks.collect(),
                ch_filtered_cross2,
                ch_perm_results,
                ch_thresholds,
                ch_study_prefix
            )

            GATHER_PEAKS.out.significant_qtls.view { "Significant QTLs identified: $it" }
            GATHER_PEAKS.out.qtl_summary.view { "QTL summary report: $it" }

            ch_significant_qtls = GATHER_PEAKS.out.significant_qtls
        } else {
            log.info "Skipping QTL identification - loading from existing files"
            ch_significant_qtls = Channel.fromPath(checkFileExists("${input_dir}/08_significant_qtls/${params.study_prefix}_significant_qtls.csv", "significant QTLs"))
        }

        // eQTL CLASSIFICATION: Classify cis vs trans eQTLs (eQTL study type only)
        if (params.study_type == 'eQTL') {
            log.info "Running eQTL cis/trans classification (study_type = eQTL)"
            CLASSIFY_CIS_TRANS_EQTLS(
                ch_significant_qtls,
                Channel.fromPath(params.gtf_file),
                ch_filtered_cross2
            )
            CLASSIFY_CIS_TRANS_EQTLS.out.summary.view { "eQTL classification summary: $it" }
        } else {
            log.info "Skipping eQTL classification (study_type = ${params.study_type})"
        }

        // MODULE 9: QTL Visualizations
        // NOTE: This module generates individual PNG plots for significant QTLs
        if (shouldRunStep('visualize', params.resume_from)) {
            log.info "Running QTL visualization (plot_coefCC for 99% significant QTLs)"

            VISUALIZE_QTLS(
                ch_alleleprob,
                ch_scan_results,
                ch_filtered_cross2,
                ch_significant_qtls,
                ch_study_prefix,
                ch_kinship_loco,
                Channel.value(effective_interactive_covar ?: "null")
            )

            VISUALIZE_QTLS.out.validation_report.view { "QTL visualization report: $it" }
        } else {
            log.info "Skipping QTL visualization - will resume from this step if needed"
        }

        // MODULE 10: TIMBR Allelic Series Analysis
        if (shouldRunStep('timbr', params.resume_from)) {
            // Load significant QTLs from file if resuming from timbr
            def ch_sig_qtls_timbr = (params.resume_from == 'timbr')
                ? Channel.fromPath(checkFileExists("${input_dir}/08_significant_qtls/${params.study_prefix}_significant_qtls.csv", "significant QTLs for TIMBR"))
                : ch_significant_qtls

            // Gate TIMBR on Module 9 completion so visualization always finishes first.
            // If Module 9 ran, use its validation_report as a dependency trigger by
            // combining it with ch_filtered_cross2 and then dropping the report path.
            // If Module 9 was skipped, pass ch_filtered_cross2 directly.
            def ch_cross2_for_timbr
            if (shouldRunStep('visualize', params.resume_from)) {
                ch_cross2_for_timbr = VISUALIZE_QTLS.out.validation_report
                    .combine(ch_filtered_cross2)
                    .map { _report, cross2 -> cross2 }
            } else {
                ch_cross2_for_timbr = ch_filtered_cross2
            }

            TIMBR_ANALYSIS(
                ch_cross2_for_timbr,
                ch_genoprob,
                ch_genetic_map,
                ch_sig_qtls_timbr,
                ch_study_prefix,
                Channel.value(params.timbr_sig_level),
                Channel.value(params.timbr_qtls_per_batch),
                Channel.value(params.timbr_samples)
            )
            TIMBR_ANALYSIS.out.master_summary.view { "TIMBR master summary: $it" }
            TIMBR_ANALYSIS.out.run_report.view     { "TIMBR run report: $it" }
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