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
params.resume_from = null  // Options: phenotype, genotype, control, cross2, prepare_scan, genome_scan, permutation, significant_qtls, qtlviewer

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
                            permutation     - Start from permutation testing
                            significant_qtls - Start from QTL identification
                            qtlviewer       - Start from QTL Viewer setup

    Other Arguments:
        --finalreport_files  Path to GeneSeek FinalReport file(s)
        --outdir            Output directory [default: results]
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
                         'genome_scan', 'permutation', 'significant_qtls', 'qtlviewer']

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
include { PERMUTATION_SETUP; PERMUTATION_BATCH; COMBINE_PERMUTATION_RESULTS } from './modules/07_permutation_testing.nf'
include { IDENTIFY_SIGNIFICANT_QTLS } from './modules/08_identify_significant_qtls.nf'
include { PREPARE_QTLVIEWER_DATA; SETUP_QTLVIEWER_DEPLOYMENT } from './modules/09_qtl_viewer.nf'

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
                    'genome_scan', 'permutation', 'significant_qtls', 'qtlviewer']

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

    log.info """
    ========================================
     QTL2_NF Pipeline Started
    ========================================
    study_prefix   : ${params.study_prefix}
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
        ch_pheno = Channel.fromPath(checkFileExists("${params.outdir}/01_phenotype_processing/${params.study_prefix}_pheno.csv", "phenotype file"))
        ch_covar = Channel.fromPath(checkFileExists("${params.outdir}/01_phenotype_processing/${params.study_prefix}_covar.csv", "covariate file"))
        ch_valid_samples = Channel.fromPath(checkFileExists("${params.outdir}/01_phenotype_processing/${params.study_prefix}_valid_samples.txt", "valid samples file"))
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
        ch_geno_files = Channel.fromPath("${params.outdir}/02_genotype_processing/${params.study_prefix}_geno*.csv").collect()
        ch_allele_codes = Channel.fromPath(checkFileExists("${params.outdir}/02_genotype_processing/GM_allelecodes.csv", "allele codes file"))
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
            ch_control_file = Channel.fromPath(checkFileExists("${params.outdir}/03_control_file_generation/${params.study_prefix}_control.json", "control file"))
            ch_founder_genos = Channel.fromPath("${params.outdir}/03_control_file_generation/GM_foundergeno*.csv").collect()
            ch_genetic_maps = Channel.fromPath("${params.outdir}/03_control_file_generation/GM_gmap*.csv").collect()
            ch_physical_maps = Channel.fromPath("${params.outdir}/03_control_file_generation/GM_pmap*.csv").collect()
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
            ch_cross2_object = Channel.fromPath(checkFileExists("${params.outdir}/04_cross2_creation/${params.study_prefix}_cross2.rds", "cross2 object"))
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

            // Make GENOME_SCAN_SETUP depend on Module 5 completion using a different output
            ch_cross2_for_setup = PREPARE_GENOME_SCAN_SETUP.out.setup_report
                .combine(ch_cross2_object)
                .map { report, cross2 -> cross2 }
        } else {
            // Load from existing files when Module 5 skipped
            ch_genoprob = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_genoprob.rds", "genotype probabilities"))
            ch_alleleprob = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_alleleprob.rds", "allele probabilities"))
            ch_kinship_loco = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_kinship_loco.rds", "kinship matrices"))

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
                Channel.value(params.lod_threshold)
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

            ch_scan_results = COMBINE_BATCH_RESULTS.out.scan_results
            ch_peaks = COMBINE_BATCH_RESULTS.out.peaks
            ch_filtered_phenotypes = COMBINE_BATCH_RESULTS.out.filtered_phenotypes

        } else {
            log.info "Skipping GENOME_SCAN_BATCH - loading from existing files"

            // Load existing combined results
            ch_scan_results = Channel.fromPath(checkFileExists("${params.outdir}/06_qtl_analysis/${params.study_prefix}_scan_results.rds", "scan results"))
            ch_peaks = Channel.fromPath(checkFileExists("${params.outdir}/06_qtl_analysis/${params.study_prefix}_all_peaks.csv", "peaks file"))
            ch_filtered_phenotypes = Channel.fromPath(checkFileExists("${params.outdir}/06_qtl_analysis/${params.study_prefix}_filtered_phenotypes.txt", "filtered phenotypes"))
        }

        // Display results (conditional based on whether batch was run)
        if (shouldRunStep('genome_scan', params.resume_from)) {
            COMBINE_BATCH_RESULTS.out.batch_summary.view { "Batch processing summary: $it" }
            COMBINE_BATCH_RESULTS.out.peaks.view { "Peaks found: $it" }
        }

        // MODULE 7: Permutation Testing (Batched Approach)
        if (shouldRunStep('permutation', params.resume_from)) {
            log.info "Running PERMUTATION_BATCH processes"

            // Setup permutation batching
            PERMUTATION_SETUP(
                ch_filtered_phenotypes,
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            PERMUTATION_SETUP.out.summary.view { "Permutation chunking summary: $it" }

            // Create batch IDs from batch file
            PERMUTATION_SETUP.out.batch_file
                .splitText()
                .map { it.trim() }
                .filter { !it.startsWith('batch_id') }
                .map { it.split('\t')[0] }
                .unique()
                .set { ch_perm_batch_ids }

            // Process permutation batches in parallel
            PERMUTATION_BATCH(
                ch_perm_batch_ids,
                ch_cross2_object,
                ch_genoprob,
                ch_kinship_loco,
                PERMUTATION_SETUP.out.chunk_file,
                PERMUTATION_SETUP.out.batch_file,
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            // Combine all permutation batch results
            COMBINE_PERMUTATION_RESULTS(
                PERMUTATION_BATCH.out.batch_perm_results.collect(),
                PERMUTATION_BATCH.out.batch_thresholds.collect(),
                PERMUTATION_BATCH.out.batch_log.collect(),
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            ch_perm_results = COMBINE_PERMUTATION_RESULTS.out.perm_results
            ch_thresholds = COMBINE_PERMUTATION_RESULTS.out.thresholds

        } else {
            log.info "Skipping PERMUTATION_BATCH - loading from existing files"

            // Load existing permutation results
            ch_perm_results = Channel.fromPath(checkFileExists("${params.outdir}/07_permutation_testing/${params.study_prefix}_permutation_results.rds", "permutation results"))
            ch_thresholds = Channel.fromPath(checkFileExists("${params.outdir}/07_permutation_testing/${params.study_prefix}_significance_thresholds.csv", "significance thresholds"))
        }

        // Display results (conditional based on whether permutation was run)
        if (shouldRunStep('permutation', params.resume_from)) {
            COMBINE_PERMUTATION_RESULTS.out.validation_report.view { "Permutation test report: $it" }
            COMBINE_PERMUTATION_RESULTS.out.thresholds.view { "Significance thresholds: $it" }
        }

        // MODULE 8: Significant QTL Identification
        IDENTIFY_SIGNIFICANT_QTLS(
            ch_cross2_object,
            ch_scan_results,
            ch_perm_results,
            ch_thresholds,
            ch_study_prefix
        )

        IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls.view { "Significant QTLs identified: $it" }
        IDENTIFY_SIGNIFICANT_QTLS.out.qtl_summary.view { "QTL summary report: $it" }

        // MODULE 9: QTL Viewer Integration
        if (shouldRunStep('prepare_scan', params.resume_from)) {
            ch_genetic_map = PREPARE_GENOME_SCAN_SETUP.out.genetic_map
        } else {
            ch_genetic_map = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_genetic_map.rds", "genetic map"))
        }

        PREPARE_QTLVIEWER_DATA(
            ch_cross2_object,
            ch_genoprob,
            ch_kinship_loco,
            ch_genetic_map,
            ch_scan_results,
            IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls,
            ch_study_prefix
        )

        SETUP_QTLVIEWER_DEPLOYMENT(
            PREPARE_QTLVIEWER_DATA.out.qtlviewer_data,
            ch_study_prefix
        )

        PREPARE_QTLVIEWER_DATA.out.validation_report.view { "QTL Viewer conversion report: $it" }
        SETUP_QTLVIEWER_DEPLOYMENT.out.instructions.view { "QTL Viewer instructions: $it" }

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