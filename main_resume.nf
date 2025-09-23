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
include { PREPARE_GENOME_SCAN } from './modules/05_prepare_genome_scan.nf'
include { GENOME_SCAN; GENOME_SCAN_CHUNK; COMBINE_SCAN_RESULTS } from './modules/06_qtl_analysis.nf'
include { PERMUTATION_TEST } from './modules/07_permutation_testing.nf'
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

        // MODULE 5: Genome Scan Preparation
        if (shouldRunStep('prepare_scan', params.resume_from)) {
            log.info "Running PREPARE_GENOME_SCAN"

            PREPARE_GENOME_SCAN(
                ch_cross2_object,
                ch_study_prefix
            )

            ch_genoprob = PREPARE_GENOME_SCAN.out.genoprob
            ch_kinship = PREPARE_GENOME_SCAN.out.kinship
            ch_genetic_map = PREPARE_GENOME_SCAN.out.genetic_map

            PREPARE_GENOME_SCAN.out.validation_report.view { "Genome scan preparation report: $it" }
        } else {
            log.info "Skipping PREPARE_GENOME_SCAN - loading from existing files"
            ch_genoprob = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_genoprob.rds", "genotype probabilities"))
            ch_kinship = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_kinship_loco.rds", "kinship matrices"))
            ch_genetic_map = Channel.fromPath(checkFileExists("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_genetic_map.rds", "genetic map"))
        }

        // MODULE 6: HPC Array Genome Scanning
        if (shouldRunStep('genome_scan', params.resume_from)) {
            log.info "Running HPC Array Genome Scan with massive parallelization"

            // Step 1: Coordinator - Create chunk information
            GENOME_SCAN(
                ch_cross2_object,
                ch_genoprob,
                ch_kinship,
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            // Step 2: Parse chunk info and create array jobs
            chunk_info = GENOME_SCAN.out.chunk_info
                .splitCsv(header: true)
                .map { row ->
                    tuple(row.chunk_id as Integer, row.pheno_start as Integer, row.pheno_end as Integer)
                }

            // Step 3: Launch array jobs - one per chunk
            GENOME_SCAN_CHUNK(
                ch_cross2_object.first(),
                ch_genoprob.first(),
                ch_kinship.first(),
                ch_study_prefix.first(),
                chunk_info.map { it[0] },  // chunk_id
                chunk_info.map { it[1] },  // pheno_start
                chunk_info.map { it[2] }   // pheno_end
            )

            // Step 4: Combine all chunk results
            COMBINE_SCAN_RESULTS(
                GENOME_SCAN_CHUNK.out.chunk_results.collect(),
                ch_cross2_object.first(),
                ch_study_prefix.first(),
                Channel.value(params.lod_threshold)
            )

            ch_scan_results = COMBINE_SCAN_RESULTS.out.scan_results
            ch_filtered_phenotypes = COMBINE_SCAN_RESULTS.out.filtered_phenotypes

            GENOME_SCAN.out.setup_report.view { "Array setup report: $it" }
            COMBINE_SCAN_RESULTS.out.validation_report.view { "HPC array results: $it" }
            COMBINE_SCAN_RESULTS.out.peaks.view { "Preliminary peaks found: $it" }
        } else {
            log.info "Skipping HPC Array Genome Scan - loading from existing files"
            ch_scan_results = Channel.fromPath(checkFileExists("${params.outdir}/06_qtl_analysis/${params.study_prefix}_scan_results.rds", "scan results"))
            ch_filtered_phenotypes = Channel.fromPath(checkFileExists("${params.outdir}/06_qtl_analysis/${params.study_prefix}_filtered_phenotypes.txt", "filtered phenotypes"))
        }

        // MODULE 7: Permutation Testing
        if (shouldRunStep('permutation', params.resume_from)) {
            log.info "Running PERMUTATION_TEST"

            PERMUTATION_TEST(
                ch_cross2_object,
                ch_genoprob,
                ch_kinship,
                ch_filtered_phenotypes,
                ch_study_prefix,
                Channel.value(params.lod_threshold)
            )

            ch_perm_results = PERMUTATION_TEST.out.perm_results
            ch_thresholds = PERMUTATION_TEST.out.thresholds

            PERMUTATION_TEST.out.validation_report.view { "Permutation test report: $it" }
            PERMUTATION_TEST.out.thresholds.view { "Significance thresholds: $it" }
        } else {
            log.info "Skipping PERMUTATION_TEST - loading from existing files"
            ch_perm_results = Channel.fromPath(checkFileExists("${params.outdir}/07_permutation_testing/${params.study_prefix}_permutation_results.rds", "permutation results"))
            ch_thresholds = Channel.fromPath(checkFileExists("${params.outdir}/07_permutation_testing/${params.study_prefix}_significance_thresholds.csv", "significance thresholds"))
        }

        // MODULE 8: Significant QTL Identification
        if (shouldRunStep('significant_qtls', params.resume_from)) {
            log.info "Running IDENTIFY_SIGNIFICANT_QTLS"

            IDENTIFY_SIGNIFICANT_QTLS(
                ch_cross2_object,
                ch_scan_results,
                ch_perm_results,
                ch_thresholds,
                ch_study_prefix
            )

            ch_significant_qtls = IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls

            IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls.view { "Significant QTLs identified: $it" }
            IDENTIFY_SIGNIFICANT_QTLS.out.qtl_summary.view { "QTL summary report: $it" }
        } else {
            log.info "Skipping IDENTIFY_SIGNIFICANT_QTLS - loading from existing files"
            ch_significant_qtls = Channel.fromPath(checkFileExists("${params.outdir}/08_significant_qtls/${params.study_prefix}_significant_qtls.csv", "significant QTLs"))
        }

        // MODULE 9: QTL Viewer Integration
        if (shouldRunStep('qtlviewer', params.resume_from)) {
            log.info "Running QTL Viewer setup"

            PREPARE_QTLVIEWER_DATA(
                ch_cross2_object,
                ch_genoprob,
                ch_kinship,
                ch_genetic_map,
                ch_scan_results,
                ch_significant_qtls,
                ch_study_prefix
            )

            SETUP_QTLVIEWER_DEPLOYMENT(
                PREPARE_QTLVIEWER_DATA.out.qtlviewer_data,
                ch_study_prefix
            )

            PREPARE_QTLVIEWER_DATA.out.validation_report.view { "QTL Viewer conversion report: $it" }
            SETUP_QTLVIEWER_DEPLOYMENT.out.instructions.view { "QTL Viewer instructions: $it" }
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