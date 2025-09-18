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
        --finalreport_file  Path to GeneSeek FinalReport file (enables full QTL pipeline)
        --outdir            Output directory [default: results]
        --auto_prefix_samples    Automatically prefix sample IDs [default: false]
        --test_mode         Positive control: chr 2 only, coat_color as phenotype [default: false]
        --lod_threshold     LOD threshold for filtering QTLs before permutation testing [default: 7.0]
        --sample_filter     JSON filter for sample subsetting by covariates (e.g., '{"Sex":["male"],"Diet":["hc"]}') [default: null]
    
    Pipeline Modules:
        Module 1: Phenotype processing and validation
        Module 2: Genotype processing (requires --finalreport_file)
        Module 3: Control file generation and cross2 object creation
        Module 4: Genome scanning and permutation testing (1000 permutations)
        Module 5: QTL Viewer data preparation and Docker deployment
    
    Example:
        nextflow run main.nf \\
            --phenotype_file Data/formatted_phenotypes_for_nextflow.csv \\
            --finalreport_file Data/FinalReport.txt \\
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
    IMPORT MODULES
========================================================================================
*/

include { PHENOTYPE_PROCESS } from './modules/phenotype_process.nf'
include { GENOTYPE_PROCESS  } from './modules/genotype_process.nf'
include { GENERATE_CONTROL_FILE; CREATE_CROSS2_OBJECT } from './modules/control_cross2.nf'
include { PREPARE_GENOME_SCAN; GENOME_SCAN; PERMUTATION_TEST; IDENTIFY_SIGNIFICANT_QTLS } from './modules/scan_perm.nf'
include { PREPARE_QTLVIEWER_DATA; SETUP_QTLVIEWER_DEPLOYMENT } from './modules/qtl_viewer.nf'

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
    outdir         : ${params.outdir}
    test_mode      : ${params.test_mode ? 'true (chr 2 only)' : 'false'}
    lod_threshold  : ${params.lod_threshold}
    ========================================
    """.stripIndent()

    // Create input channels
    ch_phenotype_file = Channel.fromPath(params.phenotype_file, checkIfExists: true)
    ch_study_prefix   = Channel.value(params.study_prefix)
    
    // MODULE 1: Phenotype Processing
    PHENOTYPE_PROCESS(
        ch_phenotype_file,
        ch_study_prefix,
        Channel.value(params.sample_filter ?: "null")
    )
    
    // MODULE 2: Genotype Processing (if FinalReport file provided)
    if (params.finalreport_file) {
        ch_finalreport = Channel.fromPath(params.finalreport_file, checkIfExists: true)
        
        GENOTYPE_PROCESS(
            ch_finalreport,
            PHENOTYPE_PROCESS.out.valid_samples,
            ch_study_prefix
        )
        
        // MODULE 3: Control File and Cross2 Object Creation
        GENERATE_CONTROL_FILE(
            PHENOTYPE_PROCESS.out.pheno,
            PHENOTYPE_PROCESS.out.covar,
            GENOTYPE_PROCESS.out.geno_files.collect(),
            GENOTYPE_PROCESS.out.allele_codes,
            ch_study_prefix
        )
        
        CREATE_CROSS2_OBJECT(
            GENERATE_CONTROL_FILE.out.control_file,
            GENERATE_CONTROL_FILE.out.founder_genos.mix(
                GENERATE_CONTROL_FILE.out.genetic_maps,
                GENERATE_CONTROL_FILE.out.physical_maps,
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
        
        // MODULE 4: Genome Scanning and Permutation Testing
        PREPARE_GENOME_SCAN(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            ch_study_prefix
        )
        
        GENOME_SCAN(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            PREPARE_GENOME_SCAN.out.genoprob,
            PREPARE_GENOME_SCAN.out.kinship,
            ch_study_prefix,
            Channel.value(params.lod_threshold)
        )
        
        PERMUTATION_TEST(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            PREPARE_GENOME_SCAN.out.genoprob,
            PREPARE_GENOME_SCAN.out.kinship,
            GENOME_SCAN.out.filtered_phenotypes,
            ch_study_prefix,
            Channel.value(params.lod_threshold)
        )
        
        IDENTIFY_SIGNIFICANT_QTLS(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            GENOME_SCAN.out.scan_results,
            PERMUTATION_TEST.out.perm_results,
            PERMUTATION_TEST.out.thresholds,
            ch_study_prefix
        )

        // MODULE 5: QTL Viewer Integration
        PREPARE_QTLVIEWER_DATA(
            CREATE_CROSS2_OBJECT.out.cross2_object,
            PREPARE_GENOME_SCAN.out.genoprob,
            PREPARE_GENOME_SCAN.out.kinship,
            PREPARE_GENOME_SCAN.out.genetic_map,
            GENOME_SCAN.out.scan_results,
            IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls,
            ch_study_prefix
        )

        SETUP_QTLVIEWER_DEPLOYMENT(
            PREPARE_QTLVIEWER_DATA.out.qtlviewer_data,
            ch_study_prefix
        )

        // Display results for Module 4
        PREPARE_GENOME_SCAN.out.validation_report.view { "Genome scan preparation report: $it" }
        GENOME_SCAN.out.validation_report.view { "Genome scan report: $it" }
        GENOME_SCAN.out.peaks.view { "Preliminary peaks found: $it" }
        PERMUTATION_TEST.out.validation_report.view { "Permutation test report: $it" }
        PERMUTATION_TEST.out.thresholds.view { "Significance thresholds: $it" }
        IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls.view { "Significant QTLs identified: $it" }
        IDENTIFY_SIGNIFICANT_QTLS.out.qtl_summary.view { "QTL summary report: $it" }

        // Display results for Module 5
        PREPARE_QTLVIEWER_DATA.out.validation_report.view { "QTL Viewer conversion report: $it" }
        SETUP_QTLVIEWER_DEPLOYMENT.out.docker_compose.view { "Docker Compose config created: $it" }
        SETUP_QTLVIEWER_DEPLOYMENT.out.startup_script.view { "QTL Viewer startup script: $it" }
        SETUP_QTLVIEWER_DEPLOYMENT.out.instructions.view { "QTL Viewer instructions: $it" }
        
    } else {
        log.info "Skipping genotype processing - no FinalReport file specified"
        log.info "To include genotype processing, add: --finalreport_file Data/YourFinalReport.txt"
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
    Results  : ${params.outdir}
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