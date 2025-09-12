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
        --finalreport_file  Path to GeneSeek FinalReport file
        --outdir            Output directory [default: results]
        --auto_prefix_samples    Automatically prefix sample IDs [default: false]
        --test_mode         Process chromosome 19 only (for development) [default: false]
    
    Example:
        nextflow run main.nf \\
            --phenotype_file Data/formatted_phenotypes_for_nextflow.csv \\
            --finalreport_file Data/FinalReport.txt \\
            --study_prefix DOChln \\
            --outdir results
    
    Profiles:
        -profile standard   Local execution with Docker (default)
        -profile test       Run with test data (chr 19 only)
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

// Future module imports (commented for now)
// include { CREATE_CROSS      } from './modules/create_cross.nf'
// include { QTL_SCAN          } from './modules/qtl_scan.nf'

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
    test_mode      : ${params.test_mode ? 'true (chr 19 only)' : 'false'}
    ========================================
    """.stripIndent()

    // Create input channels
    ch_phenotype_file = Channel.fromPath(params.phenotype_file, checkIfExists: true)
    ch_study_prefix   = Channel.value(params.study_prefix)
    
    // MODULE 1: Phenotype Processing
    PHENOTYPE_PROCESS(
        ch_phenotype_file,
        ch_study_prefix
    )
    
    // MODULE 2: Genotype Processing (if FinalReport file provided)
    if (params.finalreport_file) {
        ch_finalreport = Channel.fromPath(params.finalreport_file, checkIfExists: true)
        
        GENOTYPE_PROCESS(
            ch_finalreport,
            ch_study_prefix
        )
        
        // Display results for Module 2
        GENOTYPE_PROCESS.out.geno_files.view { "Genotype files created: $it" }
        GENOTYPE_PROCESS.out.validation_report.view { "Genotype validation report: $it" }
        GENOTYPE_PROCESS.out.allele_codes.view { "Reference allele codes downloaded: $it" }
    } else {
        log.info "Skipping genotype processing - no FinalReport file specified"
        log.info "To include genotype processing, add: --finalreport_file Data/YourFinalReport.txt"
    }
    
    // Display results for Module 1
    PHENOTYPE_PROCESS.out.covar.view { "Covariate file created: $it" }
    PHENOTYPE_PROCESS.out.pheno.view { "Phenotype file created: $it" }
    PHENOTYPE_PROCESS.out.validation_report.view { "Validation report: $it" }
    
    // Future workflow steps (commented for now)
    /*
    // MODULE 3: Create Cross Object (future)
    // CREATE_CROSS(
    //     PHENOTYPE_PROCESS.out.pheno,
    //     PHENOTYPE_PROCESS.out.covar,
    //     GENOTYPE_PROCESS.out.geno_files,
    //     GENOTYPE_PROCESS.out.allele_codes
    // )
    
    // MODULE 4: QTL Scanning (future)
    // QTL_SCAN(CREATE_CROSS.out.cross_object)
    */
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
    Test Mode: ${params.test_mode ? 'Enabled (chr 19 only)' : 'Disabled'}
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