#!/usr/bin/env nextflow

/*
========================================================================================
    QTL2_NF: Resume from GENOME_SCAN step
========================================================================================
    Resume pipeline execution starting from the genome scan step using existing outputs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Default parameters
params {
    study_prefix       = 'DOchln'
    lod_threshold      = 7.0
    outdir             = 'results'
    help               = false
}

// Print help message
def helpMessage() {
    log.info"""
    ========================================
     QTL2_NF Resume from GENOME_SCAN
    ========================================

    Usage:
        nextflow run resume_from_genome_scan.nf [options]

    Optional Arguments:
        --study_prefix      Study identifier prefix [default: DOchln]
        --lod_threshold     LOD threshold for filtering QTLs [default: 7.0]
        --outdir            Output directory [default: results]

    This script resumes execution from the GENOME_SCAN step using existing:
        - Cross2 object from 04_cross2_creation/
        - Genotype probabilities from 05_genome_scan_preparation/
        - Kinship matrices from 05_genome_scan_preparation/

    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { GENOME_SCAN; PERMUTATION_TEST; IDENTIFY_SIGNIFICANT_QTLS } from './modules/scan_perm.nf'
include { PREPARE_QTLVIEWER_DATA; SETUP_QTLVIEWER_DEPLOYMENT } from './modules/qtl_viewer.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    log.info """
    ========================================
     QTL2_NF Resume from GENOME_SCAN
    ========================================
    study_prefix   : ${params.study_prefix}
    outdir         : ${params.outdir}
    lod_threshold  : ${params.lod_threshold}
    ========================================
    """.stripIndent()

    // Check that required input files exist
    cross2_file = file("${params.outdir}/04_cross2_creation/${params.study_prefix}_cross2.rds")
    genoprob_file = file("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_genoprob.rds")
    kinship_file = file("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_kinship_loco.rds")
    genetic_map_file = file("${params.outdir}/05_genome_scan_preparation/${params.study_prefix}_genetic_map.rds")

    if (!cross2_file.exists()) {
        error "Cross2 file not found: ${cross2_file}. Run the pipeline through CREATE_CROSS2_OBJECT first."
    }
    if (!genoprob_file.exists()) {
        error "Genoprob file not found: ${genoprob_file}. Run the pipeline through PREPARE_GENOME_SCAN first."
    }
    if (!kinship_file.exists()) {
        error "Kinship file not found: ${kinship_file}. Run the pipeline through PREPARE_GENOME_SCAN first."
    }

    log.info "✓ All required input files found"

    // Create input channels
    ch_cross2 = Channel.fromPath(cross2_file)
    ch_genoprob = Channel.fromPath(genoprob_file)
    ch_kinship = Channel.fromPath(kinship_file)
    ch_genetic_map = Channel.fromPath(genetic_map_file)
    ch_study_prefix = Channel.value(params.study_prefix)

    // MODULE 4: Genome Scanning and Permutation Testing
    log.info "Starting GENOME_SCAN with 1TB memory and 8 cores..."

    GENOME_SCAN(
        ch_cross2,
        ch_genoprob,
        ch_kinship,
        ch_study_prefix,
        Channel.value(params.lod_threshold)
    )

    PERMUTATION_TEST(
        ch_cross2,
        ch_genoprob,
        ch_kinship,
        GENOME_SCAN.out.filtered_phenotypes,
        ch_study_prefix,
        Channel.value(params.lod_threshold)
    )

    IDENTIFY_SIGNIFICANT_QTLS(
        ch_cross2,
        GENOME_SCAN.out.scan_results,
        PERMUTATION_TEST.out.perm_results,
        PERMUTATION_TEST.out.thresholds,
        ch_study_prefix
    )

    // MODULE 5: QTL Viewer Integration
    PREPARE_QTLVIEWER_DATA(
        ch_cross2,
        ch_genoprob,
        ch_kinship,
        ch_genetic_map,
        GENOME_SCAN.out.scan_results,
        IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls,
        ch_study_prefix
    )

    SETUP_QTLVIEWER_DEPLOYMENT(
        PREPARE_QTLVIEWER_DATA.out.qtlviewer_data,
        ch_study_prefix
    )

    // Display results
    GENOME_SCAN.out.validation_report.view { "Genome scan report: $it" }
    GENOME_SCAN.out.peaks.view { "Preliminary peaks found: $it" }
    PERMUTATION_TEST.out.validation_report.view { "Permutation test report: $it" }
    PERMUTATION_TEST.out.thresholds.view { "Significance thresholds: $it" }
    IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls.view { "Significant QTLs identified: $it" }
    IDENTIFY_SIGNIFICANT_QTLS.out.qtl_summary.view { "QTL summary report: $it" }

    // Display QTL Viewer results
    PREPARE_QTLVIEWER_DATA.out.validation_report.view { "QTL Viewer conversion report: $it" }
    SETUP_QTLVIEWER_DEPLOYMENT.out.docker_compose.view { "Docker Compose config created: $it" }
    SETUP_QTLVIEWER_DEPLOYMENT.out.startup_script.view { "QTL Viewer startup script: $it" }
    SETUP_QTLVIEWER_DEPLOYMENT.out.instructions.view { "QTL Viewer instructions: $it" }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================
     QTL2_NF Resume Completed
    ========================================
    Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Time     : ${workflow.duration}
    Workdir  : ${workflow.workDir}
    Results  : ${params.outdir}
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