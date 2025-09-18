#!/usr/bin/env nextflow

// DO Sample Mixup QC Pipeline
// Standalone QC tool for detecting sample mix-ups in Diversity Outbred mouse data
// Based on methodology from Broman et al. (2015) G3

nextflow.enable.dsl = 2

// Import the mixup QC module
include { DO_MIXUP_QC } from './modules/do_mixup_qc'

// Print banner
log.info """
========================================
 DO Sample Mixup QC Pipeline
========================================
phenotype_file : ${params.phenotype_file}
cross2_object  : ${params.cross2_object}
qtl_results    : ${params.qtl_results ?: 'None (will perform cis-eQTL analysis)'}
study_prefix   : ${params.study_prefix}
outdir         : ${params.outdir}
========================================
"""

workflow {
    // Input validation
    if (!params.phenotype_file) {
        error "Error: --phenotype_file parameter is required"
    }
    if (!params.cross2_object) {
        error "Error: --cross2_object parameter is required"
    }

    // Prepare input files
    phenotype_file = file(params.phenotype_file)
    cross2_object = file(params.cross2_object)

    // QTL results are optional - if not provided, will perform simplified eQTL analysis
    qtl_results = params.qtl_results ? file(params.qtl_results) : file("NO_FILE")

    // Run DO mixup QC analysis
    DO_MIXUP_QC(
        phenotype_file,
        cross2_object,
        qtl_results,
        params.study_prefix
    )

    // Emit results for potential downstream use
    emit:
    mixup_report = DO_MIXUP_QC.out.mixup_report
    problems = DO_MIXUP_QC.out.mixup_problems
    distances = DO_MIXUP_QC.out.distances
    summary = DO_MIXUP_QC.out.summary
    corrected_mapping = DO_MIXUP_QC.out.corrected_mapping
}

workflow.onComplete {
    log.info """
    ========================================
     DO Mixup QC Pipeline Completed
    ========================================
    Results: ${params.outdir}/06_mixup_qc/

    Key outputs:
    - ${params.study_prefix}_mixup_report.html
    - ${params.study_prefix}_mixup_problems.csv
    - ${params.study_prefix}_corrected_sample_mapping.csv
    ========================================
    """
}