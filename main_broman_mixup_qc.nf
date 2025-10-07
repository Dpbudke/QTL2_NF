#!/usr/bin/env nextflow

// Standalone Broman Mixup QC Pipeline
// Detects sample mix-ups using expression vs genotype concordance
// Based on Broman et al. (2015) G3 methodology

nextflow.enable.dsl=2

// Import the Broman mixup QC module
include { BROMAN_MIXUP_QC } from './modules/broman_mixup_qc.nf'

// Main workflow
workflow {
    // Required input files
    cross2_file = file(params.cross2_file)
    genoprob_file = file(params.genoprob_file)
    expr_data_file = file(params.expr_data_file)
    expr_peaks_file = file(params.expr_peaks_file)

    // Run Broman mixup QC
    BROMAN_MIXUP_QC(
        cross2_file,
        genoprob_file,
        expr_data_file,
        expr_peaks_file,
        params.study_prefix,
        params.n_top_eqtl
    )
}

// Workflow completion handler
workflow.onComplete {
    println """
    ================================================
    Broman Mixup QC Analysis Complete
    ================================================
    Study: ${params.study_prefix}
    Top eQTL analyzed: ${params.n_top_eqtl}
    Output directory: ${params.outdir}/mixup_qc
    ================================================
    """.stripIndent()
}
