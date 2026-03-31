#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    eQTL Cis vs Trans Classification Workflow
========================================================================================
    Classifies significant eQTLs as cis-acting (local) or trans-acting (distal)
    based on distance between QTL peak and gene transcription start site (TSS).

    Usage:
        nextflow run workflow_cis_trans_classification.nf \
            --qtl_file Results_final/results/08_significant_qtls/DOchln_significant_qtls.csv \
            --gtf_file Data/Mus_musculus.GRCm39.105.gtf.gz \
            --outdir results_cis_trans \
            -profile standard

    Parameters:
        --qtl_file        : Path to significant QTLs CSV file
        --gtf_file        : Path to GTF annotation file (gzipped)
        --cis_window_mb   : Distance threshold in Mb for cis classification (default: 2.0, creates 4 Mb total window)
        --outdir          : Output directory (default: results_cis_trans)
========================================================================================
*/

// Parameters
params.qtl_file = null
params.gtf_file = "Data/Mus_musculus.GRCm39.105.gtf.gz"
params.cis_window_mb = 2.0
params.outdir = "results_cis_trans"

// Validate required parameters
if (!params.qtl_file) {
    error "ERROR: --qtl_file parameter is required!"
}

// Print parameters
log.info """
========================================
 eQTL Cis/Trans Classification
========================================
QTL file     : ${params.qtl_file}
GTF file     : ${params.gtf_file}
Cis window   : +/- ${params.cis_window_mb} Mb (${params.cis_window_mb * 2} Mb total)
Output dir   : ${params.outdir}
========================================
"""

// Include the module
include { CLASSIFY_CIS_TRANS_EQTLS } from './modules/analyses/classify_cis_trans_eqtls.nf'

workflow {
    // Create input channels
    qtl_ch = Channel.fromPath(params.qtl_file, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf_file, checkIfExists: true)

    // Run classification
    CLASSIFY_CIS_TRANS_EQTLS(qtl_ch, gtf_ch)
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed!
    ========================================
    Results directory: ${params.outdir}/09_eqtl_classification

    Output files:
      - eqtl_cis_trans_classification.csv  : Full classification results
      - eqtl_classification_summary.txt    : Summary statistics
    ========================================
    """
}
