#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    eQTL Cis vs Trans Classification Workflow
========================================================================================
    Classifies significant eQTLs as cis-acting (local) or trans-acting (distal)
    based on distance between QTL peak and gene transcription start site (TSS).

    IMPORTANT: QTL positions in r/qtl2 are in centiMorgans (cM), not Megabases (Mb).
    This workflow requires the cross2 object to convert QTL positions from cM to Mb
    using the genetic and physical maps before comparing to gene positions.

    Usage:
        nextflow run workflow_cis_trans_classification.nf \
            --qtl_file Results_final/Diet_additive/08_significant_qtls/DOchln_significant_qtls.csv \
            --cross2_file Results_final/Diet_additive/04_cross2_creation/DOchln_cross2.rds \
            --gtf_file Data/Mus_musculus.GRCm39.113.gtf.gz \
            --outdir results_cis_trans \
            -profile standard

    Parameters:
        --qtl_file        : Path to significant QTLs CSV file
        --cross2_file     : Path to cross2 RDS file (for cM to Mb conversion)
        --gtf_file        : Path to GTF annotation file (gzipped)
        --cis_window_mb   : Distance threshold in Mb for cis classification (default: 2.0, creates 4 Mb total window)
        --outdir          : Output directory (default: results_cis_trans)
========================================================================================
*/

// Parameters
params.qtl_file = null
params.cross2_file = null
params.gtf_file = "Data/Mus_musculus.GRCm39.113.gtf.gz"
params.cis_window_mb = 2.0
params.outdir = "results_cis_trans"

// Validate required parameters
if (!params.qtl_file) {
    error "ERROR: --qtl_file parameter is required!"
}
if (!params.cross2_file) {
    error "ERROR: --cross2_file parameter is required for cM to Mb conversion!"
}

// Print parameters
log.info """
========================================
 eQTL Cis/Trans Classification
========================================
QTL file     : ${params.qtl_file}
Cross2 file  : ${params.cross2_file}
GTF file     : ${params.gtf_file}
Cis window   : +/- ${params.cis_window_mb} Mb (${params.cis_window_mb * 2} Mb total)
Output dir   : ${params.outdir}
========================================
"""

// Include the module
include { CLASSIFY_CIS_TRANS_EQTLS } from './modules/analyses/classify_cis_trans_eqtls.nf'

workflow {
    // Create input channels
    qtl_ch    = Channel.fromPath(params.qtl_file, checkIfExists: true)
    gtf_ch    = Channel.fromPath(params.gtf_file, checkIfExists: true)
    cross2_ch = Channel.fromPath(params.cross2_file, checkIfExists: true)

    CLASSIFY_CIS_TRANS_EQTLS(qtl_ch, gtf_ch, cross2_ch)
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed!
    ========================================
    Results directory: ${params.outdir}/00_analyses

    Output files:
      - eqtl_cis_trans_classification.csv  : Full classification results
      - eqtl_classification_summary.txt    : Summary statistics
    ========================================
    """
}
