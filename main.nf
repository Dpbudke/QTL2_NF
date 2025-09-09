#!/usr/bin/env nextflow

// QTL2_NF: DO_pipe converted to Nextflow
// A pipeline for QTL analysis using r/qtl2

nextflow.enable.dsl = 2

// Include modules
include { PHENOTYPE_PROCESS } from './modules/phenotype_process'

// Main workflow
workflow {
    // Input parameters validation
    if (!params.phenotype_file) {
        error "Please provide a phenotype file with --phenotype_file"
    }
    if (!params.prefix) {
        error "Please provide a study prefix with --prefix"
    }

    // Create input channels
    phenotype_ch = Channel.fromPath(params.phenotype_file)
    
    // Run phenotype processing
    PHENOTYPE_PROCESS(phenotype_ch, params.prefix)
    
    // Emit results for potential downstream processes
    PHENOTYPE_PROCESS.out.covar.view { "Covariate file: $it" }
    PHENOTYPE_PROCESS.out.pheno.view { "Phenotype file: $it" }
}

workflow.onComplete {
    println "Workflow completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
