#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Cross-model eQTL Summary Workflow
========================================================================================
    Generates summary figures that span ALL analysis models at once (unlike the
    per-model figures in modules/analyses/plot_eqtl_map.nf). Currently produces a
    stacked bar chart of distinct cis/trans eQTL counts (95% GW threshold) per model,
    so diet- and model-dependent shifts in eQTL yield are directly comparable.

    Run once, after all models have produced their cis/trans classification CSV at
    Results_final/eQTL/<model>/00_analyses/eqtl_cis_trans_classification.csv.

    Usage:
        nextflow run summary_eqtl.nf
        nextflow run summary_eqtl.nf --study_prefix DOchln --results_base Results_final/eQTL

    Parameters:
        --study_prefix   : Output file prefix (default: DOchln)
        --results_base   : Base dir holding the per-model result dirs (default: Results_final/eQTL)
        --models         : Comma-separated model names, in plot order
                           (default: AIN76_only,HC_only,Diet_additive,Diet_interactive)
    Outputs land in: <results_base>/Summary/
========================================================================================
*/

params.study_prefix   = 'DOchln'
params.results_base    = 'Results_final/eQTL'
params.models          = 'AIN76_only,HC_only,Diet_additive,Diet_interactive'
params.summary_outdir  = params.results_base   // publishDir base -> <summary_outdir>/Summary

include { SUMMARY_EQTL_PLOTS } from './modules/analyses/summary_eqtl_plots.nf'

// Resolve results_base to an absolute path (the process runs with its CWD set to the
// task work dir, so a relative base would not resolve there).
def results_base_abs = file(params.results_base).toAbsolutePath().toString()

// Validate each model's classification CSV exists up front, with a clear error
// naming any model whose results aren't ready yet.
def model_list = params.models.tokenize(',')*.trim()
def missing = model_list.findAll { m ->
    !file("${results_base_abs}/${m}/00_analyses/eqtl_cis_trans_classification.csv").exists()
}
if (missing) {
    error "Missing eQTL classification CSV for model(s): ${missing.join(', ')} " +
          "(expected at ${params.results_base}/<model>/00_analyses/eqtl_cis_trans_classification.csv). " +
          "Has each model finished?"
}

log.info """
========================================
 Cross-model eQTL Summary
========================================
Study prefix : ${params.study_prefix}
Results base : ${params.results_base}
Models       : ${model_list.join(', ')}
Output dir   : ${params.summary_outdir}/Summary
========================================
"""

workflow {
    SUMMARY_EQTL_PLOTS(
        Channel.value(params.study_prefix),
        Channel.value(results_base_abs),
        Channel.value(model_list.join(','))
    )
    SUMMARY_EQTL_PLOTS.out.barchart.view { "eQTL yield bar chart: $it" }
}

workflow.onComplete {
    log.info "Cross-model eQTL summary complete. Results: ${params.summary_outdir}/Summary"
}
