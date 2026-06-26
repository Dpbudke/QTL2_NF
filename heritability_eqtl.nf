#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    eQTL Heritability & Effect-Size Workflow
========================================================================================
    Per eQTL (additive model, 95% GW): narrow-sense heritability h2 (est_herit, single
    overall kinship) and phenotypic variance explained by the peak marker (Haley-Knott
    dR2). Produces 4 figures + a per-eQTL results CSV. Figure D stratifies effect size by
    a diet-dependence category (constitutive / AIN-specific / HC-specific / interactive)
    derived from cross-model eQTL detection.

    Run once, after the additive model's scan-prep + classification are complete.

    Usage:
        nextflow run heritability_eqtl.nf --study_prefix DOchln

    Parameters:
        --study_prefix    : Output prefix (default: DOchln) — also the additive-model file prefix
        --results_base    : Base dir of per-model result dirs (default: Results_final/eQTL)
        --additive_model  : Model providing expression/kinship/peaks (default: Diet_additive)
        --diet_models     : Comma-separated models for the diet categories
                            (default: AIN76_only,HC_only,Diet_interactive)
    Outputs land in: <results_base>/Summary/Heritability/
========================================================================================
*/

params.study_prefix   = 'DOchln'
params.results_base    = 'Results_final/eQTL'
params.additive_model  = 'Diet_additive'
params.diet_models     = 'AIN76_only,HC_only,Diet_interactive'
params.summary_outdir  = params.results_base

include { COMPUTE_HERIT_VAREXP; PLOT_HERIT_EFFECTSIZE } from './modules/analyses/heritability_effect_size.nf'

// Resolve to absolute (the process CWD is the task work dir, so relative paths would not resolve)
def results_base_abs = file(params.results_base).toAbsolutePath().toString()
def add_dir          = "${results_base_abs}/${params.additive_model}"

// Additive-model compute inputs
def cross2_f = file("${add_dir}/04_cross2_creation/${params.study_prefix}_cross2.rds")
def aprob_f  = file("${add_dir}/05_genome_scan_preparation/${params.study_prefix}_alleleprob.rds")
def gmap_f   = file("${add_dir}/05_genome_scan_preparation/${params.study_prefix}_genetic_map.rds")
def clsf_f   = file("${add_dir}/00_analyses/eqtl_cis_trans_classification.csv")

// Validate all required inputs up front (additive inputs + every diet model's classification CSV)
def diet_list = params.diet_models.tokenize(',')*.trim()
def required = [cross2_f, aprob_f, gmap_f, clsf_f] +
    diet_list.collect { m -> file("${results_base_abs}/${m}/00_analyses/eqtl_cis_trans_classification.csv") }
def missing = required.findAll { !it.exists() }
if (missing) {
    error "Missing required input(s):\n  " + missing.join("\n  ") +
          "\n(Has the additive model's scan-prep/classification and each diet model finished?)"
}

log.info """
========================================
 eQTL Heritability & Effect Size
========================================
Study prefix   : ${params.study_prefix}
Additive model : ${params.additive_model}
Diet models    : ${diet_list.join(', ')}
Output dir     : ${params.summary_outdir}/Summary/Heritability
========================================
"""

workflow {
    pfx = Channel.value(params.study_prefix)

    COMPUTE_HERIT_VAREXP(pfx, cross2_f, aprob_f, gmap_f, clsf_f)

    PLOT_HERIT_EFFECTSIZE(
        pfx,
        COMPUTE_HERIT_VAREXP.out.results,
        Channel.value(results_base_abs),
        Channel.value(diet_list.join(','))
    )
    PLOT_HERIT_EFFECTSIZE.out.fig_a.view { "Heritability figures -> $it (+3 more)" }
}

workflow.onComplete {
    log.info "Heritability/effect-size analysis complete. Results: ${params.summary_outdir}/Summary/Heritability"
}
