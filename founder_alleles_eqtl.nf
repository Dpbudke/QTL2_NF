#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Founder Allele Contributions to eQTL Workflow
========================================================================================
    Founder-resolution view of the additive-model eQTL catalog (95% GW). Per eQTL:
    qtl2 BLUP founder effects at the peak marker (single-marker scan1blup, LOCO kinship)
    and TIMBR allele-collapsed posterior effects. Produces 7 figures + 3 per-eQTL CSVs:
        A   per-founder center-scaled BLUP effect boxplots (Kruskal-Wallis + Dunn post-hoc)
        B   circular genome plot, 8 founder tracks of center-scaled effects
        C1  TIMBR singleton frequency (founder resolved as its own allele)
        C2  TIMBR high-effect-group frequency
        D1  qtl2 BLUP vs TIMBR collapsed effects - highest-LOD eQTL
        D2  qtl2 BLUP vs TIMBR collapsed effects - CAST/PWK resolved as singleton
        D3  qtl2 BLUP vs TIMBR collapsed effects - spanning allele-group counts

    Run once, after the additive model's scan-prep, classification, and TIMBR are complete.

    Usage:
        nextflow run founder_alleles_eqtl.nf --study_prefix DOchln

    Parameters:
        --study_prefix    : Output prefix (default: DOchln) — also the additive-model file prefix
        --results_base    : Base dir of per-model result dirs (default: Results_final/eQTL)
        --additive_model  : Model providing expression/kinship/peaks/TIMBR (default: Diet_additive)
        --sig_level       : eQTL significance level to analyze (default: 95%) — matches the 95% GW TIMBR rerun
    Outputs land in: <results_base>/Summary/FounderAlleles/
========================================================================================
*/

params.study_prefix   = 'DOchln'
params.results_base    = 'Results_final/eQTL'
params.additive_model  = 'Diet_additive'
params.sig_level       = '95%'
params.summary_outdir  = params.results_base

include { COMPUTE_FOUNDER_EFFECTS; PLOT_FOUNDER_CONTRIBUTIONS } from './modules/analyses/founder_allele_contributions.nf'

// Resolve to absolute (the process CWD is the task work dir, so relative paths would not resolve)
def results_base_abs = file(params.results_base).toAbsolutePath().toString()
def add_dir          = "${results_base_abs}/${params.additive_model}"
def timbr_dir        = "${add_dir}/10_timbr"

// Additive-model compute inputs
def cross2_f   = file("${add_dir}/04_cross2_creation/${params.study_prefix}_cross2.rds")
def aprob_f    = file("${add_dir}/05_genome_scan_preparation/${params.study_prefix}_alleleprob.rds")
def gmap_f     = file("${add_dir}/05_genome_scan_preparation/${params.study_prefix}_genetic_map.rds")
def kloco_f    = file("${add_dir}/05_genome_scan_preparation/${params.study_prefix}_kinship_loco.rds")
def clsf_f     = file("${add_dir}/00_analyses/eqtl_cis_trans_classification.csv")
def tmaster_f  = file("${timbr_dir}/${params.study_prefix}_timbr_master_summary.csv")

// Validate all required inputs up front
def required = [cross2_f, aprob_f, gmap_f, kloco_f, clsf_f, tmaster_f]
def missing = required.findAll { !it.exists() }
if (missing) {
    error "Missing required input(s):\n  " + missing.join("\n  ") +
          "\n(Has the additive model's scan-prep, classification, and TIMBR finished?)"
}

log.info """
========================================
 Founder Allele Contributions to eQTL
========================================
Study prefix   : ${params.study_prefix}
Additive model : ${params.additive_model}
Significance   : ${params.sig_level}
TIMBR dir      : ${timbr_dir}
Output dir     : ${params.summary_outdir}/Summary/FounderAlleles
========================================
"""

workflow {
    pfx = Channel.value(params.study_prefix)

    COMPUTE_FOUNDER_EFFECTS(
        pfx, cross2_f, aprob_f, gmap_f, kloco_f, clsf_f, tmaster_f,
        Channel.value(timbr_dir),
        Channel.value(params.sig_level)
    )

    PLOT_FOUNDER_CONTRIBUTIONS(
        pfx,
        COMPUTE_FOUNDER_EFFECTS.out.blup,
        COMPUTE_FOUNDER_EFFECTS.out.timbr,
        COMPUTE_FOUNDER_EFFECTS.out.merged
    )
    PLOT_FOUNDER_CONTRIBUTIONS.out.fig_a.view { "Founder figures -> $it (A,B,C1; TIMBR per-founder figs held)" }
}

workflow.onComplete {
    log.info "Founder allele contributions analysis complete. Results: ${params.summary_outdir}/Summary/FounderAlleles"
}
