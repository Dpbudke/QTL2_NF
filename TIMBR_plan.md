# Plan: Module 10 — TIMBR Allelic Series Analysis

**Save this plan to the project CWD when starting implementation:**
```bash
cp /home/dawson.budke/.claude/plans/woolly-jumping-patterson.md \
   /project/do2_projects/DO_Choline/QTL2_NF/TIMBR_plan.md
```

---

## Context

TIMBR (Tree-based Inference of Multiallelism via Bayesian Regression; Crouse et al. 2020, *Genetics* 216:957–983) answers a key question after genome scanning: at a significant QTL, how many distinct functional alleles exist among the 8 DO founder strains? Rather than assuming all 8 founders carry unique alleles, TIMBR partitions founders into functionally equivalent groups (the "allelic series") using a Bayesian MCMC approach on diplotype probabilities at the QTL peak.

The pipeline currently ends at module 09 (visualization). Module 10 will take the significant QTLs from module 08 and run TIMBR on each, producing per-QTL allelic series inference, posterior effect plots, and a master summary CSV.

**TIMBR is not currently in the container — a Dockerfile rebuild is required before execution.**

---

## Critical Files

| File | Action |
|------|--------|
| `modules/10_timbr.nf` | **Create new** — three processes + sub-workflow |
| `main.nf` | Add include + opt-in workflow block after VISUALIZE_QTLS |
| `nextflow.config` | Add 4 params + 3 withName resource blocks |
| `Dockerfile` | Add `ape` + TIMBR GitHub install |

**Reference/input files (read-only):**
- `modules/08_identify_significant_qtls.nf` — significant_qtls.csv column names
- `modules/09_visualize.nf` — module structure template
- `Results_final/eQTL/Diet_additive/08_significant_qtls/{prefix}_significant_qtls.csv`
- `Results_final/eQTL/Diet_additive/05_genome_scan_preparation/{prefix}_genoprob.rds`

---

## Architecture: Three-Process + Sub-workflow

```
TIMBR_SETUP  →  TIMBR_BATCH (×N parallel SLURM jobs)  →  TIMBR_AGGREGATE
```

This mirrors the Module 06 and 07 batch patterns.

---

## Step 1 — Dockerfile: Add TIMBR

In `Dockerfile`, after the Bioconductor install block, add:

```dockerfile
# Module 10: TIMBR allelic series analysis
RUN R -e "install.packages('ape', repos='https://cran.rstudio.com/', dependencies=TRUE)"
RUN R -e "devtools::install_github('wesleycrouse/TIMBR', dependencies=TRUE, upgrade='never')"
RUN R -e "library(TIMBR); cat('TIMBR version:', as.character(packageVersion('TIMBR')), '\n')"
```

After adding, rebuild and push to Docker Hub, then re-pull to Singularity cache:
```bash
docker build -t dpbudke/qtl2-pipeline:latest .
docker push dpbudke/qtl2-pipeline:latest
# Then on HPC:
singularity pull --force singularity_cache/dpbudke-qtl2-pipeline-latest.img \
  docker://dpbudke/qtl2-pipeline:latest
```

---

## Step 2 — nextflow.config: Add Params and Resources

In the `params {}` block (after `run_perm_benchmark`):

```groovy
// Module 10: TIMBR Allelic Series Analysis
run_timbr            = false  // opt-in: requires rebuilt container
timbr_sig_level      = "95%"  // filter: "63%", "90%", "95%", or "99%"
timbr_qtls_per_batch = 50     // QTLs per SLURM batch job
timbr_samples        = 10000  // MCMC samples per TIMBR() call
```

In the `process {}` block (after VISUALIZE_QTLS entry):

```groovy
withName: TIMBR_SETUP {
    cpus = 4; memory = '32 GB'; time = '30m'
    clusterOptions = '--account=do2_projects --partition=ceres'
}
withName: TIMBR_BATCH {
    cpus = 4; memory = '32 GB'; time = '4h'  // 50 QTLs × ~2 min + 10-12 GB genoprob
    clusterOptions = '--account=do2_projects --partition=ceres'
    errorStrategy = 'retry'; maxRetries = 1
}
withName: TIMBR_AGGREGATE {
    cpus = 4; memory = '16 GB'; time = '30m'
    clusterOptions = '--account=do2_projects --partition=ceres'
}
```

---

## Step 3 — main.nf: Include + Opt-in Block

Add include (after module 09 include):
```nextflow
include { TIMBR_ANALYSIS } from './modules/10_timbr.nf'
```

Add workflow block (after VISUALIZE_QTLS, still inside the `if (params.finalreport_files)` block):
```nextflow
// MODULE 10: TIMBR Allelic Series Analysis (opt-in)
if (params.run_timbr) {
    TIMBR_ANALYSIS(
        CHUNKED_PERMUTATION_TESTING.out.filtered_cross2,
        PREPARE_GENOME_SCAN_SETUP.out.genoprob,
        PREPARE_GENOME_SCAN_SETUP.out.genetic_map,
        IDENTIFY_SIGNIFICANT_QTLS.out.significant_qtls,
        ch_study_prefix,
        Channel.value(params.timbr_sig_level),
        Channel.value(params.timbr_qtls_per_batch),
        Channel.value(params.timbr_samples)
    )
} else {
    log.info "Skipping TIMBR (--run_timbr not set)"
}
```

---

## Step 4 — modules/10_timbr.nf: Full Implementation

### 4a. TIMBR_SETUP Process

**Inputs:** `cross2_file`, `significant_qtls_file`, `val(prefix)`, `val(sig_level)`, `val(qtls_per_batch)`
**Outputs:** `filtered_qtls`, `batch_file`, `addcovar_rds`, `genetic_map_file` (passed through), `setup_log`, `stdout` (batch IDs)

**R logic:**
1. Load cross2 and significant_qtls CSV
2. Filter to QTLs at or above `sig_level` (order: 63% < 90% < 95% < 99%); deduplicate by `(lodcolumn, chr, pos)` keeping highest significance row
3. Create batch assignments: `batch_id = ceiling(row_number / qtls_per_batch)`; write to tab-delimited file
4. Pre-compute addcovar matrix using same logic as modules 06/07 (strip `coat_color`, `model.matrix()`); save as RDS
5. Write batch count to stdout (one line per batch ID) for channel splitting
6. Write setup_log with QTL counts, batch counts, filter applied

### 4b. TIMBR_BATCH Process

**Inputs:** `each batch_id`, `genoprob_file` (8.5 GB), `genetic_map_file`, `filtered_qtls_file`, `batch_file`, `cross2_file`, `addcovar_rds`, `val(prefix)`, `val(timbr_samples)`
**Outputs:** `batch_summary` (CSV), `batch_log` (TXT), `timbr_plots` (optional PNGs), `timbr_rds` (optional RDS)

**R logic per QTL (loop over `this_batch_rows`):**
```r
# 1. Find nearest pseudomarker to peak position
gmap_chr    <- genetic_map[[chr]]          # from genetic_map.rds (pseudomarker grid)
nearest_mk  <- names(which.min(abs(gmap_chr - pos_cM)))

# 2. Extract P matrix [individuals × 36 diplotype states] at that marker
P <- genoprob[[chr]][, , nearest_mk]

# 3. Build prior.D
prior.D <- list(
    P           = P[common_samples, ],
    A           = TIMBR::additive.design(8, "qtl2"),
    fixed.diplo = FALSE
)

# 4. CRP prior.M (computed once before loop)
hp <- TIMBR::calc.concentration.prior(8, 0.05, 0.01)
prior.M <- list(model.type="crp", prior.alpha.type="gamma",
                prior.alpha.shape=hp[1], prior.alpha.rate=hp[2])

# 5. Align y, prior.D$P, and Z to common samples
y <- cross2$pheno[common_samples, pheno_id]
Z <- addcovar[common_samples, ]   # NULL if no covariates

# 6. Run TIMBR
result <- TIMBR::TIMBR(y, prior.D, prior.M, Z=Z, samples=timbr_samples, calc.lnBF=TRUE)

# 7. Save RDS + PNG
saveRDS(result, "timbr_rds_batch{N}/{safe_id}_chr{C}_{pos}cM.rds")
png("timbr_plots_batch{N}/{safe_id}_chr{C}_{pos}cM.png", width=800, height=600)
TIMBR::TIMBR.plot.haplotypes(result); dev.off()

# 8. Summary row: ln_BF, n_functional_alleles (from MAP allelic series),
#    top_series_pattern (e.g. "1,1,2,1,1,3,1,2"), top_series_prob, status
```

**Key notes:**
- Set `set.seed(batch_id * 12345)` at top for reproducibility
- Use `tryCatch` per QTL so a single failure doesn't abort the batch; write `status = "error: ..."` to summary row
- `n_functional_alleles` = number of unique integers in MAP allelic series pattern
- Use `genetic_map.rds` (pseudomarker grid from module 05) NOT `cross2$gmap` (original markers only) — genoprob was computed on the pseudomarker grid

### 4c. TIMBR_AGGREGATE Process

**Inputs:** `path(batch_summaries)` (collected), `path(batch_logs)` (collected), `val(prefix)`
**Outputs:** `master_summary` CSV, `run_report` TXT

**R logic:**
1. `rbind` all `*_batch*_summary.csv` files (sort by batch number before binding)
2. Sort master by chr (1-19, X), then pos
3. Write `{prefix}_timbr_master_summary.csv`
4. Write run report: total QTLs attempted, success/fail counts, distribution of `n_functional_alleles` (how many biallelic vs. 3-allelic, etc.), summary by significance level

### 4d. Sub-workflow TIMBR_ANALYSIS

```nextflow
workflow TIMBR_ANALYSIS {
    take:
    cross2_ch; genoprob_ch; genetic_map_ch
    significant_qtls_ch; prefix_ch
    sig_level_ch; qtls_per_batch_ch; timbr_samples_ch

    main:
    TIMBR_SETUP(cross2_ch, significant_qtls_ch, prefix_ch, sig_level_ch, qtls_per_batch_ch)

    TIMBR_SETUP.out.batch_ids
        .splitText().map { it.trim() }.filter { it != "" }
        .set { ch_batch_ids }

    TIMBR_BATCH(
        ch_batch_ids,
        genoprob_ch.first(),
        genetic_map_ch.first(),
        TIMBR_SETUP.out.filtered_qtls.first(),
        TIMBR_SETUP.out.batch_file.first(),
        cross2_ch.first(),
        TIMBR_SETUP.out.addcovar_rds.first(),
        prefix_ch,
        timbr_samples_ch
    )

    TIMBR_AGGREGATE(
        TIMBR_BATCH.out.batch_summary.collect(),
        TIMBR_BATCH.out.batch_log.collect(),
        prefix_ch
    )

    emit:
    master_summary = TIMBR_AGGREGATE.out.master_summary
    run_report     = TIMBR_AGGREGATE.out.run_report
    timbr_plots    = TIMBR_BATCH.out.timbr_plots
    setup_log      = TIMBR_SETUP.out.setup_log
}
```

---

## Output Layout

```
{outdir}/10_timbr/
├── {prefix}_timbr_qtls.csv              ← filtered QTL list (SETUP)
├── {prefix}_timbr_setup_log.txt
├── {prefix}_timbr_master_summary.csv    ← PRIMARY RESULT (AGGREGATE)
├── {prefix}_timbr_run_report.txt
└── batches/
    ├── {prefix}_timbr_batch1_summary.csv
    ├── {prefix}_timbr_batch1_log.txt
    ├── timbr_plots_batch1/
    │   └── {pheno_id}_chr{C}_{pos}cM.png
    ├── timbr_rds_batch1/
    │   └── {pheno_id}_chr{C}_{pos}cM.rds
    └── ... (repeated per batch)
```

**Master summary CSV columns:**
```
prefix, lodcolumn, chr, pos, lod, significance_level, nearest_marker,
ln_BF, n_functional_alleles, top_series_pattern, top_series_prob,
n_samples_used, status
```

---

## Parallelization Estimates

| Study | QTLs at 95% | Batches (50/batch) | Est. wall time per batch |
|-------|------------|-------------------|--------------------------|
| eQTL (10,350 pheno) | ~12,000 | ~240 SLURM jobs | ~2 hr (50 × ~2 min) |
| LCMS (10 pheno) | ~5–20 | 1 SLURM job | ~10 min |

Memory per TIMBR_BATCH job: ~12 GB (genoprob) + ~2 GB (cross2) + overhead = 32 GB comfortable.

---

## Verification

1. First, test with LCMS (small: ~10 phenotypes, likely ≤20 significant QTLs):
   ```bash
   nextflow run main.nf \
     --study_type LCMS \
     --analysis_type Diet_additive \
     --phenotype_file Data/QTL2_NF_meta_pheno_LCMS_input.csv \
     --study_prefix DOchln \
     --run_timbr \
     --timbr_sig_level "95%" \
     -profile standard -bg
   ```
2. Confirm `Results_final/LCMS/Diet_additive/10_timbr/` contains master_summary.csv
3. Spot-check 1–2 RDS files in R:
   ```r
   r <- readRDS("timbr_rds_batch1/TMAO_..._chr5_....rds")
   r$ln.BF          # should be a number
   head(r$p.M.given.y, 5)   # allelic series posterior
   ```
4. If successful on LCMS, run on eQTL (larger scale test):
   ```bash
   nextflow run main_resume.nf \
     --analysis_type Diet_additive \
     --resume_from timbr \
     --study_prefix DOchln \
     --run_timbr \
     -profile standard -bg
   ```
   *(Note: `main_resume.nf` will need `timbr` added as a valid `--resume_from` option)*

---

## Implementation Order

1. **Dockerfile** — Add TIMBR; rebuild container; update Singularity image
2. **nextflow.config** — Add params + resource blocks
3. **modules/10_timbr.nf** — Create new file (TIMBR_SETUP, TIMBR_BATCH, TIMBR_AGGREGATE, TIMBR_ANALYSIS)
4. **main.nf** — Add include + opt-in TIMBR_ANALYSIS call
5. **Test on LCMS** — Small scale validation
6. **main_resume.nf** — Add `timbr` as valid resume_from option (optional, for future use)
