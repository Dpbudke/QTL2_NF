# Analysis Modules

This directory contains specialized analysis modules for optional downstream QTL analyses. These are standalone analyses that run after the main pipeline completes.

## Available Modules

| Module | Purpose |
|--------|---------|
| `classify_cis_trans_eqtls.nf` | Classify eQTLs as cis-acting or trans-acting |
| `plot_eqtl_map.nf` | Generate genome-wide eQTL map and cis-resolution plot |
| `broman_mixup_qc.nf` | Detect sample mix-ups using expression vs genotypes |
| `heritability_effect_size.nf` | Compute per-eQTL narrow-sense heritability and peak variance explained |
| `summary_eqtl_plots.nf` | Cross-model eQTL yield summary figures |
| `founder_allele_contributions.nf` | Compute per-eQTL qtl2 BLUP + TIMBR founder effects and generate 7 summary figures |

---

## Broman Sample Mixup QC

### Overview
Detects sample mix-ups in Diversity Outbred mouse data by comparing observed gene expression with genotype-predicted expression. Based on methodology from [Broman et al. (2015) G3](https://kbroman.org/qtl2/assets/vignettes/do_mixups.html) and Westra et al. (2011).

### How It Works

1. **Select Top eQTLs**: Uses the N highest-LOD eQTLs from your significant QTL results
2. **Calculate Predicted Expression**: Uses genotype probabilities at each eQTL to predict expression levels
3. **Compare Distances**: Calculates RMS distance between observed and predicted expression
4. **Identify Mismatches**: Flags samples where another sample's genotypes better predict the expression

### Usage

```groovy
include { BROMAN_MIXUP_QC } from './modules/analyses/broman_mixup_qc.nf'

workflow {
    BROMAN_MIXUP_QC(
        cross2_file,           // Cross2 RDS object
        genoprob_file,         // Genotype probabilities RDS
        expr_data_file,        // Expression matrix (RDS or CSV)
        expr_peaks_file,       // Significant eQTLs CSV (from Module 8)
        study_prefix,          // Study identifier
        100                    // Number of top eQTLs to use
    )
}
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `cross2_file` | Path | Cross2 object RDS file |
| `genoprob_file` | Path | Genotype probabilities RDS file |
| `expr_data_file` | Path | Expression data (samples x genes matrix, RDS or CSV) |
| `expr_peaks_file` | Path | Significant QTLs CSV from Module 8 |
| `study_prefix` | String | Study identifier for output naming |
| `n_top_eqtl` | Integer | Number of top eQTLs to analyze (default: 100) |

### Output Files

| File | Description |
|------|-------------|
| `{prefix}_mixup_report.html` | Summary HTML report |
| `{prefix}_mixup_problems.csv` | Samples with potential mix-ups |
| `{prefix}_distance_matrix.csv` | Full RMS distance matrix |
| `{prefix}_mixup_plots.jpg` | Self vs minimum distance plot |
| `mixup_qc_log.txt` | Analysis log |

### Interpreting Results

- **No problems**: All samples' self-distances are lower than their minimum distance to other samples
- **Potential mix-up**: Red points where `min_distance < self_distance` suggest the sample's expression better matches another sample's genotypes
- Follow up with wet-lab verification for flagged samples

### Sample Mixup Remediation Workflow

If flagged samples need to be removed and the analysis rerun, follow these steps:

1. **Review** `{prefix}_mixup_problems.csv` — samples with `is_problem == TRUE` are candidates for removal. Investigate the `best_match` column to confirm likely mix-ups with wet-lab records before removing.

2. **Return to** `Metadata_phenotype_QTL2_NF.Rmd` — set `mixup_problems_file` to the path of your `*_mixup_problems.csv` output (e.g., `"results/mixup_qc/DOchln_mixup_problems.csv"`).

3. **Reknit the Rmd** — this regenerates `Data/QTL2_NF_meta_pheno_input.csv` with the flagged samples excluded.

4. **Rerun the pipeline from Module 1** using the updated input CSV:
   ```bash
   nextflow run main.nf --phenotype_file Data/QTL2_NF_meta_pheno_input.csv --prefix DOchln
   ```

> **Note**: Broman mixup QC requires eQTL results (Module 8 output), so the full pipeline must complete successfully before running mixup QC. The remediation loop is: full pipeline run → mixup QC → filter Rmd → regenerate CSV → rerun pipeline.

---

## Cis vs Trans eQTL Classification

### Overview
Classifies significant eQTLs as **cis-acting** (local) or **trans-acting** (distal) based on the genomic distance between the QTL peak position and the gene's transcription start site (TSS).

### Classification Logic

1. **Extract Gene Information from GTF**:
   - Parse GTF annotation to extract gene positions and strand information
   - Determine TSS based on strand:
     - Forward strand (+): TSS = gene start position
     - Reverse strand (-): TSS = gene end position

2. **Calculate Distance**:
   - Match QTL gene IDs with GTF gene IDs
   - Calculate absolute distance between QTL peak position and gene TSS
   - Distance is measured in Megabases (Mb)

3. **Classification Rules**:
   - **cis**: QTL and gene on same chromosome AND distance ≤ cis_window_mb from TSS
     - Default: 4 Mb on each side = 8 Mb total window (consistent with Keller et al. 2018)
   - **trans**: QTL and gene on different chromosomes OR distance > cis_window_mb
   - **unknown_gene**: Gene ID not found in GTF annotation

### Usage

#### Standalone Workflow
```bash
nextflow run workflow_cis_trans_classification.nf \
    --qtl_file Results_final/results/08_significant_qtls/DOchln_significant_qtls.csv \
    --gtf_file Data/Mus_musculus.GRCm39.113.gtf.gz \
    --cis_window_mb 4.0 \
    --outdir results_cis_trans \
    -profile standard
```

#### As a Module (in another workflow)
```groovy
include { CLASSIFY_CIS_TRANS_EQTLS } from './modules/analyses/classify_cis_trans_eqtls.nf'

workflow {
    qtl_ch = Channel.fromPath("path/to/qtls.csv")
    gtf_ch = Channel.fromPath("path/to/annotation.gtf.gz")

    CLASSIFY_CIS_TRANS_EQTLS(qtl_ch, gtf_ch)
}
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `qtl_file` | Path | Required | Path to significant QTLs CSV file |
| `gtf_file` | Path | Required | Path to GTF annotation file (can be gzipped) |
| `cis_window_mb` | Float | 4.0 | Distance threshold in Mb for cis classification (creates 8 Mb total window: +/- 4 Mb from TSS; consistent with Keller et al. 2018) |
| `outdir` | Path | results_cis_trans | Output directory |

### Input Files

#### QTL File Format
CSV file with columns:
- `lodcolumn`: Gene ID (e.g., ENSMUSG00000000001)
- `chr`: Chromosome of QTL peak
- `pos`: Position of QTL peak in Mb
- `lod`: LOD score
- Additional QTL metadata

#### GTF Annotation File
Standard GTF format (gzipped or uncompressed):
- Gene-level entries are used
- Must contain gene_id in attributes field
- Requires chromosome, start, end, and strand information

### Output Files

1. **eqtl_cis_trans_classification.csv**
   - Full classification results for all QTLs
   - Columns include:
     - QTL information (chr, position, LOD score)
     - Gene information (chr, start, end, strand, TSS)
     - Distance to TSS (in Mb)
     - Classification (cis/trans/unknown_gene)

2. **eqtl_classification_summary.txt**
   - Summary statistics:
     - Count and percentage of cis vs trans eQTLs
     - Mean LOD scores by classification
     - Median distances for each category
     - Breakdown by significance level

### Example Output

```
===========================================
eQTL CIS vs TRANS CLASSIFICATION SUMMARY
===========================================

Analysis Parameters:
  Cis window: +/- 4 Mb from TSS (8 Mb total window)
  Total QTLs analyzed: 1245

Classification Results:
  eqtl_type  count  mean_lod  median_distance_mb  percentage
  cis          342    15.2            0.45            27.5%
  trans        898    12.8           35.20            72.1%
  unknown        5    10.1              NA             0.4%
```

### Biological Interpretation

- **Cis-eQTLs**: Likely affect gene expression through local regulatory mechanisms (e.g., promoter variants, enhancers)
- **Trans-eQTLs**: Suggest regulation through distant factors (e.g., transcription factors, chromatin remodeling)
- Higher proportion of trans-eQTLs is common in complex traits
- Cis-eQTLs typically have higher LOD scores due to stronger local effects

### References

- Classic definition: cis-eQTLs within 1-5 Mb of gene TSS
- Pipeline default: 4 Mb (±4 Mb = 8 Mb total window), consistent with Keller et al. 2018; adjustable via `cis_window_mb` parameter

---

## Genome-wide eQTL Map (`plot_eqtl_map.nf`)

### Overview
Generates two per-model publication-quality figures from the cis/trans classification output: (1) a genome-wide eQTL map (peak position on x-axis, gene TSS position on y-axis, cis points in blue on the diagonal, trans points in black off-diagonal) and (2) a cis-eQTL resolution plot showing signed peak-to-TSS distance versus LOD score for every same-chromosome eQTL.

This module is called automatically by `main_resume.nf` and `main.nf` immediately after `CLASSIFY_CIS_TRANS_EQTLS` whenever `params.study_type == 'eQTL'`. It can also be regenerated in isolation with `--resume_from eqtl_map`.

### Process: `PLOT_EQTL_MAP`

**Resources**: 2 CPUs, 16 GB RAM, 30 minutes

#### Inputs

| Input | Type | Description |
|-------|------|-------------|
| `classification_csv` | Path | `eqtl_cis_trans_classification.csv` from `CLASSIFY_CIS_TRANS_EQTLS` |
| `position_map_rds` | Path | `qtl2_position_map.rds` (gmap/pmap list) from `CLASSIFY_CIS_TRANS_EQTLS` |
| `study_prefix` | String | Study identifier for output naming and plot titles |

#### Output Files

Published to `{params.outdir}/00_analyses/`:

| File | Description |
|------|-------------|
| `{prefix}_eqtl_map.png` | Genome-wide eQTL map (10 x 10 in, 300 DPI); cis = blue, trans = black; alternating chromosome shading |
| `{prefix}_cis_resolution.png` | Signed peak-to-TSS distance (Mb) vs LOD for same-chromosome eQTLs; +/-4 Mb cis threshold marked; pseudo-log x-axis |
| `{prefix}_eqtl_map_log.txt` | Processing log with plotted cis/trans counts and mapping resolution statistics |

#### Notes

- Both plots filter to `significance_level == "95%"` (the distinct, complete 95%-significant set; eQTLs passing the 99% threshold also carry a "95%" row, so this filter captures all eQTLs significant at 95% GW or stronger)
- The cis-resolution plot uses a pseudo-log (symlog) x-axis (`sigma = 2`) so the dense cis cluster near zero and the rare same-chromosome trans eQTLs at chromosome scale are both legible
- The `analysis_type` label in plot titles comes from `params.analysis_type` if set, falling back to `params.study_type` (default `eQTL`)

---

## Heritability and Effect Size (`heritability_effect_size.nf`)

### Overview
Defines two processes that together compute and visualize per-eQTL narrow-sense heritability (h2) and the phenotypic variance explained by the peak marker (Haley-Knott dR2) for all 95% GW eQTLs from the additive model. This module is included by the standalone entry-point launcher `heritability_eqtl.nf` at the repository root.

### Process: `COMPUTE_HERIT_VAREXP`

Computes two statistics per eQTL (95% GW, cis and trans) from the additive model's scan-prep outputs:

- **Narrow-sense heritability (h2)**: `qtl2::est_herit()` with a single overall kinship matrix; one h2 value per unique gene, shared across all eQTLs for that gene.
- **Variance explained (dR2)**: Haley-Knott dR2 = R2 of the full model (covariates + 7 founder allele probabilities at peak marker) minus R2 of the null model (covariates only); one value per eQTL peak.

**Resources**: 32 CPUs, 64 GB RAM, 2 hours

#### Inputs

| Input | Type | Description |
|-------|------|-------------|
| `study_prefix` | String | Study identifier |
| `cross2_rds` | Path | Cross2 object (`04_cross2_creation/{prefix}_cross2.rds`) |
| `alleleprob_rds` | Path | Allele probabilities (`05_genome_scan_preparation/{prefix}_alleleprob.rds`) |
| `genetic_map_rds` | Path | Genetic map at pseudomarker resolution (`05_genome_scan_preparation/{prefix}_genetic_map.rds`) |
| `classification_csv` | Path | `eqtl_cis_trans_classification.csv` from the additive model's `00_analyses/` |

#### Output Files

| File | Description |
|------|-------------|
| `{prefix}_herit_varexp.csv` | Per-eQTL table: `lodcolumn, gene_name, eqtl_type, qtl_chr, qtl_pos_cM, qtl_pos_mb, lod, h2, varexp` |
| `{prefix}_herit_compute_log.txt` | Log with median h2 and dR2 values and per-step timing |

---

### Process: `PLOT_HERIT_EFFECTSIZE`

Reads `{prefix}_herit_varexp.csv` and the classification CSVs from the three diet models to assign each gene a diet-dependence category and produce four figures.

**Resources**: 2 CPUs, 16 GB RAM, 30 minutes

#### Inputs

| Input | Type | Description |
|-------|------|-------------|
| `study_prefix` | String | Study identifier |
| `herit_csv` | Path | `{prefix}_herit_varexp.csv` from `COMPUTE_HERIT_VAREXP` |
| `results_base` | String | Absolute path to the per-model results base (e.g. `Results_final/eQTL`) |
| `diet_models_str` | String | Comma-separated model names for diet categories (e.g. `AIN76_only,HC_only,Diet_interactive`) |

#### Diet-Dependence Categories

Each gene is assigned one mutually exclusive category based on 95% GW eQTL detection in the single-diet and interactive models:

| Category | Detection criterion |
|----------|---------------------|
| `constitutive` | 95% eQTL in both AIN76_only AND HC_only |
| `AIN-specific` | 95% eQTL in AIN76_only but NOT HC_only |
| `HC-specific` | 95% eQTL in HC_only but NOT AIN76_only |
| `interactive` | 95% eQTL in Diet_interactive (and not already constitutive/specific) |
| `uncategorized` | Does not meet any of the above criteria |

#### Output Files

Published to `{params.summary_outdir}/Summary/Heritability/`:

| File | Description |
|------|-------------|
| `{prefix}_heritability_distribution.png` | Histogram of h2 across unique eQTL genes; median marked with red dashed line |
| `{prefix}_herit_vs_varexp.png` | Scatter: h2 vs dR2 per eQTL, colored by LOD (log10 scale); Spearman rho annotated |
| `{prefix}_varexp_cis_trans.png` | Boxplot comparing dR2 for cis vs trans eQTLs; Wilcoxon rank-sum p-value annotated |
| `{prefix}_varexp_by_diet_category.png` | Boxplot of dR2 stratified by diet-dependence category; Kruskal-Wallis p-value annotated |
| `{prefix}_diet_categories.csv` | `herit_varexp.csv` rows with an added `category` column |
| `{prefix}_herit_summary.txt` | Plain-text digest of all summary statistics behind Figures A–D: h2 distribution (mean/median/IQR/range), Spearman rho + p (h2 vs dR2), cis vs trans n/medians/means + Wilcoxon, per-category n/medians/means + Kruskal-Wallis, and overall category counts |
| `{prefix}_herit_figures_log.txt` | Log with per-category counts and cis/trans median dR2 values |

---

## Cross-model eQTL Summary (`summary_eqtl_plots.nf`)

### Overview
Defines a single process that reads the cis/trans classification CSV from each analysis model and produces a stacked bar chart comparing eQTL yield (distinct 95% GW cis and trans counts) across all models. This module is included by the standalone entry-point launcher `summary_eqtl.nf` at the repository root.

### Process: `SUMMARY_EQTL_PLOTS`

**Resources**: 2 CPUs, 16 GB RAM, 30 minutes

#### Inputs

| Input | Type | Description |
|-------|------|-------------|
| `study_prefix` | String | Study identifier for output naming and plot title |
| `results_base` | String | Absolute path to the per-model results base (e.g. `Results_final/eQTL`) |
| `models_str` | String | Comma-separated model names in plot order (e.g. `AIN76_only,HC_only,Diet_additive,Diet_interactive`) |

#### Count Logic

Reads `<results_base>/<model>/00_analyses/eqtl_cis_trans_classification.csv` for each model and filters to `significance_level == "95%"` with `eqtl_type %in% c("cis","trans")`. Because the classification table nests significance levels (a QTL passing 95% carries rows at 63%, 90%, 95%, and 99% if it passes that threshold too), `significance_level == "95%"` selects the distinct, complete set of 95%-significant eQTLs.

#### Output Files

Published to `{params.summary_outdir}/Summary/`:

| File | Description |
|------|-------------|
| `{prefix}_eqtl_cis_trans_barchart.png` | Stacked bar chart: cis (blue) and trans (black) eQTL counts per model, with per-segment and total count labels |
| `{prefix}_eqtl_counts_by_model.csv` | Wide-format table with columns `model`, `cis`, `trans` |
| `{prefix}_eqtl_summary_log.txt` | Processing log with per-model eQTL counts |

---

## Founder Allele Contributions to eQTL (`founder_allele_contributions.nf`)

### Overview
Defines two processes that together compute and visualize the per-founder contribution to every 99% GW additive-model eQTL. For each eQTL, `COMPUTE_FOUNDER_EFFECTS` derives qtl2 BLUP founder effects at the peak marker (via single-marker `scan1blup` using the precomputed LOCO kinship) and reads each per-locus TIMBR RDS to extract allele-collapsed posterior effects. `PLOT_FOUNDER_CONTRIBUTIONS` reads the resulting CSVs and produces 7 publication-quality figures plus a plain-text summary. This module is included by the standalone entry-point launcher `founder_alleles_eqtl.nf` at the repository root.

**Founder code → strain (qtl2 order)**:

| Code | Strain |
|------|--------|
| A | A/J |
| B | C57BL/6J |
| C | 129S1/SvImJ |
| D | NOD/ShiLtJ |
| E | NZO/HlLtJ |
| F | CAST/EiJ |
| G | PWK/PhJ |
| H | WSB/EiJ |

---

### Process: `COMPUTE_FOUNDER_EFFECTS`

For each 99% GW additive cis/trans eQTL, computes two complementary sets of per-founder effects:

- **qtl2 BLUP effects**: `scan1blup()` at the single nearest-marker position on the chromosome, using the precomputed LOCO kinship matrix. Results are additionally center-scaled per locus (z-score across the 8 founders within each eQTL row) to facilitate cross-locus comparison.
- **TIMBR allele-collapsed effects**: Reads the per-locus TIMBR RDS (from `10_timbr/chr*/...`) and extracts `colMeans(post.hap.effects)`, giving a per-founder posterior mean effect that reflects TIMBR's allele-grouping model.

Diet is treated as binary (ain76 vs hc); samples labeled `ain76a` are recoded to `ain76` before model matrix construction, mirroring the covariate encoding used in Module 6.

**Resources**: 32 CPUs, 64 GB RAM, 4 hours

#### Inputs

| Input | Type | Description |
|-------|------|-------------|
| `study_prefix` | String | Study identifier |
| `cross2_rds` | Path | Cross2 object (`04_cross2_creation/{prefix}_cross2.rds`) |
| `alleleprob_rds` | Path | Allele probabilities (`05_genome_scan_preparation/{prefix}_alleleprob.rds`) |
| `genetic_map_rds` | Path | Genetic map at pseudomarker resolution (`05_genome_scan_preparation/{prefix}_genetic_map.rds`) |
| `kinship_loco_rds` | Path | LOCO kinship matrices (`05_genome_scan_preparation/{prefix}_kinship_loco.rds`) |
| `classification_csv` | Path | `eqtl_cis_trans_classification.csv` from the additive model's `00_analyses/`; filtered to `significance_level == "99%"` |
| `timbr_master_csv` | Path | `{prefix}_timbr_master_summary.csv` from `10_timbr/` |
| `timbr_dir` | String | Absolute path to `<additive_model>/10_timbr` (for per-locus RDS discovery) |
| `sig_level` | String | Significance threshold string, e.g. `"99%"` |

#### Output Files

| File | Description |
|------|-------------|
| `{prefix}_founder_blup.csv` | Per-eQTL qtl2 BLUP effects: `lodcolumn, gene_name, eqtl_type, chr, pos, peak_Mb, lod, blup_A..blup_H, z_A..z_H` |
| `{prefix}_timbr_effects.csv` | Per-locus TIMBR allele-collapsed effects: `lodcolumn, chr, pos, peak_Mb, lod, n_functional_alleles, top_series_pattern, timbr_A..timbr_H` |
| `{prefix}_founder_effects.csv` | BLUP + TIMBR merged on `lodcolumn + chr`; all columns from both tables |
| `{prefix}_founder_compute_log.txt` | Compute log with per-chromosome locus counts and non-NA effect counts |

---

### Process: `PLOT_FOUNDER_CONTRIBUTIONS`

Reads the three CSVs from `COMPUTE_FOUNDER_EFFECTS` and produces 7 figures and a plain-text digest. Implements a CC/DO consensus founder color palette keyed by founder code. The Dunn post-hoc test is implemented manually (BH correction via `p.adjust`) because `dunn.test`/`FSA` are not available in the container.

**Resources**: 2 CPUs, 16 GB RAM, 30 minutes

#### Inputs

| Input | Type | Description |
|-------|------|-------------|
| `study_prefix` | String | Study identifier |
| `founder_blup_csv` | Path | `{prefix}_founder_blup.csv` from `COMPUTE_FOUNDER_EFFECTS` |
| `timbr_effects_csv` | Path | `{prefix}_timbr_effects.csv` from `COMPUTE_FOUNDER_EFFECTS` |
| `founder_effects_csv` | Path | `{prefix}_founder_effects.csv` from `COMPUTE_FOUNDER_EFFECTS` |

#### Output Files

Published to `{params.summary_outdir}/Summary/FounderAlleles/`:

| File | Figure | Description |
|------|--------|-------------|
| `{prefix}_founder_effect_boxplots.png` | A | Per-founder boxplots of center-scaled BLUP z-scores; founders ordered by decreasing median; CC/DO palette; Kruskal-Wallis p-value annotated; BH-corrected Dunn post-hoc pairs for CAST(F)/PWK(G) reported in the summary text |
| `{prefix}_founder_effects_circos.png` | B | circlize circular genome plot with 8 founder-specific radial tracks; each point is one eQTL peak colored by founder; y-axis is the center-scaled BLUP z-score (1st–99th percentile range) |
| `{prefix}_timbr_singleton_frequency.png` | C1 | Bar chart: per-founder fraction of TIMBR loci at which the founder is its own (singleton) allele group in the MAP allelic series pattern |
| `{prefix}_timbr_higheffect_frequency.png` | C2 | Bar chart: per-founder fraction of loci at which the founder belongs to the allele group with the largest mean absolute TIMBR effect |
| `{prefix}_blup_vs_timbr_topLOD.png` | D1 | Paired-bar comparison of mean-centered qtl2 BLUP vs TIMBR collapsed effects for the up-to-6 highest-LOD eQTL with both BLUP and TIMBR data; each locus is a facet panel |
| `{prefix}_blup_vs_timbr_wildderived.png` | D2 | Paired-bar comparison for up-to-6 eQTL where CAST (F) or PWK (G) is resolved as a singleton allele in the MAP pattern |
| `{prefix}_blup_vs_timbr_spanalleles.png` | D3 | Paired-bar comparison for up-to-8 eQTL, one representative per distinct `n_functional_alleles` value, spanning the range of allelic complexity |
| `{prefix}_founder_summary.txt` | — | Plain-text digest: Fig A per-founder median z-scores + Kruskal-Wallis + Dunn pairs for CAST/PWK; Fig C1/C2 per-founder frequencies; Fig D1–D3 selected gene names |
| `{prefix}_founder_figures_log.txt` | — | Figures log with completion status and key statistics for each figure |
