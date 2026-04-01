# Analysis Modules

This directory contains specialized analysis modules for optional downstream QTL analyses. These are standalone analyses that run after the main pipeline completes.

## Available Modules

| Module | Purpose |
|--------|---------|
| `classify_cis_trans_eqtls.nf` | Classify eQTLs as cis-acting or trans-acting |
| `broman_mixup_qc.nf` | Detect sample mix-ups using expression vs genotypes |

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
     - Default: 2 Mb on each side = 4 Mb total window
   - **trans**: QTL and gene on different chromosomes OR distance > cis_window_mb
   - **unknown_gene**: Gene ID not found in GTF annotation

### Usage

#### Standalone Workflow
```bash
nextflow run workflow_cis_trans_classification.nf \
    --qtl_file Results_final/results/08_significant_qtls/DOchln_significant_qtls.csv \
    --gtf_file Data/Mus_musculus.GRCm39.113.gtf.gz \
    --cis_window_mb 2.0 \
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
| `cis_window_mb` | Float | 2.0 | Distance threshold in Mb for cis classification (creates 4 Mb total window: +/- 2 Mb from TSS) |
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
  Cis window: +/- 2 Mb from TSS (4 Mb total window)
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
- Commonly used threshold: 2 Mb (adjustable via `cis_window_mb` parameter)
