# QTL2_NF Pipeline

A Nextflow pipeline for multiparental mouse QTL analysis using r/qtl2, designed for high-throughput quantitative trait loci mapping in Diversity Outbred (DO) and other multiparental populations. This pipeline supports comprehensive multi-phenotype analyses including clinical phenotype QTLs, expression QTLs (eQTLs), metabolite QTLs (mQTLs), and other high-dimensional molecular trait mapping applications.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-brightgreen.svg)](https://www.docker.com/)

## Overview

QTL2_NF is a comprehensive bioinformatics pipeline that processes phenotype and genotype data to identify quantitative trait loci (QTLs) with rigorous statistical significance testing. The pipeline is optimized for large-scale datasets and provides interactive visualization through QTL Viewer integration. Building upon methodologies from the DO_Pipe project (https://github.com/exsquire/DO_pipe), this pipeline implements automated workflows for multiparental population analysis with enhanced scalability and reproducibility.

## Input Data Requirements

### Phenotype File Format
The pipeline requires a specially formatted CSV file with header row indicators. This standardized format enables automated parsing and supports diverse phenotype types from clinical measurements to high-dimensional molecular data:

**File Structure Requirements**:
The pipeline uses a specialized CSV format that facilitates automated parsing of metadata and phenotype information:

- **Row 1**: Column type indicators
  - Column A: empty (reserved for sample IDs)
  - Column B: "covariate" (marks the beginning of covariate columns)
  - Column where phenotypes begin: "phenotype" (marks transition from covariates to phenotype data)
  - All other columns in Row 1: empty (continuation of covariate or phenotype sections)
- **Row 2**: Actual column names
  - Column A: "SampleID" (sample identifier column)
  - Covariate columns: Sex, generation, batch, Diet, coat_color, Age, etc.
  - Phenotype columns: transcript names (for eQTL), metabolite IDs (for mQTL), clinical trait names
- **Row 3 and beyond**: Sample data with numeric phenotype values and categorical covariate values

**Required Columns**:
- **SampleID**: Must match genotype sample IDs exactly (e.g., DOchln_001, DOchln_002)
- **Sex**: "male" or "female" (required for r/qtl2 DO crosses)
- **generation**: Numeric generation number (required for DO population structure)

**Optional Covariates**:
- **batch**: Experimental batch identifier for batch effect correction
- **Diet**: Treatment group (e.g., "hc", "lc" for high/low cholesterol)
- **coat_color**: Coat color measurements for validation analyses
- **Age**: Age at measurement or sacrifice
- Additional experimental factors as needed

**Single-Sex Studies**:
For studies with only one sex (all-female or all-male cohorts), the `sex` column must still be present in the covariate file because r/qtl2 uses it as cross metadata for X chromosome handling — it is not only a statistical covariate. Do not omit this column from the input file. If the `sex` column is absent, use `single_sex = 'female'` (or `'male'`) in `nextflow.config` and the pipeline will add it automatically. See the `--single_sex` parameter documentation in the Configuration section for full details.

**Phenotype Data Types Supported**:
- **Clinical phenotypes**: Body weight, organ weights, biochemical measurements
- **Gene expression**: RNA-seq derived expression levels (log2 CPM + inverse rank normalized)
- **Metabolites**: Metabolomics measurements from various platforms
- **Proteins**: Quantitative proteomics data
- **Any numeric trait**: Suitable for linear mixed model analysis

Example structure for multi-phenotype analysis:
```csv
,covariate,,,phenotype,,,
SampleID,Sex,generation,batch,Slc2a2,Apob,liver_weight
DOchln_001,male,1,B1,-0.245,1.432,1.24
DOchln_002,female,1,B1,0.891,-0.672,0.98
```

**Data Preprocessing Workflow**:
For molecular phenotypes (eQTL, mQTL analyses), the preparatory R Markdown workflow (`Metadata_phenotype_QTL2_NF.Rmd`) provides comprehensive automated preprocessing including:

**RNA-seq Expression Data (eQTL Analysis)**:
- **Filtering**: Low-expression transcript removal using edgeR::filterByExpr() to retain genes with sufficient expression across samples
- **Normalization**: Library size normalization via counts-per-million (CPM) transformation with prior count adjustment
- **Transformation**: Gene-wise inverse rank normalization (IRN) for optimal linear modeling and normal distribution approximation
- **Quality Control**: Sample correlation analysis and batch effect assessment

**Sample Processing**:
- **ID Standardization**: Automatic conversion of variant sample naming conventions (Dochln, DOCHLN, DOChln) to standardized format (DOchln_###)
- **Duplicate Removal**: Identification and removal of duplicate samples based on predefined exclusion lists
- **Consistency Validation**: Cross-validation of sample IDs between phenotype and genotype datasets
- **Missing Data Handling**: Appropriate encoding and documentation of missing values for downstream r/qtl2 analysis

**Metadata Integration**:
- **Covariate Validation**: Verification of required covariates (Sex, generation) for DO cross analysis
- **Batch Integration**: Optional batch effect covariates for experimental design control
- **Format Conversion**: Automated generation of r/qtl2-compatible CSV structure with proper header indicators

### Genotype File Format
**GeneSeek FinalReport Files**:
- **File Type**: Standard Illumina genotyping array output (FinalReport.txt)
- **Supported Arrays**: GigaMUGA, MiniMUGA, MEGA-MUGA platforms
- **Format**: Tab-delimited with standard GeneSeek column headers
- **Multiple Files**: Supports glob patterns for batch processing (e.g., 'Data/FinalReport*.txt')
- **Sample Matching**: Automatic cross-validation with phenotype sample IDs
- **Quality Control**: Built-in marker and sample filtering with comprehensive validation

## Pipeline Architecture

The pipeline consists of 9 numbered modules (plus optional Module 6b) organized into a sequential workflow that matches the results directory structure. This modular design enables precise resume capabilities and optimal resource allocation for large-scale QTL analysis:

### Module 1: Phenotype Processing (`01_phenotype_process.nf`)
**Purpose**: Comprehensive validation and preparation of phenotype/covariate data according to r/qtl2 specifications with robust support for diverse molecular data types and experimental designs

**Key Functions**:
- **Data Structure Validation**: Automated parsing of specially formatted CSV files with header row indicators for seamless data type recognition
- **Sample ID Standardization**: Cross-validation of sample identifiers between phenotype and genotype datasets with automatic format harmonization
- **Numeric Matrix Enforcement**: Conversion of phenotype data to numeric matrices with intelligent handling of missing data patterns and outlier detection
- **Required Covariate Validation**: Verification of essential covariates (Sex, generation) as mandated by r/qtl2 for Diversity Outbred cross analysis
- **Optional Covariate Integration**: Flexible incorporation of experimental factors (Diet, batch, age, treatment groups) with automatic categorical encoding
- **Multi-Phenotype Support**: Scalable validation for high-dimensional datasets including thousands of molecular traits (transcripts, metabolites, proteins)
- **Quality Control Diagnostics**: Comprehensive outlier detection, batch effect assessment, and correlation analysis with detailed visualization reports
- **Format Compliance**: Generation of r/qtl2-compliant data structures with proper metadata integration
- **Missing Data Management**: Intelligent handling of missing values with appropriate encoding strategies for downstream statistical analysis
- **Single-Sex Study Support**: When `--single_sex` is set, auto-adds the sex column if missing and protects `sex`, `ngen`, and `generation` from covariate exclusion; these columns are always written to the covariate file because r/qtl2 requires them as cross metadata (sex chromosome handling, DO cross creation), regardless of whether they vary across samples

**Input Format Requirements**:
The module processes specially formatted CSV files with the following convention:
- **Row 1**: Type indicators distinguishing covariate and phenotype columns
- **Row 2**: Column headers with standardized naming conventions
- **Data Rows**: Sample-wise measurements with automatic data type validation
- **Supported Scales**: Handles diverse phenotype scales from clinical measurements to high-dimensional molecular data

**Technical Implementation**: Uses advanced R data processing libraries (dplyr, data.table) with comprehensive error checking and logging for reliable data preparation in high-throughput environments.

**Resources**: 1 CPU, 2GB RAM, 1 hour

### Module 2: Genotype Processing (`02_genotype_process.nf`)
**Purpose**: Comprehensive conversion of GeneSeek FinalReport files to r/qtl2-compatible format with robust quality control and validation

**Key Functions**:
- **GeneSeek FinalReport Parsing**: Advanced parsing engine supporting multiple input files with automated batch processing and error recovery
- **Genotype Call Conversion**: Systematic conversion of SNP calls to A/H/B encoding scheme specifically optimized for r/qtl2 Diversity Outbred cross analysis
- **Chromosome-Specific Organization**: Automated generation of chromosome-specific genotype files (autosomes 1-19, sex chromosomes X and Y, mitochondrial chromosome M) with proper genetic map integration
- **Sample ID Harmonization**: Cross-validation and standardization of sample identifiers ensuring perfect consistency with phenotype datasets
- **Missing Data Management**: Intelligent handling of missing genotype calls with appropriate encoding strategies for downstream QTL analysis
- **Marker Quality Control**: Comprehensive assessment including call rate evaluation, Hardy-Weinberg equilibrium testing, and marker integrity validation
- **Reference Integration**: Automated download and integration of MGI (Mouse Genome Informatics) reference files and genetic maps for accurate genomic positioning
- **Array Platform Support**: Native compatibility with standard mouse genotyping arrays (GigaMUGA, MiniMUGA, MEGA-MUGA) with platform-specific optimization
- **Allele Code Processing**: Sophisticated allele coding system handling founder strain contributions with proper heterozygote detection
- **Validation Reporting**: Detailed quality control reports with summary statistics, marker distribution analysis, and validation metrics

**Technical Implementation**: Leverages high-performance R data processing with memory-optimized data structures for handling large-scale genotyping datasets (>65K markers, >1000 samples) efficiently.

**Resources**: 2 CPUs, 64GB RAM, 4 hours

### Module 3: Control File Generation (`03_control_file_generation.nf`)
**Purpose**: Create r/qtl2 control files and metadata structures for multiparental cross analysis

**Key Functions**:
- JSON control file generation with complete metadata for multiparental crosses
- File path specification and validation for execution environment compatibility
- Genotype encoding mappings (A/H/B) and sex chromosome handling for DO populations
- Cross type specification with founder strain definitions (8 founder strains for DO)
- Reference file integration including genetic maps, founder genotypes, and marker annotations
- Cross information structure creation compatible with r/qtl2 requirements
- Automated cross2 generation methodologies for robust object creation

**Resources**: 2 CPUs, 4GB RAM, 2 hours

### Module 4: Cross2 Object Creation (`04_cross2_creation.nf`)
**Purpose**: Generate validated r/qtl2 cross2 objects

**Key Functions**:
- Cross2 object creation using read_cross2()
- Data integrity validation and diagnostics
- X chromosome encoding for DO crosses
- Quality control reporting
- Object serialization for downstream analysis

**Resources**: 2 CPUs, 4GB RAM, 2 hours

### Module 5: Genome Scan Preparation (`05_prepare_genome_scan.nf`)
**Purpose**: Prepare high-resolution genotype probabilities and kinship matrices

**Key Functions**:
- Pseudomarker insertion at 0.5 cM resolution
- Genotype probability calculation (error_prob=0.002)
- LOCO (Leave One Chromosome Out) kinship matrix generation
- Genetic map preparation for high-resolution mapping
- Memory-efficient data structure optimization

**Resources**: 16 CPUs, 256GB RAM, 12 hours

### Module 6: HPC Batch-Based QTL Analysis (`06_qtl_analysis.nf`)
**Purpose**: Scalable batch processing architecture for genome-wide QTL scanning with intelligent resource management

**Key Functions**:
- **Batch Processing Architecture**: Three-tier workflow (Setup → Batch Processing → Results Combination)
- **Intelligent Chunking**: Organizes phenotypes into chunks of ~200, batched into groups of 10 for optimal parallelization
- **Setup Process**: Creates chunk and batch assignments, calculates optimal resource allocation (4 CPUs, 32GB, 1h)
- **Batch Processing**: Processes 10 chunks per batch using alleleprob for efficient memory usage (96 CPUs, 1.6TB, 12h per batch)
- **Results Integration**: Combines all batch outputs into comprehensive genome scan results (16 CPUs, 256GB, 2h)
- **Smart LOD Filtering**: Pre-filters phenotypes using find_peaks() before expensive permutation testing
- **Multiple Peak Detection**: Identifies ALL QTL peaks per phenotype (drop=1.5 LOD threshold) for comprehensive analysis
- **Comprehensive Logging**: Detailed per-batch logs with timestamps, node information, and processing statistics

**Single-Sex Covariate Handling**: Before calling `model.matrix()`, module 06 explicitly filters out any covariate columns that have fewer than two unique non-missing values (constant columns). This step is required because passing a constant factor such as `Sex` in an all-female study directly to `model.matrix()` crashes R with `contrasts can be applied only to factors with 2 or more levels`. The removed covariate names are written to the batch log so the filtering is fully auditable. After the constant-column filter, a secondary guard checks whether the resulting covariate data frame is empty; if so, `addcovar` is set to `NULL` and `scan1` runs without additive covariates. Sex remains in the cross object as required metadata for X chromosome handling but is dropped from the statistical model explicitly in this step rather than silently.

**Technical Innovation**: Enables processing of datasets with 10,000+ phenotypes through batched parallel execution, preventing infrastructure bottlenecks while maintaining memory efficiency

**Resources**: Setup (4 CPUs, 32GB), Batch Processing (96 CPUs, 1.6TB), Combination (16 CPUs, 256GB)

### Module 6b: Regional QTL Filtering (`06b_filter_peaks_by_region.nf`) **[OPTIONAL]**
**Purpose**: Dramatically reduce permutation testing computational load by filtering QTL peaks to specific genomic regions

**Key Functions**:
- **GTF-Based Gene Mapping**: Cross-references phenotype IDs with Ensembl GTF annotations to identify gene chromosomal locations
- **Flexible Region Specification**: Supports single regions (`chr12:81-91`), whole chromosomes (`chr12`), or multiple regions (`chr12:81-91;chr2:50-60`)
- **Intelligent Overlap Detection**: Identifies phenotypes whose genes overlap with target regions (Mbp coordinates)
- **Multi-Region Support**: Processes multiple genomic regions simultaneously for complex analyses
- **Workload Reduction Reporting**: Quantifies computational savings (typically 70-95% reduction in phenotypes for permutation)
- **Comprehensive Logging**: Detailed reports on matching statistics, unmatched phenotypes, and per-region summaries
- **Seamless Integration**: Filtered phenotype list automatically passed to Module 7 when `--qtl_region` parameter is specified

**Use Cases**:
- **Targeted eQTL Analysis**: Focus on cis-eQTLs for genes in specific chromosomal regions
- **Candidate Region Investigation**: Deep analysis of known QTL hotspots or disease-associated loci
- **Trans-band Exploration**: Examine trans-regulation from specific genomic regions
- **Resource-Constrained Studies**: Maximize statistical power while minimizing computational time

**Example**: For a chr12:81-91 Mbp region filter on 13,374 phenotypes, typically reduces to ~300-500 target phenotypes, achieving 96-98% workload reduction for permutation testing

**Resources**: 4 CPUs, 16GB RAM, 30 minutes

### Module 7: Coordinator-Based Permutation Testing (`07_permutation_testing.nf`)
**Purpose**: Reliable permutation testing for large phenotype sets using a coordinator architecture that survives long-running HPC jobs without depending on the Nextflow driver session remaining active

**Architecture**: `PERM_SETUP → PERM_COORDINATOR → PERM_AGGREGATE`

The previous design submitted individual SLURM batch jobs directly through Nextflow's job manager. If the Nextflow driver process died (for example, when a VS Code session hit a 10-hour wall time), all job monitoring stopped mid-run. The coordinator architecture eliminates this failure mode.

**Process: PERM_SETUP**
- Loads the full cross2 object and filters it to the phenotypes that passed the LOD threshold from Module 6 (or regional filtering from Module 6b)
- Writes the filtered cross2 object and the phenotype list to the published output directory
- Calculates and logs the full batch configuration: 200 permutations × 5 batches per phenotype = 1,000 total permutations per phenotype

**Process: PERM_COORDINATOR**
- A **single SLURM job with a 7-day wall time** (1 CPU, 4 GB) that runs without a container so it has access to host SLURM tools (`sbatch`, `sacct`)
- Writes the shared R worker script (`perm_batch_worker.R`) to the published directory so all compute nodes can access it over the shared filesystem
- Submits all permutation batch jobs directly via `sbatch` (48 CPUs, 200 GB, 2-hour wall time per batch job)
- Each batch worker runs inside the Singularity container: `apptainer exec --bind /project:/project singularity_cache/dpbudke-qtl2-pipeline-latest.img`
- Polls job status every 5 minutes via `sacct` and retries failed batches up to 2 times automatically
- Has internal resume logic: skips any batch whose output file already exists and is non-empty
- Writes a `COORDINATOR_COMPLETE` sentinel file when all batches finish, which triggers `PERM_AGGREGATE`
- Because the coordinator is itself a 7-day SLURM job, it is completely independent of the Nextflow driver session; driver interruption does not stop the batch work

**Process: PERM_AGGREGATE**
- Triggered by the `COORDINATOR_COMPLETE` sentinel file
- Reads batch result files directly from the published batch directory on the shared filesystem rather than staging them as Nextflow inputs (avoids staging thousands of individual files)
- Combines all batch permutation matrices into a single 1,000-permutation-per-phenotype result matrix
- Calculates multi-level empirical significance thresholds (63%, 90%, 95%, 99% quantiles)
- Writes the list of phenotypes passing the LOD threshold at the 99% level for use in Module 8

**Key Properties**:
- **Fork-based parallelism**: Each batch worker uses `scan1perm()` with `cores=48` and fork (not PSOCK), which is efficient on Linux
- **Automatic resume**: The coordinator skips batches whose output files already exist; re-submitting the coordinator job after an interruption safely continues from where work left off
- **Error handling**: Batches in FAILED, TIMEOUT, CANCELLED, NODE_FAIL, or OUT_OF_MEMORY states are retried up to 2 times; permanently failed batches are logged and the coordinator exits with a non-zero status

**Efficiency Achievement**: LOD threshold of 7.0 typically reduces the permutation phenotype set by 70-90% relative to the full scan; regional filtering (Module 6b) can achieve 95-98% reduction

**Resources**: Setup (4 CPUs, 32 GB, 10 min), Coordinator (1 CPU, 4 GB, 7-day wall time), Batch Workers (48 CPUs, 200 GB, 2 h each), Aggregation (4 CPUs, 32 GB, 1 h)

### Module 8: Significant QTL Identification (`08_identify_significant_qtls.nf`)
**Purpose**: Statistical QTL calling using empirical permutation thresholds with multi-level significance classification

**Key Functions**:
- **Filtered Cross2 Integration**: Uses phenotype-filtered cross2 object from Module 7 for precise QTL identification
- **Threshold Application**: Applies phenotype-specific empirical significance thresholds from permutation testing
- **Multi-level Classification**: Identifies QTLs at four significance levels: 63% (suggestive), 90% (α=0.10), 95% (α=0.05), 99% (α=0.01)
- **Multiple Peaks Per Phenotype**: Captures all significant peaks using find_peaks() with drop=1.5 LOD criterion
- **Comprehensive Annotation**: Annotates each QTL with all four threshold values for the corresponding phenotype
- **Chromosome Distribution**: Summarizes QTL counts across all chromosomes
- **High-Priority Flagging**: Identifies exceptional QTLs with LOD ≥ 10 for immediate follow-up
- **Publication-Ready Output**: CSV format with phenotype, chromosome, position, LOD score, and significance annotations

**Output Files**:
- **significant_qtls.csv**: Complete QTL catalog with multi-level significance annotations
- **qtl_summary.txt**: Statistical overview, chromosome distribution, and top QTL highlights
- **validation_report.txt**: Processing log with data alignment verification and QTL counts

**Resources**: 8 CPUs, 128GB RAM, 4 hours

### Module 9: QTL Visualizations (`09_visualize.nf`)
**Purpose**: Generate publication-quality QTL coefficient plots for highly significant QTLs

**Key Functions**:
- **High-Confidence QTL Focus**: Visualizes only 99% significant QTLs from Module 8 for maximum stringency
- **Founder Allele Effect Plots**: Generates LOD score plots, scan1coef() founder allele coefficient plots, and BLUP effect plots for each significant QTL
- **Coefficient Calculation**: Computes scan1coef() founder allele effects and scan1blup() shrinkage estimates using alleleprobs
- **Chromosome Organization**: Automatically organizes plots into chromosome-specific subdirectories (chr1/ through chr19/ and chrX/)
- **Batch Processing**: Efficiently processes thousands of QTLs with progress monitoring
- **Comprehensive Naming**: Each plot filename includes phenotype, chromosome, position (Mb), LOD score, and plot type for easy identification
- **Validation Reporting**: Detailed summary of plots generated per chromosome with success/failure tracking

**Output Organization**:
- **chr{1-19,X}/ subdirectories**: Plots organized by chromosome for easy navigation
- **Filename format**: `{phenotype}_chr{chr}_{pos}Mb_LOD{lod}_{type}.png`
  - `{type}` is one of: `lod`, `coef_effects`, `blup_effects`
- **Example**: `TMAO_base_ln_chr12_86.1Mb_LOD27.92_lod.png`
- **GxE analyses** (`--interactive_covar`): Plots are additionally stratified by covariate level (e.g., `ain76a_coef_effects.png`, `hc_coef_effects.png`, `full_coef_effects.png`)
- **Setup files**: `viz_batch_qtls.csv`, `viz_batch_ids.txt`, `viz_setup_log.txt`

**Use Cases**:
- **Publication Figures**: High-quality visualizations for manuscripts and presentations
- **Allele Effect Interpretation**: Visual assessment of founder strain contributions to trait variation
- **QTL Validation**: Quick visual confirmation of strong, biologically meaningful QTL signals
- **Follow-up Prioritization**: Identify QTLs with interesting allele effect patterns for further investigation

**Resources**: SETUP 2 CPUs / 8 GB / 30 min; BATCH 16 CPUs / 128 GB (scales on retry) / 2 h; AGGREGATE 2 CPUs / 8 GB / 30 min

### Module 10: TIMBR Allelic Series Analysis (`10_timbr.nf`)
**Purpose**: Infer the number of functional alleles among the 8 DO founder strains at each significant QTL using Bayesian MCMC (TIMBR; Crouse et al. 2020, *Genetics* 216:957–983)

**Key Functions**:
- **Allelic Series Inference**: Determines how many of the 8 founder haplotypes carry functionally distinct alleles at each QTL
- **Tree-Informed Prior**: Uses the DO/CC founder phylogenetic tree as a Dirichlet Process prior, improving power over flat priors
- **Configurable Significance Threshold**: Filters QTLs at `--timbr_sig_level` (default: 99%; options: 63%, 90%, 95%, 99%)
- **Batch MCMC**: Parallelizes TIMBR() calls across SLURM array jobs; default 10,000 MCMC samples per QTL
- **Four Diagnostic Plots per QTL**: Haplotype effects, allele count distribution, top allelic series patterns, and TIMBR-vs-qtl2 founder effect comparison
- **Master Summary Table**: Aggregates ln(BF), MAP allelic series pattern, functional allele count, and posterior probability across all analyzed QTLs

**Process Workflow**:
1. **TIMBR_SETUP**: Filters and deduplicates significant QTLs, assigns batches, pre-computes additive covariate matrix
2. **TIMBR_BATCH**: Runs TIMBR() per QTL; extracts MAP allelic series, counts functional allele groups, saves optional RDS object and 4 PNG plots
3. **TIMBR_AGGREGATE**: Combines batch CSVs into master summary; reports functional allele distribution and median ln(BF)

**Key Parameters**:
- `--timbr_sig_level` — Minimum significance threshold (default: `"99%"`)
- `--timbr_qtls_per_batch` — QTLs per SLURM batch (default: `20`, matched to 20 CPUs for mclapply)
- `--timbr_samples` — MCMC iterations per QTL (default: `10000`)
- `--run_timbr` — Set to `false` to skip this module (default: `true`)

**Output Files**:
- `{prefix}_timbr_master_summary.csv` — All QTL results: `lodcolumn, chr, pos, peak_Mb, lod, significance_level, nearest_marker, ln_BF, n_functional_alleles, top_series_pattern, top_series_prob, n_samples_used, status`
- `{prefix}_timbr_run_report.txt` — Functional allele distribution, results by significance level, median ln(BF)
- `chr{N}/{phenotype}_chr{N}_{pos}Mb/{phenotype}_chr{N}_{pos}Mb.rds` — Full TIMBR result object (optional)
- `chr{N}/{phenotype}_chr{N}_{pos}Mb/plots/01_haplotype_effects.png` — Posterior haplotype effects
- `chr{N}/{phenotype}_chr{N}_{pos}Mb/plots/02_n_alleles_distribution.png` — Posterior functional allele count
- `chr{N}/{phenotype}_chr{N}_{pos}Mb/plots/03_top_allelic_series.png` — Top 10 allelic series patterns
- `chr{N}/{phenotype}_chr{N}_{pos}Mb/plots/04_allelic_series_vs_founder_effects.png` — TIMBR vs qtl2 ROP comparison

**Resources**: SETUP 4 CPUs / 32 GB / 30 min; BATCH 20 CPUs / 64 GB (128 GB on retry) / 2 h; AGGREGATE 4 CPUs / 16 GB / 30 min

### Optional: Broman Mixup QC (`broman_mixup_qc.nf`)
**Purpose**: Detect sample mix-ups by comparing observed expression patterns with genotype-predicted expression using eQTL signatures

**Scientific Methodology**:
Based on Broman et al. (2015) G3 and Westra et al. (2011), this module identifies sample mislabeling by testing whether observed gene expression matches the expression pattern predicted from an individual's genotype at strong eQTL positions.

**Process Workflow**:

1. **eQTL Selection - Creating the Genotype-Expression Signature**:
   - Takes significant QTLs from Module 8 and sorts by LOD score (highest first)
   - Retains **only one peak per unique gene** to avoid redundancy (strongest eQTL per gene)
   - Selects top N genes (default: 100) with the strongest genotype-expression relationships
   - **Rationale**: Strong eQTL positions are where genotype robustly predicts expression, providing reliable "fingerprint" information
   - If fewer than N unique genes are available, uses all available eQTL
   - LOD range of selected eQTL is logged for quality assessment

2. **Observed Expression Extraction - The Measured Data**:
   - Creates an observed expression matrix: **samples × 100 genes**
   - Each value represents the **actual measured expression** from RNA-seq for that gene in that sample
   - This serves as the "ground truth" molecular phenotype data
   - Matrix dimensions validated and logged (e.g., "200 samples × 100 genes")
   - **Example**: If Sample001 has normalized expression of 2.34 for Gene1, this exact value is stored in the observed matrix

3. **Predicted Expression Calculation - Genotype-Based Prediction**:

   **For each of the 100 selected genes:**

   a. **Extract genotype probabilities at eQTL position** using `qtl2::pull_genoprobpos()`:
      - For a gene with eQTL on chr12 at 85.3 Mbp, extracts each sample's **founder haplotype probabilities** at that exact genomic position
      - DO mice have 8 founder strains (A/J, C57BL/6J, 129S1, NOD, NZO, CAST, PWK, WSB), yielding 8 probabilities per sample
      - **Example**: Sample001 might have [A/J: 0.05, C57BL/6J: 0.70, 129S1: 0.03, NOD: 0.01, NZO: 0.02, CAST: 0.15, PWK: 0.02, WSB: 0.02]
      - Interpretation: This sample is 70% likely to carry the C57BL/6J allele and 15% likely to carry the CAST allele at this eQTL

   b. **Fit single-QTL model to estimate founder allele effects** using `qtl2::fit1()`:
      - Performs linear regression: `Expression ~ Genotype` using **all samples' data**
      - Parameter `zerosum=FALSE` preserves the original expression scale
      - Estimates **8 founder allele effects** representing the mean expression level for each founder haplotype
      - **Example output**: [A/J: 1.23, C57BL/6J: 2.87, 129S1: 0.45, NOD: 1.12, NZO: 0.89, CAST: -1.34, PWK: 0.67, WSB: 1.05]
      - **Biological interpretation**: At this eQTL, carrying the C57BL/6J allele increases expression to 2.87, while carrying the CAST allele decreases it to -1.34

   c. **Calculate predicted expression** via matrix multiplication:
      - For each sample, multiplies their genotype probabilities by the founder allele effects
      - Formula: `predicted_expression = genotype_probabilities %*% founder_effects`
      - **Worked example for Sample001**:
        ```
        Predicted = (0.05 × 1.23) + (0.70 × 2.87) + (0.03 × 0.45) + (0.01 × 1.12) +
                    (0.02 × 0.89) + (0.15 × -1.34) + (0.02 × 0.67) + (0.02 × 1.05)
                  = 0.06 + 2.01 + 0.01 + 0.01 + 0.02 - 0.20 + 0.01 + 0.02
                  = 1.94
        ```
      - This process repeats for all 100 genes, creating a **predicted expression matrix: samples × 100 genes**
      - **Key concept**: Each sample's predicted expression is a weighted average based on which founder alleles they carry

4. **Distance Matrix Computation - Comparing Observed vs Predicted Patterns**:
   - Uses `lineup2::dist_betw_matrices()` to calculate **RMS (root mean square) distances**
   - Compares each sample's **observed expression** (100 real measurements) with every sample's **predicted expression** (100 genotype-based predictions)
   - Creates a **samples × samples distance matrix** where entry (i,j) represents how well Sample i's observed expression matches Sample j's genotype-predicted expression
   - **Distance formula**: `RMS = sqrt(mean((observed_i - predicted_j)²))` across all 100 genes
   - **Matrix interpretation**:
     - **Diagonal elements** (e.g., Sample001 obs vs Sample001 pred): "self-distance" = how well that sample's genotype predicts its own expression
     - **Off-diagonal elements** (e.g., Sample001 obs vs Sample002 pred): distance to other samples' predicted patterns
   - **Example distance matrix**:
     ```
                     Sample001_pred  Sample002_pred  Sample003_pred
     Sample001_obs        0.23           1.87           2.45
     Sample002_obs        1.92           0.18           2.11
     Sample003_obs        2.34           2.09           0.31
     ```
   - **Normal pattern**: Diagonal values (self-distances) should be smallest in each row
   - **Mix-up pattern**: Off-diagonal value is smaller than diagonal, indicating observed expression matches a different sample's genotype

5. **Mixup Detection - Identifying Mislabeled Samples**:
   - **Extracts self-distances** using `lineup2::get_self()`: Diagonal values showing how well each sample's genotype predicts its own expression
   - **Identifies minimum distances** using `lineup2::get_best()`: For each sample, finds the smallest distance to ANY sample's predicted expression (including self)
   - **Flags potential mix-ups**: Samples where `self_distance > min_distance`
   - **Interpretation**:
     - **Normal sample**: Self-distance ≈ minimum distance (genotype correctly predicts own expression)
     - **Mix-up detected**: Self-distance > minimum distance (observed expression matches a different sample's genotype better than own genotype)
   - **Example normal sample**: Sample001 has self_distance = 0.23, min_distance = 0.23, best_match = Sample001
   - **Example mix-up**: Sample042 has self_distance = 3.12, min_distance = 0.21, best_match = Sample087 → Sample labels likely swapped
   - All per-sample statistics exported to CSV for manual review and investigation

**Input Requirements**:
- **cross2_file**: Cross2 object from Module 4 (`results/04_cross2_creation/{prefix}_cross2_object.rds`)
- **genoprob_file**: Genotype probabilities from Module 5 (`results/05_genome_scan_preparation/{prefix}_genoprob.rds`)
- **expr_data_file**: Expression data matrix (RDS or CSV format: samples × genes, typically inverse-rank-normalized log2 CPM)
  - RDS: Matrix object with samples as rows, genes as columns
  - CSV: Comma-separated file with sample IDs as row names, gene IDs as column headers
- **expr_peaks_file**: eQTL peaks CSV from Module 8 (`results/08_significant_qtls/{prefix}_significant_qtls.csv`)
- **n_top_eqtl**: Number of top eQTL to analyze (default: 100, adjustable based on dataset size)

**Output Files** (`results/mixup_qc/`):
- **{prefix}_mixup_report.html**: Minimal HTML report displaying count of samples with potential problems
- **{prefix}_mixup_problems.csv**: Per-sample analysis with columns:
  - `sample`: Sample identifier
  - `self_distance`: RMS distance between observed and own predicted expression
  - `min_distance`: Minimum RMS distance to any sample's predicted expression
  - `best_match`: Sample ID with closest predicted expression match
  - `is_problem`: Boolean flag (TRUE if self_distance > min_distance, indicating potential mixup)
- **{prefix}_distance_matrix.csv**: Full RMS distance matrix (observed samples × predicted samples)
- **{prefix}_mixup_plots.jpg**: Diagnostic scatter plot (self distance vs minimum distance, 300 DPI JPEG format)
  - Points below diagonal (red) indicate potential mix-ups where self-distance exceeds minimum distance
  - Gray points represent samples with consistent genotype-expression relationships
  - Blue dashed line represents y = x diagonal (threshold for problem detection)
- **mixup_qc_log.txt**: Detailed analysis log with timestamps and processing statistics

**Technical Implementation Notes**:
- Module automatically installs `lineup2` package from CRAN if not available in the R environment
- Uses writable R library path in Nextflow work directory to handle package installation
- Handles both RDS and CSV input formats for expression data with automatic format detection
- Extracts genetic map from cross2 object (gmap preferred, falls back to pmap if unavailable)

**Interpretation Guidelines**:
- **Normal samples**: Self-distance ≈ minimum distance (points near diagonal on scatter plot)
- **Potential mix-ups**: Self-distance > minimum distance (points below diagonal, colored red)
- **Action required**: Investigate `best_match` column to identify likely correct sample identity
- **Threshold considerations**: Small deviations may reflect biological variation; focus on substantial differences where self-distance substantially exceeds minimum distance

**Usage**: Run as standalone QC analysis after Module 8 when expression data is available

**Resources**: 8 CPUs, 64GB RAM, 4 hours

## Installation and Setup

### Quick Setup

**1. Clone the repository**
```bash
git clone <repo-url> QTL2_NF
cd QTL2_NF
```

**2. Pull the container image (required)**

The `singularity_cache/` directory is gitignored because the `.img` file is a large binary. Create the directory and pull the image manually:
```bash
mkdir singularity_cache
apptainer pull singularity_cache/dpbudke-qtl2-pipeline-latest.img \
    docker://dpbudke/qtl2-pipeline:latest
```

This is the only required setup step. No configuration changes are needed.

**3. No path configuration required**

All paths in `nextflow.config` use `${projectDir}`, a Nextflow built-in variable that automatically resolves to the directory containing `main.nf`. The pipeline works correctly regardless of what the cloned directory is named or where it lives on the filesystem.

**Note on gitignored directories**: `nextflow_work/` (Nextflow work directory) and `singularity_cache/` (container image) are both gitignored and created locally. They are not tracked in version control.

### Prerequisites

**Software Requirements**:
- **Nextflow** ≥ 23.04.0
- **Apptainer/Singularity** (for HPC environments)
- **Docker** (for local development)
- **SLURM** workload manager with array job support

**Container Requirements**:
- Access to Docker Hub or pre-pulled container: `dpbudke/qtl2-pipeline:latest`
- Singularity cache directory with sufficient storage (≥50GB)

### Quick Start Guide

1. **Environment Setup**:
```bash
# Configure HPC settings
# Edit nextflow.config to match your SLURM account and resource allocations
```

2. **Data Preparation**:
```bash
# Create data directory structure
mkdir -p Data/

# Place input files:
# - Phenotype CSV (specially formatted with header row indicators) → Data/
# - GeneSeek FinalReport files → Data/
# - Verify file naming matches glob patterns
ls Data/  # Should show your phenotype.csv and FinalReport*.txt files
```

3. **Pipeline Execution**:
```bash
# Full pipeline execution
nextflow run main.nf \
  --phenotype_file Data/your_phenotypes.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix YourStudy \
  -profile standard

# Monitor progress
tail -f .nextflow.log
```

## Usage Examples

### Full Pipeline Execution

**Standard Multi-Phenotype QTL Analysis**:
```bash
# Complete eQTL analysis with comprehensive parameter set
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/*FinalReport*.txt' \
  --study_prefix DOchln \
  --lod_threshold 7.0 \
  --outdir results/ \
  -profile standard
```

**Single-Sex Studies**:
```bash
# All-female study: sex column is retained in the cross object for X chromosome handling
# but is filtered out of the statistical model by module 06 before model.matrix() is called
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_females \
  --lod_threshold 7.0 \
  --single_sex female \
  -profile standard

# All-male study
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_males \
  --lod_threshold 7.0 \
  --single_sex male \
  -profile standard

# Single-sex study where the sex column is absent from the input file:
# --single_sex auto-adds the column with the correct value before cross object creation
nextflow run main.nf \
  --phenotype_file Data/females_only_no_sex_col.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_females_nosexcol \
  --lod_threshold 7.0 \
  --single_sex female \
  -profile standard
```

**Sample Filtering Applications**:
```bash
# Male-specific analysis on high-cholesterol diet
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_males_hc \
  --lod_threshold 7.0 \
  --sample_filter '{"Sex":["male"],"Diet":["hc"]}' \
  --outdir results/sex_diet_specific \
  -profile standard

# Batch-specific analysis (removing batch effects)
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_batch1 \
  --lod_threshold 6.0 \
  --sample_filter '{"batch":["B1","B2"]}' \
  -profile standard

# Generation-specific analysis for population structure control
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_G6plus \
  --lod_threshold 7.0 \
  --sample_filter '{"generation":["6","7","8","9"]}' \
  -profile standard
```

**Multiple Array Processing**:
```bash
# Processing multiple FinalReport files from different batches
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport_Batch*.txt' \
  --study_prefix DOchln_multibatch \
  --lod_threshold 7.0 \
  -profile standard
```

**Regional QTL Filtering (Module 6b - Targeted Analysis)**:
```bash
# Single genomic region (chr12:81-91 Mbp)
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_chr12_cis \
  --lod_threshold 7.0 \
  --qtl_region "chr12:81-91" \
  --gtf_file Data/Mus_musculus.GRCm39.105.gtf.gz \
  -profile standard

# Multiple genomic regions
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_multiregion \
  --lod_threshold 7.0 \
  --qtl_region "chr12:81-91;chr2:50-60;chr7:10-20" \
  --gtf_file Data/Mus_musculus.GRCm39.105.gtf.gz \
  -profile standard

# Whole chromosome analysis
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_chr12_all \
  --lod_threshold 7.0 \
  --qtl_region "chr12" \
  --gtf_file Data/Mus_musculus.GRCm39.105.gtf.gz \
  -profile standard
```

### Resume Capabilities

The pipeline provides precise resume functionality using `main_resume.nf`, enabling restart from any processing module:

```bash
# Resume from genome scan preparation (skip data processing modules)
nextflow run main_resume.nf \
  --resume_from prepare_scan \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from HPC batch genome scanning (most common restart point)
nextflow run main_resume.nf \
  --resume_from genome_scan \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from permutation testing (after scan completion)
# NOTE: For permutation runs, prefer submitting via sbatch wrapper — see "Running Permutation Testing" below
nextflow run main_resume.nf \
  --resume_from permutation \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume only the PERM_AGGREGATE step after all coordinator batches have finished
nextflow run main_resume.nf \
  --resume_from perm_aggregate \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from QTL visualization only
nextflow run main_resume.nf \
  --resume_from visualize \
  --study_prefix DOchln
```

**Resume Options**:
- `phenotype` - Start from phenotype processing and validation
- `genotype` - Start from genotype processing and chromosome file generation
- `control` - Start from control file generation and metadata creation
- `cross2` - Start from cross2 object creation and validation
- `prepare_scan` - Start from genome scan preparation (genotype probabilities)
- `genome_scan` - Start from HPC batch genome scanning and peak detection
- `regional_filter` - Start from regional QTL filtering (requires `--qtl_region`)
- `permutation` - Start from permutation testing (runs full PERM_SETUP → PERM_COORDINATOR → PERM_AGGREGATE)
- `perm_aggregate` - Run only PERM_AGGREGATE, reading completed batch files from the published batch directory; use this when the coordinator has finished but the aggregation step failed or was not yet run
- `significant_qtls` - Start from significant QTL identification and summarization
- `visualize` - Start from QTL visualization (plot_coefCC for 99% significant QTLs)
- `timbr` - Start from TIMBR allelic series analysis

**Note on `perm_aggregate`**: This resume case passes the phenotype list file as a trigger rather than staging thousands of batch files as Nextflow inputs. `PERM_AGGREGATE` reads batch files directly from the published `07_permutation_testing/batches/` directory on the shared filesystem.

### Development and Testing

**Test Mode Execution**:
```bash
# Positive control: chromosome 2 coat color validation
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix CoatColorTest \
  --test_mode \
  -profile standard

# Test mode features:
# - Processes chromosome 2 only (known coat color QTL)
# - Uses coat_color as primary phenotype
# - Reduced computational time for validation
# - Expected LOD peak at ~100 cM on chromosome 2
```

### Quality Control Module (Optional)

**Standalone Sample Mixup Detection**:
```bash
# Comprehensive mixup detection with existing QTL results
nextflow run modules/do_mixup_qc.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --cross2_object results/04_cross2_creation/DOchln_cross2_object.rds \
  --qtl_results results/06_qtl_analysis/DOchln_scan_results.rds \
  --study_prefix DOchln_mixup_check \
  --outdir results/quality_control

# Basic mixup detection using cis-eQTL approach
nextflow run modules/do_mixup_qc.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --cross2_object results/04_cross2_creation/DOchln_cross2_object.rds \
  --study_prefix DOchln_cis_mixup \
  --outdir results/quality_control

# Integration with main pipeline
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln \
  --run_mixup_qc true \
  -profile standard
```

## Output Structure

Results are organized in numbered directories matching the modular pipeline structure, enabling clear progress tracking and precise resume capabilities:

```
results/
├── 01_phenotype_processing/          # Validated phenotype and covariate matrices
│   ├── {prefix}_phenotypes.csv      # r/qtl2 compatible phenotype data
│   ├── {prefix}_covariates.csv      # Covariate matrix (Sex, generation, batch)
│   ├── {prefix}_valid_samples.txt   # Sample ID validation report
│   └── phenotype_validation_report.txt
├── 02_genotype_processing/           # Chromosome-specific genotype files
│   ├── {prefix}_geno_chr{1-19,X,Y}.csv  # Per-chromosome genotype matrices
│   ├── {prefix}_founder_genos.csv   # Founder strain genotypes
│   ├── {prefix}_allele_codes.csv    # Reference allele codes
│   └── genotype_validation_report.txt
├── 03_control_file_generation/       # r/qtl2 control files and metadata
│   ├── {prefix}_control.json        # Master r/qtl2 control file
│   ├── genetic_map.csv              # Genetic map positions
│   ├── physical_map.csv             # Physical map positions
│   └── control_generation_report.txt
├── 04_cross2_creation/              # Cross2 objects and validation
│   ├── {prefix}_cross2_object.rds   # Complete r/qtl2 cross2 object
│   ├── cross2_summary.txt           # Object validation summary
│   └── cross2_validation_report.txt
├── 05_genome_scan_preparation/       # High-resolution mapping components
│   ├── {prefix}_genoprob.rds        # Genotype probabilities (0.5 cM resolution)
│   ├── {prefix}_kinship.rds         # LOCO kinship matrices
│   ├── {prefix}_genetic_map.rds     # High-resolution genetic map
│   └── preparation_validation_report.txt
├── 06_qtl_analysis/                 # HPC batch genome scan results
│   ├── batches/                      # Individual batch processing results
│   │   ├── {prefix}_batch1_results.rds      # Scan results for batch 1 (10 chunks)
│   │   ├── {prefix}_batch1_peaks.csv        # Peaks found in batch 1
│   │   ├── {prefix}_batch1_filtered_phenos.txt  # Phenotypes with peaks in batch 1
│   │   ├── {prefix}_batch1_log.txt          # Batch 1 processing log
│   │   └── ...                              # Additional batch results
│   ├── {prefix}_pheno_chunks.txt            # Chunk assignments (chunk_id, phenotype)
│   ├── {prefix}_batch_assignments.txt       # Batch assignments (batch_id, chunk_id)
│   ├── {prefix}_chunk_summary.txt           # Chunking strategy summary
│   ├── {prefix}_scan_results.rds            # Combined comprehensive scan results
│   ├── {prefix}_all_peaks.csv               # All preliminary peaks (LOD ≥ threshold)
│   ├── {prefix}_filtered_phenotypes.txt     # Phenotypes passing LOD filter
│   └── {prefix}_batch_processing_summary.txt  # Batch processing summary
├── 06b_regional_filtering/          # Optional regional QTL filtering (if --qtl_region specified)
│   ├── {prefix}_regional_filtered_phenotypes.txt  # Phenotypes in target region(s)
│   ├── {prefix}_filtered_peaks.csv          # Peaks with gene location annotations
│   └── {prefix}_regional_filtering_report.txt  # Filtering statistics and workload reduction
├── 07_permutation_testing/          # Coordinator-based permutation testing
│   ├── batches/                     # Individual batch permutation results (written by batch workers)
│   │   ├── {prefix}_{phenotype}_batch_1.rds  # Permutation batch 1 for each phenotype
│   │   ├── {prefix}_{phenotype}_batch_2.rds  # Permutation batch 2 for each phenotype
│   │   └── ...                               # 5 batches × N phenotypes (200 perms each = 1000 total)
│   ├── coordinator_job_logs/        # Per-batch SLURM job stdout/stderr logs
│   │   └── {prefix}_{phenotype}_b{N}.log     # Log from each individual batch worker job
│   ├── perm_batch_worker.R          # Shared R worker script written by coordinator
│   ├── {prefix}_coordinator_log.txt # Coordinator progress log (submissions, retries, completion)
│   ├── COORDINATOR_COMPLETE         # Sentinel file written when all batches finish
│   ├── {prefix}_filtered_cross2.rds # Phenotype-filtered cross2 object (from PERM_SETUP)
│   ├── {prefix}_phenotype_list.txt  # List of phenotypes submitted for permutation
│   ├── {prefix}_setup_log.txt       # Permutation setup log
│   ├── {prefix}_fullPerm.rds        # Complete permutation matrix (1000 perms × N phenotypes)
│   ├── {prefix}_permThresh.rds      # Multi-level significance thresholds (63%, 90%, 95%, 99%)
│   ├── {prefix}_filtered_phenotypes_lod{threshold}.txt  # Phenotypes passing final threshold
│   └── {prefix}_aggregation_log.txt # Permutation aggregation log
├── 08_significant_qtls/             # Final QTL identification and analysis
│   ├── {prefix}_significant_qtls.csv    # Comprehensive significant QTL list
│   ├── {prefix}_qtl_summary.txt     # Summary statistics and distributions
│   ├── {prefix}_high_priority_qtls.csv  # Exceptional QTLs (LOD ≥ 10)
│   └── qtl_identification_report.txt
├── 09_visualize/                    # QTL visualization plots
│   ├── viz_batch_qtls.csv           # QTLs with batch assignments
│   ├── viz_setup_log.txt            # Setup report
│   ├── chr1/
│   │   ├── {phenotype}_chr1_{pos}Mb_LOD{lod}_lod.png
│   │   ├── {phenotype}_chr1_{pos}Mb_LOD{lod}_coef_effects.png
│   │   └── {phenotype}_chr1_{pos}Mb_LOD{lod}_blup_effects.png
│   ├── ...                          # Additional chromosomes
│   ├── chrX/
│   │   └── {phenotype}_chrX_{pos}Mb_LOD{lod}_{type}.png
│   └── visualization_report.txt     # Plot generation summary
├── 10_timbr/                        # TIMBR allelic series analysis
│   ├── {prefix}_timbr_qtls.csv      # QTLs filtered for TIMBR analysis
│   ├── {prefix}_timbr_batch_map.tsv # Batch assignments
│   ├── {prefix}_addcovar.rds        # Pre-computed covariate matrix
│   ├── {prefix}_timbr_setup_log.txt # Setup log
│   ├── {prefix}_timbr_master_summary.csv  # Master results table (all QTLs)
│   ├── {prefix}_timbr_run_report.txt      # Run statistics and allele distribution
│   ├── chr{N}/                      # Per-chromosome QTL subdirectories
│   │   └── {phenotype}_chr{N}_{pos}Mb/
│   │       ├── {phenotype}_chr{N}_{pos}Mb.rds   # Full TIMBR result object
│   │       └── plots/
│   │           ├── 01_haplotype_effects.png
│   │           ├── 02_n_alleles_distribution.png
│   │           ├── 03_top_allelic_series.png
│   │           └── 04_allelic_series_vs_founder_effects.png
│   └── ...
└── pipeline_info/                   # Execution monitoring and reporting
    ├── execution_timeline.html      # Interactive execution timeline
    ├── execution_report.html        # Comprehensive resource usage report
    ├── execution_trace.txt          # Detailed process execution trace
    └── pipeline_dag.svg             # Visual workflow representation
```

## Configuration

### HPC Environment Setup
Update `nextflow.config` for your SLURM HPC environment:

```groovy
process {
    executor = 'slurm'
    clusterOptions = '--account=your_project_account'  // Update for your allocation

    // Global container for all processes
    container = 'dpbudke/qtl2-pipeline:latest'

    // Module-specific resource allocations
    withName: PREPARE_GENOME_SCAN {
        cpus   = 32
        memory = '128 GB'
        time   = '8h'
    }

    withName: GENOME_SCAN_CHUNK {
        cpus   = 16
        memory = '64 GB'
        time   = '4h'
    }

    withName: PERMUTATION_TEST {
        cpus   = 48
        memory = '1 TB'
        time   = '72h'
    }
}

// Apptainer configuration for HPC
apptainer {
    enabled = true
    autoMounts = true
    cacheDir = '/path/to/your/singularity_cache'
}
```

### Pipeline Parameters

**Required Parameters**:
- `--phenotype_file`: Path to master phenotype CSV file (specially formatted with header row indicators)
- `--study_prefix`: Unique study identifier for all output files

**Optional Parameters**:
- `--finalreport_files`: GeneSeek FinalReport file(s) with glob pattern support (e.g., 'Data/FinalReport*.txt')
- `--lod_threshold`: LOD score threshold for permutation pre-filtering (default: 7.0)
- `--outdir`: Output directory path (default: results)
- `--sample_filter`: JSON string for sample subsetting (e.g., '{"Sex":["male"],"Diet":["hc"]}')
- `--single_sex`: Set to `'female'` or `'male'` for single-sex studies (default: null). See details below.
- `--test_mode`: Development mode - process chromosome 2 only (default: false)

**Regional Filtering Parameters (Module 6b)**:
- `--qtl_region`: Genomic region(s) for targeted QTL analysis (e.g., 'chr12:81-91' or 'chr12:81-91;chr2:50-60')
- `--gtf_file`: Path to Ensembl GTF annotation for gene location mapping (default: 'Data/Mus_musculus.GRCm39.105.gtf.gz')

**Resume Parameters**:
- `--resume_from`: Resume from specific module (01-09, or module names)

**Quality Control Parameters**:
- `--run_mixup_qc`: Enable standalone sample mixup detection (default: false)

**Example Parameter Sets**:
```bash
# Full eQTL analysis
--phenotype_file Data/eQTL_phenotypes.csv \
--finalreport_files 'Data/FinalReport*.txt' \
--study_prefix DOchln_eQTL \
--lod_threshold 7.0

# Male-only high-cholesterol diet analysis
--phenotype_file Data/clinical_traits.csv \
--finalreport_files 'Data/FinalReport*.txt' \
--study_prefix DOchln_males_hc \
--sample_filter '{"Sex":["male"],"Diet":["hc"]}' \
--lod_threshold 5.0

# All-female single-sex study
--phenotype_file Data/females_only.csv \
--finalreport_files 'Data/FinalReport*.txt' \
--study_prefix DOchln_females \
--single_sex female \
--lod_threshold 7.0
```

### Single-Sex Study Parameter (`--single_sex`)

**Background**: r/qtl2 requires a `sex` column in the covariate file to build the cross object correctly, specifically for X chromosome genotype encoding and sex chromosome handling. This requirement exists even when all animals are the same sex — the column is used as cross metadata, not only as a statistical covariate. In a standard mixed-sex study the column is always present. In a single-sex study it may be absent from the phenotype input, or a user may attempt to exclude it via `exclude_covariates`, which previously caused module 03 to fail with a "sex column not found" error.

**What `--single_sex` does** (across four pipeline components):

1. **`nextflow.config`**: Defines the parameter (default `null`). Set to `'female'` or `'male'` in the config file or pass on the command line.

2. **Module 01 (Phenotype Processing)**: If `--single_sex` is set and the `sex` column is absent from the covariate data, the module auto-adds a `Sex` column populated with the appropriate value (`Female` or `Male`). Additionally, `sex`, `ngen`, and `generation` are now permanently protected from `exclude_covariates` removal — the pipeline logs a note if a user attempts to exclude them and retains the columns regardless.

3. **Module 03 (Control File Generation)**: Receives the covariate file with the `sex` column present, so cross object creation succeeds without modification.

4. **Module 06 (QTL Analysis)**: Before calling `model.matrix()`, a constant-column filter explicitly removes any covariates that have fewer than two unique non-missing values (e.g., `Sex` in an all-female study). Passing such a covariate directly to `model.matrix()` would crash R with a contrasts error; the filter prevents this and logs the names of removed covariates. A secondary guard then checks whether the covariate data frame is empty after filtering; if so, `addcovar` is set to `NULL` and `scan1` runs without additive covariates.

**Key behaviors**:
- `sex`, `ngen`, and `generation` are always retained in the covariate file regardless of `exclude_covariates`; attempting to exclude them produces a log note rather than an error
- In a single-sex study, sex remains in the cross object as required metadata but is explicitly filtered out in module 06 before `model.matrix()` is called — it does not affect the statistical model, and the filtering is recorded in the batch log
- Users do not need to set `exclude_covariates = 'Sex'`; sex is handled correctly without it
- Valid values: `'female'`, `'male'`, `'f'`, `'m'` (case-insensitive); `null` (default) disables the feature

**When to use**:
- Your study cohort contains only females or only males
- The input phenotype CSV does not have a `Sex` column (auto-add mode)
- The input CSV has a `Sex` column but you want to ensure it cannot be accidentally excluded by `exclude_covariates`

**Configuration in `nextflow.config`**:
```groovy
params {
    single_sex = 'female'   // or 'male', or null for mixed-sex studies
}
```

**Command-line equivalent**:
```bash
nextflow run main.nf \
  --phenotype_file Data/females_only.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix MyFemaleStudy \
  --single_sex female \
  -profile standard
```

## Citation

If you use QTL2_NF in your research, please cite:

- **r/qtl2**: Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen Ś, Yandell BS, Churchill GA (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics 211:495-502.

- **Nextflow**: Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319.

- **DO_Pipe**: Please acknowledge the DO_Pipe project (https://github.com/exsquire/DO_pipe) for foundational methodologies incorporated into this pipeline.

## Acknowledgments

This pipeline builds upon the foundational work of the Churchill Laboratory's R/qtl2 software, incorporates methodology from the DO_Pipe Diversity Outbred analysis pipeline, and leverages the broader quantitative genetics community's contributions. Special thanks to:

- The r/qtl2 development team for creating robust software for multiparental population analysis
- The DO_Pipe project for providing foundational scripts and methodologies
- The Mouse Genetics community for continued contributions to quantitative genetics analysis tools
- The Nextflow development team for enabling scalable workflow management
- The broader systems genetics community for advancing high-dimensional phenotype analysis methods