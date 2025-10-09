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

### Module 7: Simplified Batch Permutation Testing (`07_permutation_testing.nf`)
**Purpose**: Efficient single-phenotype permutation testing with fork-based parallelism and intelligent pre-filtering

**Key Functions**:
- **Smart Pre-filtering**: Only processes phenotypes passing LOD threshold from Module 6 or regional filtering from Module 6b
- **Simplified Batch Architecture**: 1 batch = 1 phenotype × 50 permutations = 1 SLURM job (20 batches per phenotype = 1,000 total permutations)
- **Fork-Based Parallelism**: Uses scan1perm() with cores=48 and fork (NOT PSOCK) for efficient memory sharing on Linux
- **Memory Optimization**: Alleleprob-based (2.1GB) instead of genoprob (9.3GB) for 4.4x memory reduction
- **Real-Time Progress**: Individual batch completion enables live monitoring of permutation progress
- **Automatic Resume**: Detects existing batch results and skips completed permutation jobs
- **Filtered Cross2 Object**: Creates phenotype-filtered cross2 for downstream Module 8 analysis
- **Multi-level Thresholds**: Calculates significance thresholds at 63%, 90%, 95%, 99% quantiles
- **Comprehensive Aggregation**: Quilts together all batch results into final permutation matrix (1000 perms × N phenotypes)

**Efficiency Achievement**: LOD threshold of 7.0 typically reduces permutation workload by 70-90%; regional filtering can achieve 95-98% reduction

**Resources**: Setup (4 CPUs, 32GB, 10m), Batch Jobs (48 CPUs, 100GB, 1h each), Aggregation (4 CPUs, 32GB, 1h)

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

### Module 9: QTL Viewer Data Preparation (`09_qtl_viewer.nf`)
**Purpose**: Automated conversion to QTL Viewer RData format with complete Docker deployment package for local visualization

**Key Functions**:
- **Intelligent Cross2 Detection**: Auto-detects and prioritizes regionally-filtered cross2 (Module 7) over full cross2 (Module 4)
- **QTL Viewer Format Conversion**: Transforms r/qtl2 objects into official QTL Viewer RData structure
- **Required Elements Assembly**: Creates ensembl.version, genoprobs, K (kinship), map, markers, and dataset.phenotypes objects
- **Optimized Allele Probabilities**: Uses alleleprob (2.1GB) instead of genoprob (9.3GB) for 4.4x size reduction while maintaining all necessary information for visualization
- **Phenotype Annotations**: Generates complete phenotype metadata with data types and categories
- **Sample Annotations**: Integrates covariate information and sample identifiers
- **Covariate Matrix**: Prepares model matrix for interactive QTL Viewer covariate selection
- **LOD Peaks Integration**: Incorporates significant QTLs from Module 8 for interactive exploration
- **Docker Deployment Package**: Creates docker-compose.yml, startup script, README, and data directory structure
- **MGI Database Integration**: Automatically sources SQLite database for founder SNP annotation (tries local cache, then Figshare download)
- **Portable Deployment**: Complete self-contained package downloadable to local machine for visualization

**Deployment Features**:
- **One-Command Startup**: `./start_qtlviewer.sh` launches both qtl2rest (API) and qtl2web (interface)
- **Dual-Container Architecture**: qtl2rest (port 8001) + qtl2web (port 8000) with bridge networking
- **Interactive Web Interface**: http://localhost:8000 for LOD plots, effect plots, and SNP associations
- **REST API Access**: http://localhost:8001 for programmatic data queries
- **Cross-Platform**: Works on any system with Docker Desktop (Windows, Mac, Linux)

**Resources**: PREPARE_QTLVIEWER_DATA (16 CPUs, 256GB, 6h), SETUP_QTLVIEWER_DEPLOYMENT (2 CPUs, 16GB, 1h)

### Optional: Broman Mixup QC (`broman_mixup_qc.nf`)
**Purpose**: Detect sample mix-ups by comparing observed expression patterns with genotype-predicted expression using eQTL signatures

**Scientific Methodology**:
Based on Broman et al. (2015) G3 and Westra et al. (2011), this module identifies sample mislabeling by testing whether observed gene expression matches the expression pattern predicted from an individual's genotype at strong eQTL positions.

**Process Workflow**:
1. **eQTL Selection**: Identifies top N eQTL (default: 100) with highest LOD scores from Module 8 significant QTLs, retaining only one peak per unique gene to avoid redundancy
2. **Observed Expression Extraction**: Retrieves actual measured expression values for genes with strong eQTL
3. **Predicted Expression Calculation**:
   - Extracts genotype probabilities at eQTL marker positions using `qtl2::pull_genoprobpos()`
   - Fits single-QTL models using `qtl2::fit1()` with `zerosum=FALSE` to preserve original scale
   - Calculates predicted expression from genotype probabilities: `predicted = genoprob %*% coefficients`
4. **Distance Matrix Computation**: Uses `lineup2::dist_betw_matrices()` to calculate RMS distances between all pairs of observed and predicted expression profiles
5. **Mixup Detection**:
   - Extracts self-distances using `lineup2::get_self()` (observed vs own predicted)
   - Identifies minimum distances using `lineup2::get_best()` (observed vs any predicted)
   - Flags samples where self-distance exceeds minimum distance to any other sample's predicted expression

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
# Clone repository
git clone https://github.com/Dpbudke/QTL2_NF.git
cd QTL2_NF

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
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_eQTL \
  --lod_threshold 7.0 \
  --outdir results/eQTL_analysis \
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
  --resume_from 05 \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from HPC array genome scanning (most common restart point)
nextflow run main_resume.nf \
  --resume_from 06 \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from permutation testing (after scan completion)
nextflow run main_resume.nf \
  --resume_from 07 \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from QTL Viewer setup only (visualization-only run)
# Note: Module 9 should be run in background for long processing times
nohup nextflow run main_resume.nf \
  --resume_from qtlviewer \
  --study_prefix DOchln > nextflow_run.log 2>&1 &
```

**Background Execution for Long-Running Modules**:
For Module 9 and other pipeline executions exceeding 2 minutes, use background execution to avoid timeout issues:
```bash
nohup nextflow run main_resume.nf \
  --resume_from qtlviewer \
  --study_prefix DOchln > nextflow_run.log 2>&1 &

# Monitor progress
tail -f nextflow_run.log
```

**Resume Options (Module Numbers)**:
-  `phenotype` - Start from phenotype processing and validation
-  `genotype` - Start from genotype processing and chromosome file generation
-  `control` - Start from control file generation and metadata creation
-  `cross2` - Start from cross2 object creation and validation
-  `prepare_scan` - Start from genome scan preparation (genotype probabilities)
-  `qtl_analysis` - Start from HPC array genome scanning and peak detection
-  `permutation` - Start from permutation testing and significance thresholds
-  `significant_qtls` - Start from significant QTL identification and summarization
-  `qtlviewer` - Start from QTL Viewer data preparation and deployment

**Perfect Module Alignment**: Resume capability matches the numbered module structure exactly, enabling intuitive restart from any processing step with full compatibility.

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
├── 07_permutation_testing/          # Simplified batch permutation testing
│   ├── batches/                     # Individual batch permutation results
│   │   ├── {prefix}_phenotype1_batch_1.rds   # Permutation batch 1 for phenotype 1
│   │   ├── {prefix}_phenotype1_batch_2.rds   # Permutation batch 2 for phenotype 1
│   │   └── ...                               # Additional batch results (20 batches × N phenotypes)
│   ├── {prefix}_filtered_cross2.rds         # Phenotype-filtered cross2 object
│   ├── {prefix}_phenotype_list.txt          # List of phenotypes for permutation
│   ├── {prefix}_setup_log.txt               # Permutation setup log
│   ├── {prefix}_fullPerm.rds                # Complete permutation matrix (1000 × N phenotypes)
│   ├── {prefix}_permThresh.rds              # Multi-level significance thresholds (63%, 90%, 95%, 99%)
│   ├── {prefix}_filtered_phenotypes_lod{threshold}.txt  # Phenotypes passing final threshold
│   └── {prefix}_aggregation_log.txt         # Permutation aggregation log
├── 08_significant_qtls/             # Final QTL identification and analysis
│   ├── {prefix}_significant_qtls.csv    # Comprehensive significant QTL list
│   ├── {prefix}_qtl_summary.txt     # Summary statistics and distributions
│   ├── {prefix}_high_priority_qtls.csv  # Exceptional QTLs (LOD ≥ 10)
│   └── qtl_identification_report.txt
├── 09_qtl_viewer_data/              # QTL Viewer deployment package
│   ├── {prefix}_qtlviewer.RData    # QTL Viewer compatible RData file
│   ├── qtlviewer_conversion_report.txt  # Data conversion validation report
│   ├── docker-compose.yml          # Docker container orchestration
│   ├── start_qtlviewer.sh          # One-command deployment script (executable)
│   ├── README_qtlviewer.md         # Comprehensive deployment instructions
│   └── data/                       # QTL Viewer data directory
│       ├── rdata/                  # RData files for QTL Viewer
│       │   └── {prefix}.RData      # Study-specific QTL data
│       └── sqlite/                 # SQLite databases
│           └── ccfoundersnps.sqlite  # MGI mouse gene annotations and founder SNPs
└── pipeline_info/                   # Execution monitoring and reporting
    ├── execution_timeline.html      # Interactive execution timeline
    ├── execution_report.html        # Comprehensive resource usage report
    ├── execution_trace.txt          # Detailed process execution trace
    └── pipeline_dag.svg             # Visual workflow representation
```

### QTL Viewer Deployment

Module 9 creates a **fully portable QTL Viewer deployment package** with auto-detection of filtered vs. full cross2 objects:

```bash
# Enable QTL Viewer during pipeline run (optional - designed for local deployment)
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln \
  --run_qtlviewer \
  -profile standard

# Download entire deployment directory to local machine
scp -r user@hpc:path/to/results/09_qtl_viewer_data/ ~/my_qtl_study/
cd ~/my_qtl_study/09_qtl_viewer_data/

# One-command local deployment
./start_qtlviewer.sh

# Alternative: Direct Docker Compose deployment
docker-compose up -d
```

**Intelligent Cross2 Detection**:
- **Prioritizes filtered cross2**: Uses `{prefix}_filtered_cross2.rds` from Module 7 if regional filtering was applied
- **Falls back to full cross2**: Uses `{prefix}_cross2.rds` from Module 4 if no filtering was performed
- **Automatic selection**: No manual configuration required - pipeline detects and uses appropriate dataset

**QTL Viewer Features**:
- **Interactive LOD Score Visualization**: Chromosome-wide LOD plots with zoom and pan capabilities
- **Founder Strain Effect Plots**: Visualize allelic effects for each DO founder strain
- **High-Resolution SNP Association**: Fine-mapping with MGI founder SNP database integration
- **Phenotype Browser**: Search and filter phenotypes with interactive data tables
- **Covariate Selection**: Dynamic covariate adjustment through web interface
- **REST API Access**: Programmatic data queries at http://localhost:8001

**Deployment Package Contents**:
- **{prefix}_qtlviewer.RData**: Complete QTL Viewer data bundle with genoprobs, kinship, markers, and results
- **docker-compose.yml**: Pre-configured dual-container setup (qtl2rest + qtl2web)
- **start_qtlviewer.sh**: Automated startup script with health checks and status reporting
- **README_qtlviewer.md**: Comprehensive local deployment instructions
- **data/ directory**: Organized RData and SQLite database files for QTL Viewer containers

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
- `--test_mode`: Development mode - process chromosome 2 only (default: false)

**Regional Filtering Parameters (Module 6b)**:
- `--qtl_region`: Genomic region(s) for targeted QTL analysis (e.g., 'chr12:81-91' or 'chr12:81-91;chr2:50-60')
- `--gtf_file`: Path to Ensembl GTF annotation for gene location mapping (default: 'Data/Mus_musculus.GRCm39.105.gtf.gz')

**Resume Parameters**:
- `--resume_from`: Resume from specific module (01-09, or module names)

**Quality Control Parameters**:
- `--run_mixup_qc`: Enable standalone sample mixup detection (default: false)
- `--run_qtlviewer`: Enable QTL Viewer data preparation and deployment package creation (default: false)

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