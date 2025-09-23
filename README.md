# QTL2_NF Pipeline

A Nextflow pipeline for multiparental mouse QTL analysis using r/qtl2, designed for high-throughput quantitative trait loci mapping in Diversity Outbred (DO) and other multiparental populations. This pipeline supports comprehensive multi-phenotype analyses including clinical phenotype QTLs, expression QTLs (eQTLs), metabolite QTLs (mQTLs), and other high-dimensional molecular trait mapping applications.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-brightgreen.svg)](https://www.docker.com/)

## Overview

QTL2_NF is a comprehensive bioinformatics pipeline that processes phenotype and genotype data to identify quantitative trait loci (QTLs) with rigorous statistical significance testing. The pipeline is optimized for large-scale datasets and provides interactive visualization through QTL Viewer integration. Building upon methodologies from the DO_Pipe project (https://github.com/exsquire/DO_pipe), this pipeline implements automated workflows for multiparental population analysis with enhanced scalability and reproducibility.

### Multi-Phenotype Applications

This pipeline is specifically designed to handle diverse phenotype types commonly encountered in systems genetics:

- **Clinical Phenotype QTLs**: Traditional quantitative traits such as body weight, organ weights, biochemical measurements, and physiological parameters
- **Expression QTLs (eQTLs)**: RNA-seq derived gene expression levels with built-in normalization workflows including filtering, CPM normalization, and inverse rank normalization
- **Metabolite QTLs (mQTLs)**: Metabolomics data from mass spectrometry or NMR platforms
- **Protein QTLs (pQTLs)**: Proteomics measurements from various analytical platforms
- **Epigenetic QTLs**: DNA methylation, histone modifications, and chromatin accessibility data
- **Behavioral and Cognitive Traits**: Neurobehavioral assessments and cognitive performance metrics
- **High-Dimensional Molecular Data**: Any quantitative molecular phenotype suitable for linear modeling

The pipeline automatically handles the preprocessing requirements for molecular data types, including appropriate filtering thresholds, normalization strategies, and transformation procedures to ensure optimal QTL detection power.

## Pipeline Architecture

The pipeline consists of 9 numbered modules organized into a sequential workflow that matches the results directory structure. This modular design enables precise resume capabilities and optimal resource allocation for large-scale QTL analysis:

### Module 1: Phenotype Processing (`01_phenotype_process.nf`)
**Purpose**: Comprehensive validation and preparation of phenotype/covariate data according to r/qtl2 specifications with robust support for diverse molecular data types and experimental designs

**Key Functions**:
- **Data Structure Validation**: Automated parsing of DO_Pipe-formatted CSV files with header row indicators for seamless data type recognition
- **Sample ID Standardization**: Cross-validation of sample identifiers between phenotype and genotype datasets with automatic format harmonization
- **Numeric Matrix Enforcement**: Conversion of phenotype data to numeric matrices with intelligent handling of missing data patterns and outlier detection
- **Required Covariate Validation**: Verification of essential covariates (Sex, generation) as mandated by r/qtl2 for Diversity Outbred cross analysis
- **Optional Covariate Integration**: Flexible incorporation of experimental factors (Diet, batch, age, treatment groups) with automatic categorical encoding
- **Multi-Phenotype Support**: Scalable validation for high-dimensional datasets including thousands of molecular traits (transcripts, metabolites, proteins)
- **Quality Control Diagnostics**: Comprehensive outlier detection, batch effect assessment, and correlation analysis with detailed visualization reports
- **Format Compliance**: Generation of r/qtl2-compliant data structures following DO_Pipe formatting conventions with proper metadata integration
- **Missing Data Management**: Intelligent handling of missing values with appropriate encoding strategies for downstream statistical analysis

**Input Format Requirements**:
The module processes specially formatted CSV files with the DO_Pipe convention:
- **Row 1**: Type indicators distinguishing covariate and phenotype columns
- **Row 2**: Column headers with standardized naming conventions
- **Data Rows**: Sample-wise measurements with automatic data type validation
- **Supported Scales**: Handles diverse phenotype scales from clinical measurements to high-dimensional molecular data

**Technical Implementation**: Uses advanced R data processing libraries (dplyr, data.table) with comprehensive error checking and logging for reliable data preparation in high-throughput environments.

**Resources**: 1 CPU, 2GB RAM, 1 hour

### Module 2: Genotype Processing (`02_genotype_process.nf`)
**Purpose**: Comprehensive conversion of GeneSeek FinalReport files to r/qtl2-compatible format using enhanced DO_Pipe-derived methodologies with robust quality control and validation

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
**Purpose**: Create r/qtl2 control files and metadata structures following DO_Pipe specifications

**Key Functions**:
- JSON control file generation with complete metadata for multiparental crosses
- File path specification and validation for execution environment compatibility
- Genotype encoding mappings (A/H/B) and sex chromosome handling for DO populations
- Cross type specification with founder strain definitions (8 founder strains for DO)
- Reference file integration including genetic maps, founder genotypes, and marker annotations
- Cross information structure creation compatible with r/qtl2 requirements
- Integration of DO_Pipe cross2 generation methodologies for robust object creation

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

### Module 6: HPC Array-Based QTL Analysis (`06_qtl_analysis.nf`)
**Purpose**: Revolutionary HPC array-based genome scanning with intelligent chunking for massive datasets

**Key Functions**:
- **HPC Array Processing**: Coordinator + chunked scanning + results combination architecture
- **Intelligent Chunking**: Processes phenotypes in chunks of ~500 for optimal parallelization
- **Coordinator Process**: Creates array job configuration and chunk assignments (2 CPUs, 8GB, 1h)
- **Chunk Processing**: Individual SLURM array jobs for parallel phenotype analysis (16 CPUs, 64GB, 4h each)
- **Results Integration**: Combines chunk outputs into comprehensive genome scan results (4 CPUs, 32GB, 2h)
- **Smart LOD Filtering**: Pre-filters phenotypes before expensive permutation testing
- **Preliminary Peak Detection**: Identifies QTL peaks with LOD ≥ 3.0 for downstream analysis

**Technical Innovation**: Enables processing of datasets with 10,000+ phenotypes through massive parallelization while preventing memory bottlenecks

**Resources**: Coordinator (2 CPUs, 8GB), Per-chunk (16 CPUs, 64GB), Combination (4 CPUs, 32GB)

### Module 7: Permutation Testing (`07_permutation_testing.nf`)
**Purpose**: Computationally optimized permutation testing with intelligent pre-filtering

**Key Functions**:
- **Smart Pre-filtering**: Only processes phenotypes that pass LOD threshold from Module 6
- **Massive Efficiency Gains**: Typically 70-90% reduction in computational workload
- **Comprehensive Testing**: 1000 permutation rounds using scan1perm() for robust statistics
- **Multi-level Thresholds**: Calculates significance thresholds at 0.63, 0.10, 0.05, 0.01 levels
- **Empirical FDR Control**: False discovery rate calculation for multiple testing correction
- **Resource Optimization**: Maximum resource allocation for intensive statistical computation

**Resources**: 48 CPUs, 1TB RAM, 72 hours

### Module 8: Significant QTL Identification (`08_identify_significant_qtls.nf`)
**Purpose**: Final QTL identification and comprehensive summary generation

**Key Functions**:
- **Threshold Application**: Applies empirical significance thresholds from permutation testing
- **Multi-level Classification**: Categorizes QTLs by significance levels (0.63, 0.10, 0.05, 0.01)
- **Comprehensive Summaries**: Generates detailed QTL lists with genomic positions and effect sizes
- **High-priority Identification**: Flags exceptional QTLs with LOD ≥ 10 for priority follow-up
- **Statistical Reporting**: Summary statistics, chromosome distribution, and validation metrics
- **Publication-ready Output**: Formatted tables and visualizations for scientific reporting

**Resources**: 8 CPUs, 64GB RAM, 4 hours

### Module 9: QTL Viewer Integration (`09_qtl_viewer.nf`)
**Purpose**: Interactive QTL visualization and portable deployment system

**Key Functions**:
- **QTL Viewer Integration**: Converts results to official QTL Viewer RData format
- **Portable Deployment**: Creates complete Docker deployment package with one-command startup
- **Interactive Visualization**: LOD score plots, founder effect plots, SNP association mapping
- **Cross-platform Compatibility**: Deployable on HPC, local machines, or cloud platforms
- **Database Integration**: Self-contained database for rapid phenotype and QTL querying
- **User-friendly Interface**: Web-based interface accessible via http://localhost:8000

**Resources**: 8 CPUs, 64GB RAM, 2 hours

### Optional: DO Sample Mixup QC (`do_mixup_qc.nf`)
**Purpose**: Standalone quality control module for detecting sample mix-ups in Diversity Outbred populations

**Key Functions**:
- **Sample Identity Verification**: Cross-validation using genetic markers and phenotype correlations
- **Cis-eQTL Analysis**: Transcript-genotype consistency checking for molecular phenotypes
- **Distance-based Clustering**: Sample relationship validation using genetic similarity
- **Automated Detection**: Identifies potential sample swaps, mislabeling, or contamination
- **Corrective Mapping**: Generates corrected sample assignments for data integrity
- **Comprehensive Reporting**: HTML reports with visualizations and actionable recommendations
- **Integration Ready**: Compatible with existing QTL results for enhanced validation

**Scientific Methodology**: Implements Broman et al. (2015) G3 validation approaches for multiparental populations

**Usage**: Deployable as standalone QC or integrated into main pipeline workflow

**Resources**: 8 CPUs, 32GB RAM, 3 hours

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
# - Phenotype CSV (DO_Pipe format) → Data/
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

# Clinical phenotype QTL analysis with lower threshold
nextflow run main.nf \
  --phenotype_file Data/clinical_traits.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_clinical \
  --lod_threshold 5.0 \
  --outdir results/clinical_QTL \
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
nextflow run main_resume.nf \
  --resume_from 09 \
  --study_prefix DOchln
```

**Resume Options (Module Numbers)**:
- `01` or `phenotype` - Start from phenotype processing and validation
- `02` or `genotype` - Start from genotype processing and chromosome file generation
- `03` or `control` - Start from control file generation and metadata creation
- `04` or `cross2` - Start from cross2 object creation and validation
- `05` or `prepare_scan` - Start from genome scan preparation (genotype probabilities)
- `06` or `qtl_analysis` - Start from HPC array genome scanning and peak detection
- `07` or `permutation` - Start from permutation testing and significance thresholds
- `08` or `significant_qtls` - Start from significant QTL identification and summarization
- `09` or `qtlviewer` - Start from QTL Viewer data preparation and deployment

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

**Local Development**:
```bash
# Local execution with Docker
nextflow run main.nf \
  --phenotype_file Data/small_dataset.csv \
  --finalreport_files 'Data/test_FinalReport.txt' \
  --study_prefix DevTest \
  -profile docker
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

## Input Data Requirements

### Phenotype File Format
The pipeline requires a specially formatted CSV file following the DO_Pipe convention with header row indicators. This standardized format enables automated parsing and supports diverse phenotype types from clinical measurements to high-dimensional molecular data:

**File Structure Requirements**:
The pipeline uses a specialized CSV format developed by the DO_Pipe project that facilitates automated parsing of metadata and phenotype information:

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

## Resource Requirements and Technical Specifications

### Current Working Dataset (Example Scale)
- **Study Population**: 382 Diversity Outbred mice
- **Molecular Phenotypes**: 13,374 transcriptome traits (eQTL analysis)
- **Genotyping Platform**: GigaMUGA array (65,000+ markers)
- **Study Identifier**: DOchln (Diversity Outbred Choline study)
- **Processing Strategy**: 27 chunks of ~500 phenotypes for optimal parallelization
- **LOD Threshold**: 7.0 (balanced efficiency and sensitivity)
- **Container Platform**: dpbudke/qtl2-pipeline:latest (unified containerization)

### HPC Platform Requirements
- **Scheduler**: SLURM workload manager with array job support
- **Containerization**: Apptainer/Singularity for HPC environments
- **Storage**: High-performance filesystem for intensive I/O operations
- **Network**: High-bandwidth interconnect for parallel processing phases

### Resource Optimization Strategy
The pipeline implements intelligent resource allocation optimized for each computational phase:

- **Module 1 (Phenotype)**: 1 CPU, 2GB RAM, 1h - Lightweight data validation
- **Module 2 (Genotype)**: 2 CPUs, 64GB RAM, 4h - Memory-intensive genotype processing
- **Module 3 (Control Files)**: 2 CPUs, 4GB RAM, 2h - Metadata and configuration generation
- **Module 4 (Cross2 Creation)**: 2 CPUs, 4GB RAM, 2h - r/qtl2 object creation
- **Module 5 (Scan Preparation)**: 32 CPUs, 128GB RAM, 8h - High-performance probability calculation
- **Module 6 (QTL Analysis)**: HPC array architecture with three-tier resource allocation
- **Module 7 (Permutation)**: 48 CPUs, 1TB RAM, 72h - Maximum resources for statistical validation
- **Module 8 (QTL Identification)**: 8 CPUs, 64GB RAM, 4h - Results processing and summary
- **Module 9 (QTL Viewer)**: 8 CPUs, 64GB RAM, 2h - Visualization data preparation

### HPC Array Processing Architecture (Module 6)
```bash
# Three-tier resource allocation for massive scalability
Coordinator Job:    2 CPUs,   8GB RAM,  1h    # Array setup and configuration
Per-chunk Job:     16 CPUs,  64GB RAM,  4h    # Individual phenotype chunks (parallel)
Combination Job:    4 CPUs,  32GB RAM,  2h    # Results integration and validation

# Example: 13,374 phenotypes processed as 27 chunks of ~500 phenotypes each
# Total parallel processing: 27 simultaneous SLURM array jobs
```

### Performance Achievements
- **Modular Architecture**: Independent optimization of each processing phase
- **Bottleneck Resolution**: Module 5 resource enhancement eliminated major pipeline constraints
- **Massive Scalability**: HPC array processing enables 10,000+ phenotype datasets
- **Intelligent Pre-filtering**: 70-90% reduction in permutation computational workload
- **Memory Management**: Chunked processing prevents out-of-memory failures
- **Resource Scaling**: Precise resource allocation matched to computational complexity

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
├── 06_qtl_analysis/                 # HPC array genome scan results
│   ├── chunks/                      # Individual chunk processing results
│   │   ├── chunk_1_results.rds      # Scan results for phenotype chunk 1
│   │   ├── chunk_2_results.rds      # Scan results for phenotype chunk 2
│   │   └── ...                      # Additional chunk results
│   ├── chunk_info.txt              # Array job configuration and phenotype assignments
│   ├── {prefix}_scan_results.rds    # Combined comprehensive scan results
│   ├── {prefix}_scan_peaks.csv      # Preliminary peaks (LOD ≥ 3.0)
│   └── genome_scan_validation_report.txt
├── 07_permutation_testing/          # Statistical significance determination
│   ├── {prefix}_permutation_results.rds  # Complete permutation test results
│   ├── {prefix}_significance_thresholds.csv  # Multi-level significance thresholds
│   ├── filtered_phenotypes.txt      # Phenotypes meeting LOD threshold criteria
│   └── permutation_validation_report.txt
├── 08_significant_qtls/             # Final QTL identification and analysis
│   ├── {prefix}_significant_qtls.csv    # Comprehensive significant QTL list
│   ├── {prefix}_qtl_summary.txt     # Summary statistics and distributions
│   ├── {prefix}_high_priority_qtls.csv  # Exceptional QTLs (LOD ≥ 10)
│   └── qtl_identification_report.txt
├── 09_qtl_viewer/                   # Interactive visualization deployment
│   ├── qtlviewer_data.RData        # QTL Viewer compatible data package
│   ├── docker-compose.yml          # Container orchestration configuration
│   ├── start_qtlviewer.sh          # One-command deployment script
│   ├── instructions.txt            # Detailed deployment instructions
│   └── qtlviewer_conversion_report.txt
└── pipeline_info/                   # Execution monitoring and reporting
    ├── execution_timeline.html      # Interactive execution timeline
    ├── execution_report.html        # Comprehensive resource usage report
    ├── execution_trace.txt          # Detailed process execution trace
    └── pipeline_dag.svg             # Visual workflow representation
```

### QTL Viewer Deployment

Module 9 creates a **fully portable QTL Viewer deployment** accessible at `http://localhost:8000`:

```bash
# On HPC where pipeline completed
cd results/09_qtl_viewer/
./start_qtlviewer.sh

# Download to local machine for portable analysis
scp -r user@hpc:path/to/results/09_qtl_viewer/ ~/my_qtl_study/
cd ~/my_qtl_study/09_qtl_viewer/
./start_qtlviewer.sh

# Alternative: Direct Docker deployment
docker-compose up -d
```

**Features**:
- Interactive LOD score visualization
- Founder strain effect plots
- High-resolution SNP association mapping
- Cross-phenotype correlation analysis
- Mediation analysis capabilities

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
- `--phenotype_file`: Path to master phenotype CSV file (DO_Pipe format)
- `--study_prefix`: Unique study identifier for all output files

**Optional Parameters**:
- `--finalreport_files`: GeneSeek FinalReport file(s) with glob pattern support (e.g., 'Data/FinalReport*.txt')
- `--lod_threshold`: LOD score threshold for permutation pre-filtering (default: 7.0)
- `--outdir`: Output directory path (default: results)
- `--sample_filter`: JSON string for sample subsetting (e.g., '{"Sex":["male"],"Diet":["hc"]}')
- `--test_mode`: Development mode - process chromosome 2 only (default: false)

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
```

## Performance Optimization

### LOD Threshold Strategy
Intelligent pre-filtering dramatically reduces computational requirements:

- **Conservative (8.0-10.0)**: Maximum efficiency, fewer phenotypes in permutation testing
- **Balanced (7.0)**: Recommended default - optimal efficiency/sensitivity trade-off
- **Comprehensive (5.0-6.0)**: More inclusive analysis, increased computational time
- **Liberal (3.0-4.0)**: Maximum sensitivity, longest runtime

**Efficiency Impact**: LOD threshold of 7.0 typically reduces permutation workload by 70-90%

### Memory Management Strategy
- **Chunked Processing**: HPC array jobs prevent memory bottlenecks in large datasets
- **Progressive Resource Allocation**: Resources scale with computational complexity
- **Smart Memory Usage**: Module 5 uses 128GB for probability calculations, Module 7 uses 1TB for permutations
- **Container Optimization**: Single container approach reduces overhead and complexity

### Scalability Achievements
- **Phenotype Capacity**: Successfully processes 10,000+ molecular traits
- **Sample Capacity**: Handles 1,000+ samples with full genotype data
- **Parallel Efficiency**: 27 simultaneous array jobs for 13,374 phenotypes
- **Resource Efficiency**: Intelligent filtering reduces computational time by orders of magnitude


## Scientific Methodology and Integration

QTL2_NF represents a comprehensive automation and scaling of methodologies originally developed in the DO_Pipe project (https://github.com/exsquire/DO_pipe), which provides foundational scripts and approaches for Diversity Outbred population analysis. This pipeline transforms the manual DO_Pipe workflow into a fully automated, reproducible, and scalable Nextflow implementation while preserving the scientific rigor and methodological approaches of the original project.

**Key DO_Pipe Integrations and Enhancements**:

**Genotype Processing (Module 2)**:
- Adapts DO_Pipe's comprehensive approach for GeneSeek FinalReport file parsing and validation
- Incorporates DO_Pipe's sample filtering strategies and marker quality control methodologies
- Automates chromosome-specific file generation following DO_Pipe formatting standards
- Enhances original scripts with improved error handling and resource management for HPC environments

**Cross2 Object Creation (Module 3)**:
- Directly implements DO_Pipe's cross2 generation logic with enhanced automation
- Preserves DO_Pipe's metadata structure requirements and founder strain definitions (8 DO founder strains)
- Adds comprehensive validation and diagnostics beyond the original DO_Pipe implementation
- Integrates reference file handling and genetic map specifications from DO_Pipe standards

**Data Format Conventions**:
- Maintains full compatibility with DO_Pipe's phenotype file format using header row indicators
- Preserves DO_Pipe's covariate/phenotype column distinction methodology
- Extends DO_Pipe's sample ID standardization procedures for improved robustness
- Implements DO_Pipe's approach to handling missing data and categorical covariates

**Statistical Methodology**:
- Builds upon DO_Pipe's QTL mapping strategies using r/qtl2 linear mixed models
- Incorporates DO_Pipe's permutation testing framework with enhanced computational efficiency
- Extends DO_Pipe's significance threshold calculation methods with multi-level thresholding
- Preserves DO_Pipe's kinship matrix generation and LOCO (Leave One Chromosome Out) correction approaches

**Quality Control Integration**:
- Implements enhanced versions of DO_Pipe's data validation procedures
- Adds automated sample mixup detection capabilities based on established DO population analysis methods
- Incorporates DO_Pipe's outlier detection and batch effect assessment strategies

**Automation and Scalability Enhancements**:
- Transforms DO_Pipe's manual R scripts into containerized, reproducible Nextflow modules
- Adds intelligent resource management and parallel processing capabilities for large datasets
- Implements resume functionality and checkpoint recovery not available in original DO_Pipe
- Provides automated Docker-based QTL Viewer deployment extending DO_Pipe's visualization capabilities

We gratefully acknowledge the DO_Pipe project and its contributors for providing the foundational scientific methodology and computational framework that enabled the development of this automated pipeline workflow. The DO_Pipe project's rigorous approach to multiparental population analysis forms the scientific backbone of QTL2_NF's implementation.

## Citation

If you use QTL2_NF in your research, please cite:

- **r/qtl2**: Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen Ś, Yandell BS, Churchill GA (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics 211:495-502.

- **Nextflow**: Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319.

- **DO_Pipe**: Please acknowledge the DO_Pipe project (https://github.com/exsquire/DO_pipe) for foundational methodologies incorporated into this pipeline.

## Acknowledgments

This pipeline builds upon the foundational work of the Churchill Laboratory's R/qtl2 software, incorporates methodology from the DO_Pipe Diversity Outbred analysis pipeline, and leverages the broader quantitative genetics community's contributions. Special thanks to:

- The r/qtl2 development team for creating robust software for multiparental population analysis
- The DO_Pipe project contributors for providing foundational scripts and methodologies
- The Mouse Genetics community for continued contributions to quantitative genetics analysis tools
- The Nextflow development team for enabling scalable workflow management
- The broader systems genetics community for advancing high-dimensional phenotype analysis methods