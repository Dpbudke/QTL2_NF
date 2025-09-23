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

The pipeline consists of 9 main modules plus an optional quality control module, organized into a sequential workflow, incorporating proven methodologies from the DO_Pipe project with enhanced automation and scalability:

### Module 1: Phenotype Processing (`phenotype_process.nf`)
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

### Module 2: Genotype Processing (`genotype_process.nf`)
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

### Module 3: Control File Generation (`control_cross2.nf`)
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

### Module 4: Cross2 Object Creation (`control_cross2.nf`)
**Purpose**: Generate validated r/qtl2 cross2 objects

**Key Functions**:
- Cross2 object creation using read_cross2()
- Data integrity validation and diagnostics
- X chromosome encoding for DO crosses
- Quality control reporting
- Object serialization for downstream analysis

**Resources**: 2 CPUs, 4GB RAM, 2 hours

### Module 5: Genome Scan Preparation (`scan_perm.nf`)
**Purpose**: Prepare high-resolution genotype probabilities and kinship matrices

**Key Functions**:
- Pseudomarker insertion at 0.5 cM resolution
- Genotype probability calculation (error_prob=0.002)
- LOCO (Leave One Chromosome Out) kinship matrix generation
- Genetic map preparation for high-resolution mapping
- Memory-efficient data structure optimization

**Resources**: 16 CPUs, 256GB RAM, 12 hours

### Module 6: Genome Scanning (`scan_perm.nf`) - NEW SPLIT MODULE
**Purpose**: Perform comprehensive genome-wide QTL scanning

**Key Functions**:
- Phenotype-agnostic scanning of ALL traits simultaneously
- Linear mixed model analysis with LOCO kinship correction
- Covariate handling (Sex, Generation, Batch)
- Smart LOD threshold filtering for computational efficiency
- Preliminary peak identification (LOD ≥ 3.0)

**Resources**: 32 CPUs, 1TB RAM, 24 hours
**Optimization**: Reduced from cores=0 to cores=8 for memory management

### Module 7: Permutation Testing (`scan_perm.nf`)
**Purpose**: Determine statistical significance thresholds through permutation testing

**Key Functions**:
- **Computational Optimization**: Only test phenotypes passing LOD threshold filter
- 1000 permutation rounds using scan1perm()
- **Massive Efficiency Gains**: Typically 70-90% reduction in permutation workload
- Multi-level significance threshold calculation (0.63, 0.10, 0.05, 0.01)
- Empirical false discovery rate control

**Resources**: 48 CPUs, 1TB RAM, 72 hours

### Module 8: Significant QTL Identification (`scan_perm.nf`)
**Purpose**: Apply significance thresholds and generate final QTL lists

**Key Functions**:
- Empirical significance threshold application
- Multi-level QTL classification
- Summary statistics and chromosome distribution analysis
- High-priority QTL identification (LOD ≥ 10)
- Comprehensive reporting and visualization

**Resources**: 8 CPUs, 64GB RAM, 4 hours

### Module 9: QTL Viewer Integration (`qtl_viewer.nf`)
**Purpose**: Create interactive web-based QTL visualization

**Key Functions**:
- Official QTL Viewer RData format conversion
- Portable Docker deployment package creation
- Interactive visualization setup (LOD plots, effect plots, SNP association)
- Database integration and one-command startup
- Cross-platform compatibility

**Resources**: 8 CPUs, 64GB RAM, 2 hours

### Module 10: DO Sample Mixup QC (Optional) (`do_mixup_qc.nf`)
**Purpose**: Optional quality control module for detecting sample mix-ups in Diversity Outbred mouse data

**Key Functions**:
- Sample identity verification using genetic markers and phenotype correlations
- Cis-eQTL analysis for transcript-genotype consistency checking
- Distance-based clustering analysis for sample relationship validation
- Automated detection of potential sample swaps or mislabeling
- Corrected sample mapping generation for data integrity
- Comprehensive HTML reporting with visualizations and recommendations
- Integration with existing QTL results for enhanced validation accuracy

**Methodology**: Based on Broman et al. (2015) G3 methodologies for multiparental population sample validation, this module performs systematic checks for sample integrity using both genetic and molecular phenotype data.

**Usage**: Can be run as standalone quality control or integrated into the main pipeline workflow for comprehensive data validation.

**Resources**: 8 CPUs, 32GB RAM, 3 hours

## Installation and Setup

### Prerequisites

- **Nextflow** ≥ 23.04.0
- **Docker** or **Apptainer/Singularity** (for HPC)
- **SLURM** scheduler (for HPC execution)

### Quick Start

1. **Clone the repository**:
```bash
git clone https://github.com/Dpbudke/QTL2_NF.git
cd QTL2_NF
```

2. **Prepare your data**:
   - Place phenotype CSV file in `Data/` directory
   - Place GeneSeek FinalReport.txt file(s) in `Data/` directory
   - Update SLURM configuration in `nextflow.config` for your HPC system

3. **Run the pipeline**:
```bash
nextflow run main.nf \
  --phenotype_file Data/your_phenotypes.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix YourStudy \
  -profile standard
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

### Resume Feature - NEW

The pipeline now supports resuming from any step using the `main_resume.nf` workflow:

```bash
# Resume from genome scanning (skip all earlier steps)
nextflow run main_resume.nf \
  --resume_from genome_scan \
  --study_prefix DOchln \
  --lod_threshold 5.0

# Resume from permutation testing
nextflow run main_resume.nf \
  --resume_from permutation \
  --study_prefix DOchln \
  --lod_threshold 7.0

# Resume from QTL Viewer setup only
nextflow run main_resume.nf \
  --resume_from qtlviewer \
  --study_prefix DOchln
```

**Resume Options**:
- `phenotype` - Start from phenotype processing
- `genotype` - Start from genotype processing
- `control` - Start from control file generation
- `cross2` - Start from cross2 object creation
- `prepare_scan` - Start from genome scan preparation
- `genome_scan` - Start from genome scanning
- `permutation` - Start from permutation testing
- `significant_qtls` - Start from QTL identification
- `qtlviewer` - Start from QTL Viewer setup

### Test Mode (Development)
```bash
# Positive control: chromosome 2 coat color analysis
nextflow run main.nf \
  --phenotype_file Data/test_phenotypes.csv \
  --finalreport_files 'Data/FinalReport.txt' \
  --study_prefix CoatColorTest \
  --test_mode \
  -profile standard
```

### DO Sample Mixup QC (Standalone)
The pipeline includes an optional quality control module for detecting sample mix-ups:

```bash
# Standalone mixup detection using existing QTL results
nextflow run modules/do_mixup_qc.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --cross2_object results/04_cross2_creation/DOchln_cross2_object.rds \
  --qtl_results results/06_genome_scanning/DOchln_genome_scan_results.rds \
  --study_prefix DOchln_mixup_check \
  --outdir results/quality_control

# Simplified mixup detection without existing QTL results (uses cis-eQTL)
nextflow run modules/do_mixup_qc.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --cross2_object results/04_cross2_creation/DOchln_cross2_object.rds \
  --study_prefix DOchln_simple_mixup \
  --outdir results/quality_control
```

## Input Data Requirements

### Phenotype File Format
The pipeline requires a specially formatted CSV file that follows the DO_Pipe convention with header row indicators. This format supports multi-phenotype analyses including clinical traits, molecular measurements, and high-dimensional data:

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
- **GeneSeek FinalReport.txt**: Standard Illumina genotyping array output
- **Supported arrays**: GigaMUGA, MiniMUGA, MEGA-MUGA
- **File format**: Tab-delimited with standard GeneSeek headers

## Resource Requirements

### Current Working Dataset Specifications
- **Samples**: 382 Diversity Outbred mice
- **Phenotypes**: 13,374 molecular traits
- **Study prefix**: DOchln
- **LOD threshold**: 7.0
- **Platform**: SLURM HPC with containerization

### Memory Optimizations (Updated 2024)
- **GENOME_SCAN**: Increased to 1TB memory allocation
- **PERMUTATION_TEST**: 1TB memory with intelligent filtering
- **Parallelization**: Reduced from cores=0 to cores=8 to prevent OOM errors
- **Smart filtering**: 70-90% reduction in permutation workload through LOD threshold pre-filtering

### HPC Configuration Requirements
```bash
# Minimum resource allocations
PREPARE_GENOME_SCAN: 16 CPUs, 256GB RAM, 12h
GENOME_SCAN:        32 CPUs, 1TB RAM, 24h
PERMUTATION_TEST:   48 CPUs, 1TB RAM, 72h
```

## Output Structure

Results are organized by analysis in the `results/` directory:

```
results/
├── 01_phenotype_processing/          # Validated phenotype and covariate files
├── 02_genotype_processing/           # Chromosome-specific genotype files
├── 03_control_file_generation/       # r/qtl2 control files and metadata
├── 04_cross2_creation/              # Cross2 objects and validation reports
├── 05_genome_scan_preparation/       # Genotype probabilities and kinship matrices
├── 06_genome_scanning/              # Genome scan results and preliminary peaks
├── 07_permutation_testing/          # Permutation results and significance thresholds
├── 08_significant_qtls/             # Final significant QTL identification
├── 09_qtlviewer_data/              # Interactive QTL Viewer deployment package
└── pipeline_info/                   # Execution reports and monitoring
```

### QTL Viewer Deployment

Module 9 creates a **fully portable QTL Viewer deployment** accessible at `http://localhost:8000`:

```bash
# On HPC where pipeline completed
cd results/09_qtlviewer_data/
./start_qtlviewer.sh

# Download to local machine
scp -r user@hpc:path/to/results/09_qtlviewer_data/ ~/my_qtl_study/
cd ~/my_qtl_study/09_qtlviewer_data/
./start_qtlviewer.sh
```

**Features**:
- Interactive LOD score visualization
- Founder strain effect plots
- High-resolution SNP association mapping
- Cross-phenotype correlation analysis
- Mediation analysis capabilities

## Configuration

### SLURM Configuration
Update `nextflow.config` for your HPC environment:

```groovy
process {
    executor = 'slurm'
    clusterOptions = '--account=your_project_account'

    // Adjust resource allocations based on your cluster
    withName: GENOME_SCAN {
        cpus   = 32
        memory = '1 TB'
        time   = '24h'
    }
}
```

### Pipeline Parameters

**Core Parameters**:
- `--phenotype_file`: Master phenotype CSV file path
- `--finalreport_files`: GeneSeek FinalReport file(s) (supports glob patterns)
- `--study_prefix`: Unique identifier for output files
- `--lod_threshold`: Minimum LOD score for permutation testing (default: 7.0)

**Advanced Parameters**:
- `--sample_filter`: JSON filter for sample subsetting by covariates
- `--test_mode`: Enable chromosome 2 analysis only (development)
- `--outdir`: Output directory (default: results)

**Resume Parameters**:
- `--resume_from`: Resume from specific pipeline step

## Performance Optimization

### LOD Threshold Strategy
The `--lod_threshold` parameter provides significant computational savings:

- **Higher values (8.0-10.0)**: Faster execution, fewer phenotypes tested
- **Lower values (5.0-7.0)**: More comprehensive analysis, longer runtime
- **Recommended**: 7.0 for balanced efficiency and sensitivity

### Memory Management
- **Large datasets**: Use 1TB memory allocation for scanning steps
- **Reduced parallelization**: cores=8 prevents out-of-memory errors
- **Smart filtering**: Pre-filter phenotypes before expensive permutation testing


## DO_Pipe Integration and Acknowledgments

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