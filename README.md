# QTL2_NF Pipeline

A Nextflow pipeline for multiparental mouse QTL analysis using r/qtl2, designed for high-throughput quantitative trait loci mapping in Diversity Outbred (DO) and other multiparental populations.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-brightgreen.svg)](https://www.docker.com/)

## Overview

QTL2_NF is a comprehensive bioinformatics pipeline that processes phenotype and genotype data to identify quantitative trait loci (QTLs) with rigorous statistical significance testing. The pipeline is optimized for large-scale datasets and provides interactive visualization through QTL Viewer integration.

### Key Features

- **Scalable Architecture**: Handles large datasets (382+ samples, 13,000+ phenotypes)
- **Memory Optimized**: Up to 1TB memory allocation for genome scanning
- **Resume Capability**: Start from any pipeline step using existing outputs
- **Smart Filtering**: LOD threshold-based permutation optimization
- **Interactive Visualization**: Integrated QTL Viewer deployment
- **HPC Ready**: SLURM executor with containerization support

## Pipeline Architecture

The pipeline consists of 9 main modules organized into a sequential workflow:

### Module 1: Phenotype Processing (`phenotype_process.nf`)
**Purpose**: Validate and prepare phenotype/covariate data according to r/qtl2 specifications

**Key Functions**:
- Sample ID format validation and consistency checking
- Numeric phenotype matrix enforcement
- Required covariate validation (sex, generation)
- Outlier detection and batch effect diagnostics
- r/qtl2-compliant CSV structure generation

**Resources**: 1 CPU, 2GB RAM, 1 hour

### Module 2: Genotype Processing (`genotype_process.nf`)
**Purpose**: Convert GeneSeek FinalReport to r/qtl2 format

**Key Functions**:
- GeneSeek FinalReport parsing with batch handling
- Genotype call conversion to A/H/B encoding
- Chromosome-specific file generation (1-19, X, Y, M)
- Sample ID consistency with phenotype data
- Missing data handling and integrity validation

**Resources**: 2 CPUs, 64GB RAM, 4 hours

### Module 3: Control File Generation (`control_cross2.nf`)
**Purpose**: Create r/qtl2 control files and metadata structures

**Key Functions**:
- JSON control file generation with complete metadata
- File path specification for execution environment
- Genotype and sex code mappings
- Cross information and founder strain definitions
- Reference file integration (genetic maps, founder genotypes)

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
```bash
# Complete analysis with standard parameters
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln \
  --lod_threshold 7.0 \
  -profile standard

# With sample filtering (males on high-cholesterol diet)
nextflow run main.nf \
  --phenotype_file Data/QTL2_NF_meta_pheno_input.csv \
  --finalreport_files 'Data/FinalReport*.txt' \
  --study_prefix DOchln_males_hc \
  --lod_threshold 7.0 \
  --sample_filter '{"Sex":["male"],"Diet":["hc"]}' \
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

## Input Data Requirements

### Phenotype File Format
CSV file with specific structure requirements:
- **Sample ID column**: Must match genotype sample IDs exactly
- **Numeric phenotype columns**: All phenotype data must be numeric
- **Required covariates**: Sex (male/female), generation number
- **Optional covariates**: Diet, batch, age, etc.

Example structure:
```csv
sample_id,Sex,generation,coat_color,weight_10wk,liver_weight
DO1001,male,1,0.85,25.4,1.2
DO1002,female,1,0.12,22.1,0.9
```

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

## Troubleshooting

### Common Issues

**Memory Errors**:
```bash
# Increase memory allocation in nextflow.config
withName: GENOME_SCAN {
    memory = '1 TB'
}
```

**Sample ID Mismatches**:
- Verify phenotype and genotype sample IDs match exactly
- Check for leading/trailing whitespace or special characters
- Review Module 1 validation report for specific issues

**Resume Failures**:
```bash
# Verify required files exist for resume step
ls results/05_genome_scan_preparation/  # For --resume_from genome_scan
ls results/06_genome_scanning/          # For --resume_from permutation
```

**SLURM Job Failures**:
- Check cluster resource availability
- Verify account permissions and allocations
- Review SLURM output logs in `work/` directory

### Getting Help

- **Pipeline issues**: Check execution reports in `results/pipeline_info/`
- **Resource problems**: Review SLURM logs and resource monitoring
- **Data format errors**: Examine Module 1 validation reports
- **Performance optimization**: Adjust `--lod_threshold` based on dataset size

## Citation

If you use QTL2_NF in your research, please cite:

- **r/qtl2**: Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen Ś, Yandell BS, Churchill GA (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics 211:495-502.

- **Nextflow**: Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319.

## Acknowledgments

This pipeline builds upon the foundational work of the Churchill Laboratory's R/qtl2 software and incorporates methodology from the Diversity Outbred analysis pipeline. Special thanks to the r/qtl2 development team and the Mouse Genetics community for their continued contributions to quantitative genetics analysis tools.

## License

This project is licensed under the MIT License - see the LICENSE file for details.