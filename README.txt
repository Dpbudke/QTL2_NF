# QTL2_NF Pipeline Overview

## Module 1: Data Validation & Phenotype Processing
**Purpose:** Prepare and validate phenotype/covariate data according to r/qtl2 specifications

**Inputs:** Master phenotype CSV, study prefix

**Outputs:** Validated covar.csv, pheno.csv, diagnostic plots

**Custom Scripts Used:** robustZmat.R, covCheck.R

**Key Responsibilities:**
- Validate sample ID format consistency (critical for downstream matching)
- Enforce strict numeric phenotype matrix requirement
- Ensure required covariates exist (sex, ngen) with proper encodings
- Handle both traditional phenotypes and expression data formats
- Generate diagnostic plots for outlier detection and batch effects
- Create r/qtl2-compliant CSV structures

## Module 2: Genotype File Processing
**Purpose:** Convert GeneSeek FinalReport to r/qtl2 format

**Inputs:** FinalReport.txt, allele codes CSV, study prefix

**Outputs:** Chromosome-specific genotype CSVs (prefix_geno1.csv, etc.)

**Custom Scripts Used:** Modified version of DO_pipe genotype processing logic

**Key Responsibilities:**
- Parse GeneSeek FinalReport.txt (handle variable headers, batch inconsistencies)
- Apply sample ID transformations consistently with Module 1
- Convert genotype calls to A/H/B encoding using founder allele codes
- Generate chromosome-specific files (1-19, X, Y, M)
- Handle missing data properly (encode as - not NA)
- Validate genotype encoding integrity

## Module 3: Control File Generation & Cross2 Creation
**Purpose:** Create r/qtl2 control file and generate cross2 object

**Inputs:** All processed CSV files from Modules 1-2, MUGA reference files

**Outputs:** Control JSON file, cross2.rds object, validation reports

**Custom Scripts Used:** Portions of DO_pipe cross2 generation logic

**Key Responsibilities:**
- Generate properly formatted control file (JSON/YAML) with all metadata
- Specify file paths relative to execution environment
- Define genotype codes, sex codes, cross information mappings
- Create cross2 object using read_cross2()
- Validate cross2 object integrity and report diagnostics
- Handle X chromosome encoding complexities for DO crosses

## Module 4: Genome Scanning & Permutations
**Purpose:** Perform QTL/eQTL mapping with statistical significance testing

**Inputs:** cross2 object, kinship matrices, covariate selections

**Outputs:** Scan results, permutation thresholds, intermediate files

**Custom Scripts Used:** Simplified versions of DO_pipe scanning logic without complex SLURM management

**Key Responsibilities:**
- Calculate genotype probabilities (calc_genoprob())
- Generate kinship matrices (calc_kinship())
- Perform genome scans (scan1()) with proper covariate handling
- Execute permutation testing (scan1perm()) for significance thresholds
- Handle both autosomal and X chromosome analysis appropriately
- Manage memory efficiently for large-scale eQTL analysis

## Module 5: Results Processing & QTL Viewer Integration
**Purpose:** Extract significant QTL, generate visualizations, and create QTL Viewer-compatible output

**Inputs:** Scan results, permutation thresholds, gene databases

**Outputs:** Peak lists, LOD plots, candidate gene lists, enrichment results, **QTL Viewer .RData object**

**Custom Scripts Used:** enrichRMod.R, peak finding and plotting functions, QTL Viewer data formatting

**Key Responsibilities:**
- Find significant peaks using empirical thresholds
- Generate Manhattan plots and LOD score visualizations
- Query gene databases for candidate genes within confidence intervals
- Perform pathway enrichment analysis
- Create publication-ready figures and summary tables
- Handle both traditional QTL and eQTL-specific outputs
- **Generate QTL Viewer-compatible .RData file following [Churchill Lab specifications](https://qtlviewer.jax.org/docs/RDataObject.pdf)**
- **Create Docker Compose configuration for local QTL Viewer deployment**
- **Provide bash script for launching local QTL Viewer web interface with study data**

## QTL Viewer Integration
The final module creates a comprehensive .RData object containing all QTL mapping results in the format required by the Churchill Lab's QTL Viewer web application. This enables:

- **Interactive visualization** of QTL mapping results through a web browser
- **Dynamic exploration** of LOD score plots, effect plots, and gene annotations
- **Local deployment** using Docker Compose for secure data analysis
- **Seamless integration** with the broader Jackson Laboratory QTL analysis ecosystem

Users can launch a local instance of the QTL Viewer populated with their study results, providing an intuitive interface for exploring and sharing QTL findings.

To Add:
- explain users need to upload data into Data/ (pheno csv, FinalReport.txt, ect)
- give credits to Excel DO_pipe for inspo -- ensure relevant Broman links are embedded here
- Note people will need to update their SLURM scheduler config in "nextflow.config"