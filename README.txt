# QTL2_NF Pipeline Overview

## Module 1: Data Validation & Phenotype Processing (phenotype_process.nf)
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

## Module 2: Genotype File Processing (genotype_process.nf)
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

## Module 3: Control File Generation & Cross2 Creation (control_cross2.nf)
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

## Module 4: Genome Scanning & Permutations (scan_perm.nf)
**Purpose:** Perform comprehensive QTL/eQTL mapping with statistical significance testing using r/qtl2

**Inputs:** cross2 object from Module 3

**Outputs:** Genotype probabilities, kinship matrices, genome scan results, permutation thresholds, significant QTL lists

**Four-Process Architecture:**

### Process 1: PREPARE_GENOME_SCAN (16 CPUs, 256GB, 12h)
- Insert pseudomarkers at 0.5 cM resolution for high-resolution mapping
- Calculate genotype probabilities using calc_genoprob() with error_prob=0.002
- Generate LOCO (Leave One Chromosome Out) kinship matrices using calc_kinship()
- Save intermediate files: genoprob.rds, kinship_loco.rds, genetic_map.rds

### Process 2: GENOME_SCAN (32 CPUs, 512GB, 24h)  
- Perform genome scan using scan1() with linear mixed models
- Analyze ALL phenotypes simultaneously (phenotype-agnostic pipeline)
- Apply LOCO kinship correction for population structure
- Handle covariates: Sex (required for X chr), Generation, Batch (if available)
- **NEW: Smart LOD threshold filtering** - identify phenotypes with peak LOD ≥ threshold (default: 7.0)
- Find preliminary peaks (LOD ≥ 3.0) using find_peaks()
- Outputs: scan_results.rds, preliminary peaks CSV, **filtered_phenotypes.txt**

### Process 3: PERMUTATION_TEST (48 CPUs, 1TB, 72h) - **OPTIMIZED FOR EFFICIENCY**
- **NEW: Only runs permutations on phenotypes passing LOD threshold filter**
- Execute 1000 permutation rounds using scan1perm() on **filtered subset only**
- **Massive computational savings**: Typically 70-90% reduction in permutation workload
- Use all available CPU cores (cores=0) for maximum parallelization
- Proper permutation strategy: genotype permutation with kinship matrices
- Calculate significance thresholds for 4 levels: 0.63 (suggestive), 0.10, 0.05, 0.01
- Outputs: permutation_results.rds, significance_thresholds.csv

### Process 4: IDENTIFY_SIGNIFICANT_QTLS (8 CPUs, 64GB, 4h)
- Apply empirical significance thresholds to identify significant QTLs
- Generate comprehensive QTL lists for all significance levels
- Create summary statistics and chromosome distribution reports
- Filter and report top QTLs (LOD ≥ 10) for immediate attention
- Outputs: significant_qtls.csv, qtl_summary.txt

**Output Directory Structure:**
- `05_genome_scan_preparation/`: Genotype probabilities and kinship matrices
- `06_genome_scan_results/`: Genome scan results and preliminary peaks
- `07_permutation_testing/`: Permutation results and significance thresholds
- `08_significant_qtls/`: Final significant QTL identifications and summaries
- `09_qtl_viewer_data/`: **QTL Viewer deployment package (fully portable)**

**Performance Optimizations:**
- High-performance computing resource allocation (up to 48 CPUs, 1TB RAM)
- Parallel processing using all available cores
- Memory-efficient handling of large-scale expression datasets
- LOCO kinship approach for proper mixed model analysis
- **NEW: Smart permutation filtering** - user-configurable LOD threshold (--lod_threshold) for computational efficiency

**User Parameters:**
- `--lod_threshold`: Minimum LOD score for permutation testing (default: 7.0)
  - Higher values = fewer phenotypes tested = faster runtime
  - Lower values = more comprehensive testing = longer runtime
  - Recommended range: 5.0-10.0 depending on study goals and computational resources

## Module 5: QTL Viewer Integration & Deployment (qtl_viewer.nf)
**Purpose:** Convert QTL analysis results to interactive QTL Viewer format and create portable deployment package

**Inputs:** cross2 object, genotype probabilities, kinship matrices, genetic map, scan results, significant QTLs

**Outputs:** QTL Viewer .RData object, Docker Compose deployment package, startup scripts

**Two-Process Architecture:**

### Process 1: PREPARE_QTLVIEWER_DATA (8 CPUs, 64GB, 2h)
- **Convert QTL results to official QTL Viewer RData format** following [Churchill Lab specifications](https://qtlviewer.jax.org/docs/RDataObject.pdf)
- Create all required elements: `ensembl.version`, `genoprobs`, `K`, `map`, `markers`, `dataset.phenotypes`
- Build proper dataset structure with phenotype annotations, sample metadata, covariate matrices
- Convert genetic map positions from cM to approximate Mb coordinates
- Format LOD peaks data for interactive visualization
- Generate complete validation report for format compliance

### Process 2: SETUP_QTLVIEWER_DEPLOYMENT (1 CPU, 8GB, 30min)
- **Create portable Docker deployment package** that can be downloaded and run anywhere
- Generate `docker-compose.yml` with qtl2rest API backend + qtl2web frontend services
- **Reuse existing MGI SQLite database** from Module 2/3 (no redundant downloads)
- Create `start_qtlviewer.sh` one-command startup script
- Generate comprehensive documentation and troubleshooting guide
- Structure data directory for qtl2rest container volume mounting

**Key Responsibilities:**
- **Exact QTL Viewer format compliance** - all required data elements with proper structure
- **Portable deployment** - entire package can be downloaded from HPC to local machine
- **Database reuse** - leverages existing MGI mouse gene database from earlier modules
- **One-command startup** - simple `./start_qtlviewer.sh` launches complete web application
- **Docker orchestration** - qtl2rest (port 8001) + qtl2web (port 8000) services
- **Documentation** - complete setup and usage instructions included

## QTL Viewer Integration & Deployment
Module 5 creates a **fully portable QTL Viewer deployment package** that enables:

### **Interactive Web Application Features:**
- **LOD Score Visualization**: Interactive genome-wide association plots
- **Effect Plots**: Founder strain coefficient visualization
- **SNP Association**: High-resolution association mapping
- **Correlation Analysis**: Cross-phenotype correlation exploration
- **Mediation Analysis**: Causal relationship investigation

### **Deployment Options:**
- **HPC Deployment**: Run directly on cluster where pipeline completed
- **Local Download**: Copy entire `09_qtl_viewer_data/` directory to any machine with Docker
- **Permanent Access**: QTL Viewer instance independent of HPC session status

### **Simple Usage:**
```bash
# Option 1: On HPC
cd results/09_qtl_viewer_data/
./start_qtlviewer.sh

# Option 2: Download to local machine
scp -r user@hpc:path/to/results/09_qtl_viewer_data/ ~/my_qtl_study/
cd ~/my_qtl_study/09_qtl_viewer_data/
./start_qtlviewer.sh
```

**Result**: Interactive QTL Viewer accessible at `http://localhost:8000` with complete study data

To Add:
- explain users need to upload data into Data/ (pheno csv, FinalReport.txt, ect)
- give credits to Excel DO_pipe for inspo -- ensure relevant Broman links are embedded here
- Note people will need to update their SLURM scheduler config in "nextflow.config"
- Note I'm not including the downloading of reference files in the prepatory step since we want to seamlessly integrate new refs when they are published (ie just uping the download link in Module2)
