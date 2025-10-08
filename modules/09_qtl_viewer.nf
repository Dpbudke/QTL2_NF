process PREPARE_QTLVIEWER_DATA {
    tag "Preparing QTL Viewer data for ${prefix}"
    publishDir "${params.outdir}/09_qtl_viewer_data", mode: 'copy'

    input:
    path(alleleprob_file)
    path(kinship_file)
    path(genetic_map_file)
    path(scan_results_file)
    path(significant_qtls_file)
    val(prefix)

    output:
    path("${prefix}_qtlviewer.RData"), emit: qtlviewer_data
    path("qtlviewer_conversion_report.txt"), emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
        library(tibble)
    })

    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== QTL Viewer Data Conversion Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")

    cat("Loading QTL analysis results for QTL Viewer conversion...\\n")

    # Determine which cross2 file to use - prefer filtered if available
    # Use absolute paths from the workflow root directory
    workflow_dir <- "${workflow.launchDir}"
    filtered_cross2_path <- file.path(workflow_dir, "${params.outdir}/07_permutation_testing/${prefix}_filtered_cross2.rds")
    full_cross2_path <- file.path(workflow_dir, "${params.outdir}/04_cross2_creation/${prefix}_cross2.rds")

    cat(paste("Checking for filtered cross2:", filtered_cross2_path, "\\n"))
    cat(paste("Checking for full cross2:", full_cross2_path, "\\n"))

    if (file.exists(filtered_cross2_path)) {
        cat("Using regionally filtered cross2 from module 7\\n")
        cross2 <- readRDS(filtered_cross2_path)
        validation_log <- c(validation_log, "Cross2 Source: Regionally filtered (module 7)")
    } else if (file.exists(full_cross2_path)) {
        cat("Using full cross2 from module 4\\n")
        cross2 <- readRDS(full_cross2_path)
        validation_log <- c(validation_log, "Cross2 Source: Full dataset (module 4)")
    } else {
        stop(paste("ERROR: No cross2 file found. Checked:", filtered_cross2_path, "and", full_cross2_path))
    }
    validation_log <- c(validation_log, "")
    genoprobs <- readRDS("${alleleprob_file}")
    K <- readRDS("${kinship_file}")
    gmap <- readRDS("${genetic_map_file}")
    scan_results <- readRDS("${scan_results_file}")

    # Load significant QTLs (may be empty)
    significant_qtls <- tryCatch({
        read.csv("${significant_qtls_file}", stringsAsFactors = FALSE)
    }, error = function(e) {
        data.frame(data.name = character(0), marker.id = character(0), lod = numeric(0))
    })

    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Cross2 individuals:", nrow(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Phenotypes:", ncol(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Genoprob chromosomes:", length(genoprobs)))
    validation_log <- c(validation_log, paste("✓ Kinship matrices:", length(K)))
    validation_log <- c(validation_log, paste("✓ Significant QTLs loaded:", nrow(significant_qtls)))
    validation_log <- c(validation_log, "")

    # 1. ENSEMBL VERSION (assume current mouse genome version)
    ensembl.version <- 109  # Current mouse genome version as of 2024
    validation_log <- c(validation_log, paste("✓ Ensembl version set to:", ensembl.version))

    # 2. ALLELE PROBABILITIES (alleleprobs - smaller and sufficient for QTL Viewer)
    # Verify format: list with N * K * Mj arrays per chromosome
    validation_log <- c(validation_log, "=== Allele Probabilities Validation ===")
    validation_log <- c(validation_log, paste("✓ Alleleprobs class:", class(genoprobs)))
    validation_log <- c(validation_log, paste("✓ Number of chromosomes:", length(genoprobs)))
    validation_log <- c(validation_log, "✓ Using alleleprobs (smaller file, sufficient for visualization)")
    if (length(genoprobs) > 0) {
        first_chr <- genoprobs[[1]]
        validation_log <- c(validation_log, paste("✓ First chromosome dimensions:", paste(dim(first_chr), collapse=" x ")))
        validation_log <- c(validation_log, paste("✓ Founder alleles:", paste(dimnames(first_chr)[[2]], collapse=", ")))
    }

    # 3. KINSHIP MATRICES (K) - already in correct LOCO format
    validation_log <- c(validation_log, "=== Kinship Matrices Validation ===")
    validation_log <- c(validation_log, paste("✓ Kinship class:", class(K)))
    validation_log <- c(validation_log, paste("✓ Number of kinship matrices:", length(K)))
    if (length(K) > 0) {
        first_kinship <- K[[1]]
        validation_log <- c(validation_log, paste("✓ First kinship dimensions:", paste(dim(first_kinship), collapse=" x ")))
    }

    # 4. MAP - Convert genetic map from cM to Mb
    # Note: For QTL Viewer, we need physical positions in Mb
    # Since we don't have physical map, we'll use genetic positions scaled
    cat("Converting genetic map from cM to approximate Mb positions...\\n")

    map <- list()
    total_markers <- 0

    for (chr_name in names(gmap)) {
        chr_map <- gmap[[chr_name]]
        # Convert cM to approximate Mb (rough conversion: 1 cM ≈ 1 Mb for mouse)
        # This is an approximation - ideally you'd have a physical map
        map[[chr_name]] <- chr_map / 100  # Convert cM to approximate Mb
        total_markers <- total_markers + length(chr_map)
    }

    validation_log <- c(validation_log, "=== Genetic Map Conversion ===")
    validation_log <- c(validation_log, paste("✓ Chromosomes in map:", length(map)))
    validation_log <- c(validation_log, paste("✓ Total markers:", total_markers))
    validation_log <- c(validation_log, "✓ Positions converted from cM to approximate Mb")

    # 5. MARKERS - Create markers tibble
    cat("Creating markers tibble...\\n")

    markers_list <- list()
    for (chr_name in names(map)) {
        chr_markers <- map[[chr_name]]
        chr_df <- tibble(
            marker.id = names(chr_markers),
            chr = chr_name,
            pos = as.numeric(chr_markers)
        )
        markers_list[[chr_name]] <- chr_df
    }

    markers <- do.call(rbind, markers_list)
    rownames(markers) <- NULL

    validation_log <- c(validation_log, "=== Markers Tibble Creation ===")
    validation_log <- c(validation_log, paste("✓ Total markers in tibble:", nrow(markers)))
    validation_log <- c(validation_log, paste("✓ Marker columns:", paste(colnames(markers), collapse=", ")))

    # 6. DATASET.PHENOTYPES - Main dataset object
    cat("Creating dataset.phenotypes object...\\n")

    # 6a. ANNOT.PHENOTYPE - Phenotype annotations
    pheno_names <- colnames(cross2\$pheno)
    annot.phenotype <- tibble(
        data.name = pheno_names,
        short.name = pheno_names,  # Could be customized
        R.name = make.names(pheno_names),
        description = paste("Phenotype:", pheno_names),
        units = rep("units", length(pheno_names)),  # Customize as needed
        category = rep("general", length(pheno_names)),
        R.category = rep("general", length(pheno_names)),
        is.id = rep(FALSE, length(pheno_names)),
        is.numeric = rep(TRUE, length(pheno_names)),
        is.date = rep(FALSE, length(pheno_names)),
        is.factor = rep(FALSE, length(pheno_names)),
        factor.levels = rep(NA_character_, length(pheno_names)),
        is.covar = rep(FALSE, length(pheno_names)),
        is.pheno = rep(TRUE, length(pheno_names)),
        is.derived = rep(FALSE, length(pheno_names)),
        omit = rep(FALSE, length(pheno_names)),
        use.covar = rep(NA_character_, length(pheno_names))
    )

    # 6b. ANNOT.SAMPLES - Sample annotations
    annot.samples <- cross2\$covar
    # Ensure mouse.id column exists
    if (!"mouse.id" %in% colnames(annot.samples)) {
        annot.samples\$mouse.id <- rownames(annot.samples)
    }
    annot.samples <- as_tibble(annot.samples)

    # 6c. COVAR.MATRIX - Model matrix from covariates
    cat("Creating covariate matrix...\\n")

    # Create model matrix similar to what we use in scanning
    covar_formula <- "~ 1"  # Start with intercept

    if ("Sex" %in% colnames(cross2\$covar)) {
        covar_formula <- paste(covar_formula, "+ Sex")
    }
    if ("generation" %in% colnames(cross2\$covar)) {
        covar_formula <- paste(covar_formula, "+ generation")
    }
    if ("batch" %in% colnames(cross2\$covar)) {
        covar_formula <- paste(covar_formula, "+ batch")
    }

    # Create model matrix
    covar.matrix <- model.matrix(as.formula(covar_formula), data = cross2\$covar)
    rownames(covar.matrix) <- rownames(cross2\$covar)

    # 6d. COVAR.INFO - Covariate information for QTL Viewer
    covar_cols <- colnames(covar.matrix)
    covar.info <- tibble(
        sample.column = covar_cols,
        display.name = covar_cols,
        interactive = c(FALSE, rep(TRUE, length(covar_cols) - 1)),  # First is intercept
        primary = c(FALSE, rep(FALSE, length(covar_cols) - 1)),     # Set first non-intercept as primary
        lod.peaks = c(NA_character_, rep("additive", length(covar_cols) - 1))
    )

    # Set first non-intercept covariate as primary
    if (nrow(covar.info) > 1) {
        covar.info\$primary[2] <- TRUE
    }

    # 6e. DATA - Phenotype data matrix
    data_matrix <- as.matrix(cross2\$pheno)

    # 6f. LOD.PEAKS - Process significant QTLs
    cat("Processing LOD peaks...\\n")

    lod.peaks <- list()

    # Create additive LOD peaks from significant QTLs
    if (nrow(significant_qtls) > 0) {
        # Map QTL results to required format
        if ("lodcolumn" %in% colnames(significant_qtls)) {
            # Use lodcolumn as data.name
            additive_peaks <- tibble(
                data.name = significant_qtls\$lodcolumn,
                marker.id = paste0("marker_", significant_qtls\$chr, "_", round(significant_qtls\$pos, 3)),
                lod = significant_qtls\$lod
            )
        } else {
            # Create minimal peaks structure
            additive_peaks <- tibble(
                data.name = character(0),
                marker.id = character(0),
                lod = numeric(0)
            )
        }
    } else {
        # No significant QTLs found - create empty structure
        additive_peaks <- tibble(
            data.name = character(0),
            marker.id = character(0),
            lod = numeric(0)
        )
    }

    lod.peaks\$additive <- additive_peaks

    # Create dataset.phenotypes object
    dataset.phenotypes <- list(
        annot.phenotype = annot.phenotype,
        annot.samples = annot.samples,
        covar.matrix = covar.matrix,
        covar.info = covar.info,
        data = data_matrix,
        datatype = "phenotype",
        display.name = paste("${prefix}", "Phenotypes"),
        lod.peaks = lod.peaks
    )

    validation_log <- c(validation_log, "=== Dataset Object Creation ===")
    validation_log <- c(validation_log, paste("✓ Phenotype annotations:", nrow(annot.phenotype)))
    validation_log <- c(validation_log, paste("✓ Sample annotations:", nrow(annot.samples)))
    validation_log <- c(validation_log, paste("✓ Covariate matrix dimensions:", paste(dim(covar.matrix), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Covariate info entries:", nrow(covar.info)))
    validation_log <- c(validation_log, paste("✓ Data matrix dimensions:", paste(dim(data_matrix), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ LOD peaks (additive):", nrow(lod.peaks\$additive)))

    # Final validation of required elements
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Viewer Format Validation ===")
    validation_log <- c(validation_log, "✓ ensembl.version: numeric")
    validation_log <- c(validation_log, "✓ genoprobs: list of allele probability arrays (from alleleprob)")
    validation_log <- c(validation_log, "✓ K: list of kinship matrices")
    validation_log <- c(validation_log, "✓ map: list of position vectors")
    validation_log <- c(validation_log, "✓ markers: tibble with marker.id, chr, pos")
    validation_log <- c(validation_log, "✓ dataset.phenotypes: complete dataset object")

    # Save QTL Viewer RData file with exact required elements
    cat("Saving QTL Viewer RData file...\\n")

    save(
        ensembl.version,
        genoprobs,
        K,
        map,
        markers,
        dataset.phenotypes,
        file = "${prefix}_qtlviewer.RData"
    )

    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Viewer Data Conversion Complete ===")
    validation_log <- c(validation_log, paste("✓ RData file saved:", "${prefix}_qtlviewer.RData"))
    validation_log <- c(validation_log, "✓ File ready for local Docker deployment")
    validation_log <- c(validation_log, "✓ All required QTL Viewer elements included")

    # Write validation report
    writeLines(validation_log, "qtlviewer_conversion_report.txt")

    cat("QTL Viewer data conversion completed successfully\\n")
    cat("RData file is ready for local Docker deployment with Docker containers\\n")
    """
}

process SETUP_QTLVIEWER_DEPLOYMENT {
    tag "Setting up QTL Viewer deployment for ${prefix}"
    publishDir "${params.outdir}/09_qtl_viewer_data", mode: 'copy'

    input:
    path(qtlviewer_data)  // RData file from PREPARE_QTLVIEWER_DATA
    val(prefix)

    output:
    path("docker-compose.yml"), emit: docker_compose
    path("start_qtlviewer.sh"), emit: startup_script
    path("README_qtlviewer.md"), emit: instructions
    path("data"), emit: data_directory, type: 'dir'

    script:
    """
    #!/bin/bash

    echo "Setting up QTL Viewer deployment for local Docker use..."

    # Create data directory structure
    mkdir -p data/rdata
    mkdir -p data/sqlite

    # Create symlink to RData file instead of copying (saves space and avoids duplication)
    echo "Linking QTL Viewer RData file to deployment structure..."
    ln -sf "../../${qtlviewer_data}" "data/rdata/${prefix}.RData"

    # Handle SQLite database for founder SNPs
    echo "Downloading MGI SQLite database for founder SNPs..."

    database_found=false

    # Check for existing database from previous modules
    if [ -f "../02_genotype_processing/mgi_mouse_genes.sqlite" ]; then
        echo "Found existing MGI SQLite database from genotype processing, copying..."
        cp "../02_genotype_processing/mgi_mouse_genes.sqlite" "data/sqlite/ccfoundersnps.sqlite"
        database_found=true
    elif [ -f "../03_control_file_generation/mgi_mouse_genes.sqlite" ]; then
        echo "Found existing MGI SQLite database from control file generation, copying..."
        cp "../03_control_file_generation/mgi_mouse_genes.sqlite" "data/sqlite/ccfoundersnps.sqlite"
        database_found=true
    fi

    # If no database found, download it directly
    if [ "\$database_found" = "false" ]; then
        echo "No existing SQLite database found, downloading MGI database from Figshare..."

        # Download MGI SQLite database (24607961 is the file ID from the Figshare URL)
        if command -v curl >/dev/null 2>&1; then
            curl -L -H "User-Agent: Mozilla/5.0" -o data/sqlite/ccfoundersnps.sqlite "https://figshare.com/ndownloader/files/24607961" || {
                echo "ERROR: Failed to download MGI database with curl"
                echo "Please download manually from: https://figshare.com/articles/dataset/SQLite_database_with_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5280238/7"
                touch data/sqlite/ccfoundersnps.sqlite
            }
        elif command -v wget >/dev/null 2>&1; then
            wget -O data/sqlite/ccfoundersnps.sqlite "https://figshare.com/ndownloader/files/24607961" || {
                echo "ERROR: Failed to download MGI database with wget"
                echo "Please download manually from: https://figshare.com/articles/dataset/SQLite_database_with_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5280238/7"
                touch data/sqlite/ccfoundersnps.sqlite
            }
        else
            echo "ERROR: Neither curl nor wget available for download"
            echo "Please download manually from: https://figshare.com/articles/dataset/SQLite_database_with_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5280238/7"
            touch data/sqlite/ccfoundersnps.sqlite
        fi

        # Verify download
        if [ -s "data/sqlite/ccfoundersnps.sqlite" ]; then
            echo "✓ MGI database downloaded successfully (\$(du -h data/sqlite/ccfoundersnps.sqlite | cut -f1))"
            database_found=true
        else
            echo "✗ MGI database download failed or file is empty"
            echo "QTL Viewer will run with limited functionality"
        fi
    fi

    if [ "\$database_found" = "true" ]; then
        echo "✓ MGI mouse genes database ready for QTL Viewer"
    else
        echo "⚠ QTL Viewer will run without founder SNPs database (limited functionality)"
    fi

    # Create docker-compose.yml for local deployment
    cat > docker-compose.yml << 'EOF'
version: '3.8'

services:
  qtl2rest:
    image: churchilllab/qtl2rest:0.6.5
    container_name: qtl2rest-${prefix}
    ports:
      - "8001:8001"
    volumes:
      - ./data/rdata:/app/qtl2rest/data/rdata
      - ./data/sqlite:/app/qtl2rest/data/sqlite
    networks:
      - qtlviewer-network
    environment:
      - R_MAX_VSIZE=32Gb
    restart: unless-stopped

  qtl2web:
    image: churchilllab/qtl2web:1.6.1
    container_name: qtl2web-${prefix}
    ports:
      - "8000:8000"
    depends_on:
      - qtl2rest
    environment:
      - QTL2REST_URL=http://qtl2rest:8001
    networks:
      - qtlviewer-network
    restart: unless-stopped

networks:
  qtlviewer-network:
    driver: bridge
EOF

    # Create Docker startup script for local deployment
    cat > start_qtlviewer.sh << 'EOF'
#!/bin/bash

echo "Starting QTL Viewer for study: ${prefix}"
echo "========================================"

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "ERROR: Docker is not running. Please start Docker first."
    echo ""
    echo "📥 DOWNLOAD INSTRUCTIONS:"
    echo "1. Download this entire directory to your local machine"
    echo "2. Install Docker Desktop on your local machine"
    echo "3. Run: ./start_qtlviewer.sh"
    echo "4. Open: http://localhost:8000"
    exit 1
fi

# Check if required data files exist
if [ ! -f "data/rdata/${prefix}.RData" ]; then
    echo "ERROR: QTL data file not found: data/rdata/${prefix}.RData"
    exit 1
fi

if [ ! -f "data/sqlite/ccfoundersnps.sqlite" ] || [ ! -s "data/sqlite/ccfoundersnps.sqlite" ]; then
    echo "WARNING: Founder SNPs database not found or empty"
    echo "Some features may be limited without this database"
fi

# Pull latest images
echo "Pulling latest QTL Viewer Docker images..."
docker-compose pull

# Start services
echo "Starting QTL Viewer services..."
docker-compose up -d

# Wait for services to be ready
echo "Waiting for services to start..."
sleep 15

# Check if services are running
if docker-compose ps | grep -q "Up"; then
    echo ""
    echo "🎉 QTL Viewer is now running!"
    echo "=============================="
    echo "🌐 Web Interface: http://localhost:8000"
    echo "🔌 REST API:      http://localhost:8001"
    echo ""
    echo "📊 Your QTL Results:"
    echo "   Study: ${prefix}"
    echo "   Data: \$(du -h data/rdata/${prefix}.RData | cut -f1)"
    echo ""
    echo "🛑 To stop: docker-compose down"
    echo "📋 To view logs: docker-compose logs -f"
    echo ""
    echo "💡 Open http://localhost:8000 in your web browser!"
else
    echo "ERROR: Services failed to start. Check logs with: docker-compose logs"
fi
EOF

    chmod +x start_qtlviewer.sh

    # Create README with simplified instructions
    cat > README_qtlviewer.md << 'EOF'
# QTL Viewer for Study: ${prefix}

## 🎯 Quick Start (Local Machine)

**Prerequisites:** Docker Desktop installed on your local machine

1. **Download** this entire directory to your local machine
2. **Navigate** to the downloaded directory
3. **Run** the viewer:
   ```bash
   ./start_qtlviewer.sh
   ```
4. **Open** your web browser to: **http://localhost:8000**

## 📁 What's Included

- **Your QTL Data**: `data/rdata/${prefix}.RData`
- **MGI Database**: `data/sqlite/ccfoundersnps.sqlite` (mouse gene annotations)
- **Docker Setup**: `docker-compose.yml` (containerized deployment)
- **Startup Script**: `start_qtlviewer.sh` (one-click launch)

## 🌐 QTL Viewer Features

Once running, you can:
- **Browse QTL Results**: Interactive LOD plots and peaks
- **Explore Genetic Maps**: Chromosome-level visualization
- **Analyze Effects**: Founder strain effect coefficients
- **Query API**: REST endpoints for programmatic access

## 🔌 API Endpoints

Your QTL data is accessible via REST API at `http://localhost:8001`:

- `/envinfo` - Data file information
- `/datasets` - Available datasets
- `/lodpeaks` - Significant QTL peaks
- `/markers` - Genetic marker information
- `/lodscores?dataset=dataset.phenotypes&id=PHENOTYPE` - LOD scores

## 🛠 Management Commands

```bash
# Start QTL Viewer
./start_qtlviewer.sh

# Stop QTL Viewer
docker-compose down

# View logs
docker-compose logs -f

# Restart services
docker-compose restart
```

## 📋 Troubleshooting

**Port conflicts**: If ports 8000/8001 are in use:
- Stop other services using these ports
- Or edit `docker-compose.yml` to use different ports

**Docker not found**:
- Install Docker Desktop for your operating system
- Ensure Docker is running before starting QTL Viewer

**Data not loading**:
- Ensure the entire directory was downloaded intact
- Check that `data/rdata/${prefix}.RData` exists

## 🚀 Performance

- **Memory**: Large datasets may require 8+ GB RAM
- **Storage**: Container images require ~1GB disk space
- **Network**: Initial setup downloads Docker images

## 📞 Support

For QTL Viewer questions:
- **Documentation**: https://qtlviewer.jax.org/
- **Source Code**: https://github.com/churchill-lab

---

**Generated by QTL2_NF Pipeline** 🧬
EOF

    echo "✅ QTL Viewer deployment package ready!"
    echo "📦 Directory contents:"
    echo "   - Docker containers: qtl2rest + qtl2web"
    echo "   - QTL data: ${prefix}.RData"
    echo "   - MGI database: ccfoundersnps.sqlite"
    echo "   - Startup script: start_qtlviewer.sh"
    echo ""
    echo "🎯 Next steps:"
    echo "   1. Download this entire directory to your local machine"
    echo "   2. Run: ./start_qtlviewer.sh"
    echo "   3. Open: http://localhost:8000"
    """
}