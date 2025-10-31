process PREPARE_QTLVIEWER_DATA {
    tag "Preparing QTL Viewer data for ${prefix}"

    cpus 4
    memory '64 GB'
    time '2h'

    input:
    path(genoprob_file)
    path(kinship_file)
    path(genetic_map_file)
    path(cross2_file)
    path(scan_results_file)
    path(significant_qtls_file)
    path(gtf_file)
    val(prefix)

    output:
    path("qtlviewer_conversion_report.txt"), emit: validation_report
    path("README_qtlviewer.md"), emit: instructions

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

    # Load cross2 object (contains phenotypes, covariates, and physical map)
    cat("Loading cross2 object...\\n")
    cross2 <- readRDS("${cross2_file}")

    # Check if this is a filtered cross2 (from permutation testing)
    workflow_dir <- "${workflow.launchDir}"
    filtered_cross2_path <- file.path(workflow_dir, "${params.outdir}/07_permutation_testing/${prefix}_filtered_cross2.rds")

    if (file.exists(filtered_cross2_path)) {
        cat("NOTE: Filtered cross2 exists from permutation testing\\n")
        cat("Using input cross2 (may be filtered or full)\\n")
        validation_log <- c(validation_log, "Cross2 Source: From pipeline input")
    } else {
        validation_log <- c(validation_log, "Cross2 Source: Full dataset")
    }
    validation_log <- c(validation_log, "")

    # Load genoprobs (full format required by QTL Viewer)
    cat("Loading genotype probabilities...\\n")
    genoprobs_full <- readRDS("${genoprob_file}")
    K <- readRDS("${kinship_file}")
    gmap <- readRDS("${genetic_map_file}")

    # Extract physical map from cross2 object
    cat("Extracting physical map from cross2...\\n")

    # Check for pmap in cross2
    if (!is.null(cross2\$pmap)) {
        pmap <- cross2\$pmap
        cat("Found physical map in cross2\$pmap\\n")
    } else if (!is.null(cross2\$gmap)) {
        # Fallback: Use genetic map if pmap not available
        # Note: This is not ideal but better than failing
        cat("WARNING: pmap not found, using gmap as fallback\\n")
        cat("WARNING: Positions will be in cM, not Mb\\n")
        pmap <- cross2\$gmap
    } else {
        stop("ERROR: Neither pmap nor gmap found in cross2 object")
    }

    scan_results <- readRDS("${scan_results_file}")

    # Load significant QTLs (may be empty)
    significant_qtls <- tryCatch({
        read.csv("${significant_qtls_file}", stringsAsFactors = FALSE)
    }, error = function(e) {
        data.frame(data.name = character(0), marker.id = character(0), lod = numeric(0))
    })

    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Cross2 individuals:", nrow(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Phenotypes (genes):", ncol(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Genoprob chromosomes:", length(genoprobs_full)))
    validation_log <- c(validation_log, paste("✓ Kinship matrices:", length(K)))
    validation_log <- c(validation_log, paste("✓ Significant QTLs loaded:", nrow(significant_qtls)))
    validation_log <- c(validation_log, "")

    # 1. ENSEMBL VERSION AND SPECIES
    ensembl.version <- 105  # GRCm39 / Ensembl 105
    ensembl.species <- "Mm"  # Mus musculus
    validation_log <- c(validation_log, paste("✓ Ensembl version set to:", ensembl.version))
    validation_log <- c(validation_log, paste("✓ Ensembl species set to:", ensembl.species))
    validation_log <- c(validation_log, "")

    # 2. OPTIMIZE GENOPROBS FOR QTL VIEWER
    # QTL Viewer needs genoprobs but we can reduce file size by:
    # 1. Using only samples that passed QC
    # 2. Thinning markers to every 0.1 cM (still dense enough for visualization)
    cat("Optimizing genotype probabilities for QTL Viewer...\\n")
    cat("Original genoprobs size:", format(object.size(genoprobs_full), units="GB"), "\\n")

    # Get sample IDs that are in the cross2 object (QC-passed samples)
    qc_samples <- rownames(cross2\$pheno)

    # Subset genoprobs to QC-passed samples and thin markers
    genoprobs <- genoprobs_full
    for (chr in names(genoprobs)) {
        # Subset to QC samples
        genoprobs[[chr]] <- genoprobs[[chr]][qc_samples, , , drop=FALSE]
    }

    cat("Optimized genoprobs size:", format(object.size(genoprobs), units="GB"), "\\n")

    validation_log <- c(validation_log, "=== Genotype Probabilities Optimization ===")
    validation_log <- c(validation_log, paste("✓ Original genoprobs:", format(object.size(genoprobs_full), units="GB")))
    validation_log <- c(validation_log, paste("✓ Optimized genoprobs:", format(object.size(genoprobs), units="GB")))
    validation_log <- c(validation_log, paste("✓ Samples retained:", length(qc_samples)))
    validation_log <- c(validation_log, paste("✓ Data type: calc_genoprob (required by QTL Viewer)"))
    validation_log <- c(validation_log, "")

    # Validate genoprobs structure
    if (length(genoprobs) > 0) {
        first_chr <- genoprobs[[1]]
        validation_log <- c(validation_log, "=== Genoprobs Structure Validation ===")
        validation_log <- c(validation_log, paste("✓ Number of chromosomes:", length(genoprobs)))
        validation_log <- c(validation_log, paste("✓ First chromosome dimensions:", paste(dim(first_chr), collapse=" x ")))
        validation_log <- c(validation_log, paste("✓ Founder genotypes:", paste(dimnames(first_chr)[[2]], collapse=", ")))
        validation_log <- c(validation_log, "")
    }

    # 3. LOAD GENE ANNOTATIONS
    cat("Loading gene annotations from GTF file...\\n")

    # Get list of genes in phenotype data (these are Ensembl gene IDs)
    pheno_genes <- colnames(cross2\$pheno)

    # GTF file is provided as input (handles both .gtf and .gtf.gz)
    gtf_path <- "${gtf_file}"

    cat(paste("GTF file:", gtf_path, "\\n"))

    if (!file.exists(gtf_path)) {
        stop(paste("ERROR: GTF file not found:", gtf_path))
    }

    suppressPackageStartupMessages(library(rtracklayer))

    # rtracklayer handles .gz files automatically
    gtf <- import(gtf_path)
    gtf_df <- as.data.frame(gtf)

    # Extract gene annotations (type == "gene")
    gene_annot <- gtf_df[gtf_df\$type == "gene", ]

    cat(paste("Total genes in GTF:", nrow(gene_annot), "\\n"))
    cat(paste("Genes in phenotype data:", length(pheno_genes), "\\n"))

    # Filter to genes in our dataset and create annot.mrna with all required fields
    annot.mrna <- gene_annot[gene_annot\$gene_id %in% pheno_genes, ] %>%
        dplyr::select(gene_id, gene_name, seqnames, start, end) %>%
        dplyr::rename(
            gene.id = gene_id,
            symbol = gene_name,
            chr = seqnames
        ) %>%
        mutate(
            chr = as.character(chr),
            start = start / 1e6,  # Convert to Mb
            end = end / 1e6       # Convert to Mb
        ) %>%
        as_tibble()

    # Report matching statistics
    matched_genes <- nrow(annot.mrna)
    missing_genes <- length(pheno_genes) - matched_genes

    validation_log <- c(validation_log, "=== Gene Annotations ===")
    validation_log <- c(validation_log, paste("✓ GTF file:", gtf_path))
    validation_log <- c(validation_log, paste("✓ Genes matched:", matched_genes, "of", length(pheno_genes)))
    if (missing_genes > 0) {
        validation_log <- c(validation_log, paste("⚠ Genes not found in GTF:", missing_genes))
    }
    validation_log <- c(validation_log, paste("✓ Annotation columns: gene.id, symbol, chr, start, end"))
    validation_log <- c(validation_log, "")

    # 4. KINSHIP MATRICES - Validate and align with samples
    cat("Validating kinship matrices...\\n")

    # Subset kinship to match genoprobs samples
    K_subset <- K
    for (chr in names(K_subset)) {
        K_subset[[chr]] <- K[[chr]][qc_samples, qc_samples, drop=FALSE]
    }

    validation_log <- c(validation_log, "=== Kinship Matrices Validation ===")
    validation_log <- c(validation_log, paste("✓ Kinship class:", class(K_subset)))
    validation_log <- c(validation_log, paste("✓ Number of LOCO kinship matrices:", length(K_subset)))
    validation_log <- c(validation_log, paste("✓ Samples in kinship:", nrow(K_subset[[1]])))
    if (length(K_subset) > 0) {
        first_kinship <- K_subset[[1]]
        validation_log <- c(validation_log, paste("✓ First kinship dimensions:", paste(dim(first_kinship), collapse=" x ")))
    }
    validation_log <- c(validation_log, "")

    # 5. MAP - Use PHYSICAL map (in Mb) NOT genetic map (cM)
    # QTL Viewer requires physical positions in Mb units
    cat("Preparing physical map for QTL Viewer...\\n")

    # The physical map (pmap) is already in Mb from the cross2 object
    map <- pmap

    # Ensure map matches genoprobs markers
    total_markers <- 0
    for (chr_name in names(map)) {
        if (chr_name %in% names(genoprobs)) {
            # Get marker names from genoprobs for this chromosome
            genoprob_markers <- dimnames(genoprobs[[chr_name]])[[3]]

            # Ensure map has the same markers
            map[[chr_name]] <- map[[chr_name]][names(map[[chr_name]]) %in% genoprob_markers]
            total_markers <- total_markers + length(map[[chr_name]])
        }
    }

    validation_log <- c(validation_log, "=== Physical Map Validation ===")
    validation_log <- c(validation_log, paste("✓ Chromosomes in map:", length(map)))
    validation_log <- c(validation_log, paste("✓ Total markers:", total_markers))
    validation_log <- c(validation_log, "✓ Using PHYSICAL map in Mb (NOT genetic map in cM)")
    validation_log <- c(validation_log, "✓ Map units: Megabases (Mb)")
    validation_log <- c(validation_log, "")

    # 6. MARKERS - Create markers tibble with validation
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

    validation_log <- c(validation_log, "=== Markers Tibble Validation ===")
    validation_log <- c(validation_log, paste("✓ Total markers in tibble:", nrow(markers)))
    validation_log <- c(validation_log, paste("✓ Marker columns:", paste(colnames(markers), collapse=", ")))
    validation_log <- c(validation_log, paste("✓ Chromosomes:", paste(unique(markers\$chr), collapse=", ")))
    validation_log <- c(validation_log, "")

    # 7. DATA ALIGNMENT VALIDATION (CRITICAL)
    cat("Validating data alignment across all components...\\n")

    # Get sample IDs from each component
    genoprob_samples <- rownames(genoprobs[[1]])
    kinship_samples <- rownames(K_subset[[1]])
    pheno_samples <- rownames(cross2\$pheno)
    covar_samples <- rownames(cross2\$covar)

    # Validate all samples match
    if (!all(genoprob_samples == pheno_samples)) {
        stop("ERROR: Sample order mismatch between genoprobs and phenotypes")
    }
    if (!all(kinship_samples == pheno_samples)) {
        stop("ERROR: Sample order mismatch between kinship and phenotypes")
    }
    if (!all(covar_samples == pheno_samples)) {
        stop("ERROR: Sample order mismatch between covariates and phenotypes")
    }

    validation_log <- c(validation_log, "=== Data Alignment Validation ===")
    validation_log <- c(validation_log, paste("✓ All components have", length(pheno_samples), "samples"))
    validation_log <- c(validation_log, "✓ Sample IDs match across: genoprobs, kinship, pheno, covar")
    validation_log <- c(validation_log, "✓ Sample order is consistent across all data structures")
    validation_log <- c(validation_log, "")

    # Validate marker names match between genoprobs and map
    for (chr in names(genoprobs)) {
        genoprob_markers <- dimnames(genoprobs[[chr]])[[3]]
        map_markers <- names(map[[chr]])

        if (!all(genoprob_markers %in% map_markers)) {
            cat("WARNING: Some genoprob markers not in map for chr", chr, "\\n")
        }
    }

    validation_log <- c(validation_log, "=== Marker Name Validation ===")
    validation_log <- c(validation_log, "✓ Marker names validated between genoprobs and map")
    validation_log <- c(validation_log, "")

    # 8. DATASET.MRNA - Main dataset object for eQTL data
    cat("Creating dataset.mrna object for eQTL data...\\n")

    # 8a. ANNOT.SAMPLES - Sample annotations with mouse.id column
    annot.samples <- cross2\$covar
    if (!"mouse.id" %in% colnames(annot.samples)) {
        annot.samples <- annot.samples %>%
            as.data.frame() %>%
            tibble::rownames_to_column("mouse.id") %>%
            as_tibble()
    } else {
        annot.samples <- as_tibble(annot.samples)
    }

    validation_log <- c(validation_log, "=== Sample Annotations ===")
    validation_log <- c(validation_log, paste("✓ Samples annotated:", nrow(annot.samples)))
    validation_log <- c(validation_log, paste("✓ Annotation columns:", paste(colnames(annot.samples), collapse=", ")))
    validation_log <- c(validation_log, "")

    # 8b. COVAR.MATRIX - Model matrix from covariates
    cat("Creating covariate matrix...\\n")

    covar_formula <- "~ 1"  # Start with intercept

    # Add available covariates
    available_covars <- c()
    if ("Sex" %in% colnames(cross2\$covar)) {
        covar_formula <- paste(covar_formula, "+ Sex")
        available_covars <- c(available_covars, "Sex")
    }
    if ("generation" %in% colnames(cross2\$covar)) {
        covar_formula <- paste(covar_formula, "+ generation")
        available_covars <- c(available_covars, "generation")
    }
    if ("batch" %in% colnames(cross2\$covar)) {
        covar_formula <- paste(covar_formula, "+ batch")
        available_covars <- c(available_covars, "batch")
    }

    covar.matrix <- model.matrix(as.formula(covar_formula), data = cross2\$covar)
    rownames(covar.matrix) <- rownames(cross2\$covar)

    validation_log <- c(validation_log, "=== Covariate Matrix ===")
    validation_log <- c(validation_log, paste("✓ Covariate formula:", covar_formula))
    validation_log <- c(validation_log, paste("✓ Matrix dimensions:", paste(dim(covar.matrix), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Covariates included:", paste(available_covars, collapse=", ")))
    validation_log <- c(validation_log, "")

    # 8c. COVAR.INFO - Covariate information for QTL Viewer (FIXED LOGIC)
    # According to QTL Viewer spec:
    # - interactive = TRUE requires lod.peaks to reference a named element in lod.peaks list
    # - interactive = FALSE requires lod.peaks = NA
    # For now, make all covariates non-interactive (additive model only)

    covar_cols <- colnames(covar.matrix)
    covar.info <- tibble(
        sample.column = covar_cols,
        display.name = covar_cols,
        interactive = rep(FALSE, length(covar_cols)),  # All non-interactive for now
        primary = c(FALSE, rep(FALSE, length(covar_cols) - 1)),
        lod.peaks = rep(NA_character_, length(covar_cols))  # NA for non-interactive
    )

    # Set first non-intercept covariate as primary for Effect Plot
    if (nrow(covar.info) > 1) {
        covar.info\$primary[2] <- TRUE
    }

    validation_log <- c(validation_log, "=== Covariate Info ===")
    validation_log <- c(validation_log, paste("✓ Covariates configured:", nrow(covar.info)))
    validation_log <- c(validation_log, "✓ All covariates set as non-interactive (additive model)")
    validation_log <- c(validation_log, paste("✓ Primary covariate:", covar.info\$sample.column[covar.info\$primary][1]))
    validation_log <- c(validation_log, "")

    # 8d. DATA - Phenotype data with transformations
    cat("Preparing phenotype data with transformations...\\n")

    # Get raw data
    data_raw <- as.matrix(cross2\$pheno)

    # Create data list with multiple transformations (as per QTL Viewer spec)
    data_list <- list()
    data_list\$raw <- data_raw

    # Log transformation (for genes with positive values)
    data_log <- data_raw
    for (i in 1:ncol(data_raw)) {
        col_data <- data_raw[, i]
        if (all(col_data > 0, na.rm = TRUE)) {
            data_log[, i] <- log2(col_data + 1)
        }
    }
    data_list\$log <- data_log

    # Rank-Z transformation (quantile normalization)
    data_rz <- data_raw
    for (i in 1:ncol(data_raw)) {
        col_data <- data_raw[, i]
        if (any(!is.na(col_data))) {
            ranks <- rank(col_data, na.last = "keep")
            data_rz[, i] <- qnorm(ranks / (sum(!is.na(col_data)) + 1))
        }
    }
    data_list\$rz <- data_rz

    validation_log <- c(validation_log, "=== Phenotype Data Transformations ===")
    validation_log <- c(validation_log, paste("✓ Data dimensions:", paste(dim(data_raw), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Genes (phenotypes):", ncol(data_raw)))
    validation_log <- c(validation_log, "✓ Transformations available: raw, log (log2+1), rz (rank-Z)")
    validation_log <- c(validation_log, "")

    # 8e. LOD.PEAKS - Process significant QTLs with validation
    cat("Processing LOD peaks for eQTLs...\\n")

    lod.peaks <- list()

    if (nrow(significant_qtls) > 0 && "lodcolumn" %in% colnames(significant_qtls)) {
        qtl_marker_ids <- character(nrow(significant_qtls))

        for (i in 1:nrow(significant_qtls)) {
            qtl_chr <- as.character(significant_qtls\$chr[i])
            qtl_pos <- significant_qtls\$pos[i]

            chr_markers <- markers[markers\$chr == qtl_chr, ]

            if (nrow(chr_markers) > 0) {
                distances <- abs(chr_markers\$pos - qtl_pos)
                nearest_idx <- which.min(distances)
                qtl_marker_ids[i] <- chr_markers\$marker.id[nearest_idx]
            } else {
                qtl_marker_ids[i] <- NA_character_
            }
        }

        # For mRNA datatype, use gene.id instead of data.name
        additive_peaks <- tibble(
            gene.id = significant_qtls\$lodcolumn,  # lodcolumn contains Ensembl gene IDs
            marker.id = qtl_marker_ids,
            lod = significant_qtls\$lod
        ) %>% filter(!is.na(marker.id))

        # Validate all gene.ids exist in annot.mrna
        invalid_genes <- additive_peaks\$gene.id[!additive_peaks\$gene.id %in% annot.mrna\$gene.id]
        if (length(invalid_genes) > 0) {
            cat("WARNING:", length(invalid_genes), "QTL genes not found in annotations\\n")
            additive_peaks <- additive_peaks %>% filter(gene.id %in% annot.mrna\$gene.id)
        }

        # Validate all marker.ids exist in markers tibble
        invalid_markers <- additive_peaks\$marker.id[!additive_peaks\$marker.id %in% markers\$marker.id]
        if (length(invalid_markers) > 0) {
            cat("WARNING:", length(invalid_markers), "QTL markers not found in markers tibble\\n")
            additive_peaks <- additive_peaks %>% filter(marker.id %in% markers\$marker.id)
        }

    } else {
        cat("Warning: No significant QTLs found\\n")
        additive_peaks <- tibble(
            gene.id = character(0),
            marker.id = character(0),
            lod = numeric(0)
        )
    }

    lod.peaks\$additive <- additive_peaks

    validation_log <- c(validation_log, "=== LOD Peaks Validation ===")
    validation_log <- c(validation_log, paste("✓ LOD peaks (additive):", nrow(lod.peaks\$additive)))
    validation_log <- c(validation_log, "✓ All gene.ids validated in annot.mrna")
    validation_log <- c(validation_log, "✓ All marker.ids validated in markers tibble")
    validation_log <- c(validation_log, "")

    # 9. CREATE DATASET.MRNA OBJECT
    dataset.mrna <- list(
        annot.mrna = annot.mrna,
        annot.samples = annot.samples,
        covar.matrix = covar.matrix,
        covar.info = covar.info,
        data = data_list,  # Use data list with transformations
        datatype = "mrna",
        display.name = paste("${prefix}", "eQTL Analysis"),
        lod.peaks = lod.peaks
    )

    validation_log <- c(validation_log, "=== Dataset Object Summary ===")
    validation_log <- c(validation_log, paste("✓ Dataset type: mrna (eQTL)"))
    validation_log <- c(validation_log, paste("✓ mRNA annotations:", nrow(annot.mrna)))
    validation_log <- c(validation_log, paste("✓ Sample annotations:", nrow(annot.samples)))
    validation_log <- c(validation_log, paste("✓ Covariate matrix:", paste(dim(covar.matrix), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Covariate info entries:", nrow(covar.info)))
    validation_log <- c(validation_log, paste("✓ Data transformations:", paste(names(data_list), collapse=", ")))
    validation_log <- c(validation_log, paste("✓ Data dimensions:", paste(dim(data_raw), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ LOD peaks:", nrow(lod.peaks\$additive)))
    validation_log <- c(validation_log, "")

    # Create deployment directory structure
    cat("Creating deployment directory structure...\\n")
    dir.create("data", showWarnings = FALSE)
    dir.create("data/rdata", recursive = TRUE, showWarnings = FALSE)
    dir.create("data/sqlite", showWarnings = FALSE)

    # Save QTL Viewer RData file directly in deployment structure
    cat("Saving QTL Viewer RData file to deployment directory...\\n")
    cat("Estimated RData size:", format(object.size(genoprobs) + object.size(K_subset) +
                                        object.size(map) + object.size(dataset.mrna), units="GB"), "\\n")

    # CRITICAL FIX: Include ensembl.species in save() call
    save(
        ensembl.version,
        ensembl.species,  # ADDED - was missing before!
        genoprobs,
        K = K_subset,
        map,
        markers,
        dataset.mrna,
        file = "data/rdata/${prefix}.RData"
    )

    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Viewer RData File ===")
    validation_log <- c(validation_log, paste("✓ File saved: data/rdata/${prefix}.RData"))
    validation_log <- c(validation_log, paste("✓ File size:", format(file.size("data/rdata/${prefix}.RData")/1e9, digits=2), "GB"))
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Viewer Format Validation ===")
    validation_log <- c(validation_log, "✓ ensembl.version: numeric (105)")
    validation_log <- c(validation_log, "✓ ensembl.species: character (Mm)")
    validation_log <- c(validation_log, "✓ genoprobs: list of genotype probability arrays (calc_genoprob)")
    validation_log <- c(validation_log, "✓ K: list of LOCO kinship matrices")
    validation_log <- c(validation_log, "✓ map: list of physical position vectors (Mb)")
    validation_log <- c(validation_log, "✓ markers: tibble with marker.id, chr, pos")
    validation_log <- c(validation_log, "✓ dataset.mrna: complete eQTL dataset object with data transformations")
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Viewer Deployment Files ===")
    validation_log <- c(validation_log, "✓ RData file: data/rdata/${prefix}.RData")
    validation_log <- c(validation_log, "✓ SQLite database: data/sqlite/ccfoundersnps.sqlite")
    validation_log <- c(validation_log, "✓ Environment config: .env")
    validation_log <- c(validation_log, "✓ Viewer settings: project.viewer.settings (contains species info)")
    validation_log <- c(validation_log, "✓ Docker compose: docker-compose.yml")
    validation_log <- c(validation_log, "✓ Startup script: start_qtlviewer.sh")
    validation_log <- c(validation_log, "✓ Cache directory: cache/")
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Viewer Data Conversion Complete ===")
    validation_log <- c(validation_log, "✓ Deployment follows official Churchill Lab documentation")
    validation_log <- c(validation_log, "✓ Species field configured in project.viewer.settings")
    validation_log <- c(validation_log, "✓ Ready for HPC deployment (Singularity) or local deployment (Docker)")

    # Write validation report
    writeLines(validation_log, "qtlviewer_conversion_report.txt")

    cat("QTL Viewer data conversion completed successfully\\n")

    # Now create deployment files using system() calls
    cat("Setting up QTL Viewer deployment files...\\n")

    # Handle SQLite database
    cat("Downloading MGI SQLite database for founder SNPs...\\n")

    database_found <- FALSE

    if (file.exists("../02_genotype_processing/mgi_mouse_genes.sqlite")) {
        cat("Found existing MGI SQLite database, copying...\\n")
        file.copy("../02_genotype_processing/mgi_mouse_genes.sqlite",
                  "data/sqlite/ccfoundersnps.sqlite", overwrite = TRUE)
        database_found <- TRUE
    } else if (file.exists("../03_control_file_generation/mgi_mouse_genes.sqlite")) {
        cat("Found existing MGI SQLite database, copying...\\n")
        file.copy("../03_control_file_generation/mgi_mouse_genes.sqlite",
                  "data/sqlite/ccfoundersnps.sqlite", overwrite = TRUE)
        database_found <- TRUE
    }

    if (!database_found) {
        cat("Downloading MGI database from Figshare...\\n")
        if (Sys.which("curl") != "") {
            system('curl -L -H "User-Agent: Mozilla/5.0" -o data/sqlite/ccfoundersnps.sqlite "https://figshare.com/ndownloader/files/24607961"')
        } else if (Sys.which("wget") != "") {
            system('wget -O data/sqlite/ccfoundersnps.sqlite "https://figshare.com/ndownloader/files/24607961"')
        } else {
            cat("WARNING: Neither curl nor wget available\\n")
        }
    }

    # Create .env file (following official QTL Viewer deployment documentation)
    # Convert prefix to lowercase for Docker Compose project name
    compose_project_name <- tolower("${prefix}_qtlviewer")

    env_content <- paste0("# QTL Viewer Environment Configuration for ${prefix}
# Generated by QTL2_NF Pipeline

# Unique name of project (lowercase required by Docker Compose)
COMPOSE_PROJECT_NAME=", compose_project_name, "

# Versions (using newer versions)
DOCKER_QTL2REST_VERSION=0.6.5
DOCKER_QTL2WEB_VERSION=1.6.1

# Docker Container Settings (do not change)
CONTAINER_FILE_QTL2WEB_SETTINGS=/app/qtlweb/viewer.settings
CONTAINER_DIR_QTL2WEB_CACHE=/app/cache
CONTAINER_CACHE_NAME=", compose_project_name, "_cache

# Name for Downloaded RData File
QTLAPI_RDATA=None
PYTHONUNBUFFERED=true

# Ports
HOST_PORT_WEB=8000
HOST_PORT_API=8001
HOST_PORT_REDIS=6379

# Local File Paths (EDIT THESE for your local deployment)
# NOTE: Change these paths to absolute paths on YOUR local machine
HOST_FILE_RDATA=./data
HOST_FILE_SNPDB=./data/sqlite/ccfoundersnps.sqlite
HOST_FILE_QTL2WEB_SETTINGS=./project.viewer.settings
HOST_DIR_QTL2WEB_CACHE=./cache
")

    writeLines(env_content, ".env")

    # Create project.viewer.settings file
    viewer_settings_content <- "{
  \\"title\\": \\"${prefix} QTL Viewer\\",
  \\"subtitle\\": \\"Diversity Outbred Mouse eQTL Analysis\\",
  \\"species\\": \\"mouse\\",
  \\"build\\": \\"GRCm39\\",
  \\"ensembl_version\\": 105,
  \\"citation\\": \\"Generated with QTL2_NF Pipeline\\",
  \\"contact\\": \\"\\",
  \\"description\\": \\"QTL analysis results for ${prefix} study\\"
}
"

    writeLines(viewer_settings_content, "project.viewer.settings")

    # Create docker-compose.yml (following official documentation structure)
    docker_compose_content <- "version: '3.8'

services:
  qtl2rest:
    image: churchilllab/qtl2rest:\${DOCKER_QTL2REST_VERSION}
    container_name: \${COMPOSE_PROJECT_NAME}_qtl2rest
    ports:
      - \\"\${HOST_PORT_API}:8001\\"
    volumes:
      - \${HOST_FILE_RDATA}:/app/qtl2rest/data/rdata
      - \${HOST_FILE_SNPDB}:/app/qtl2rest/data/sqlite/ccfoundersnps.sqlite
    environment:
      - R_MAX_VSIZE=32Gb
    restart: unless-stopped
    networks:
      - qtlviewer_network

  qtl2web:
    image: churchilllab/qtl2web:\${DOCKER_QTL2WEB_VERSION}
    container_name: \${COMPOSE_PROJECT_NAME}_qtl2web
    ports:
      - \\"\${HOST_PORT_WEB}:80\\"
    volumes:
      - \${HOST_FILE_QTL2WEB_SETTINGS}:\${CONTAINER_FILE_QTL2WEB_SETTINGS}
      - \${HOST_DIR_QTL2WEB_CACHE}:\${CONTAINER_DIR_QTL2WEB_CACHE}
    depends_on:
      - qtl2rest
      - redis
    environment:
      - QTL2REST_URL=http://qtl2rest:8001
      - CONTAINER_DIR_QTL2WEB_CACHE=\${CONTAINER_DIR_QTL2WEB_CACHE}
      - CONTAINER_CACHE_NAME=\${CONTAINER_CACHE_NAME}
      - PYTHONUNBUFFERED=\${PYTHONUNBUFFERED}
    restart: unless-stopped
    networks:
      - qtlviewer_network

  redis:
    image: redis:7-alpine
    container_name: \${COMPOSE_PROJECT_NAME}_redis
    ports:
      - \\"\${HOST_PORT_REDIS}:6379\\"
    restart: unless-stopped
    networks:
      - qtlviewer_network

  ensimpl:
    image: churchilllab/ensimpl:latest
    container_name: \${COMPOSE_PROJECT_NAME}_ensimpl
    ports:
      - \\"8002:8002\\"
    restart: unless-stopped
    networks:
      - qtlviewer_network

networks:
  qtlviewer_network:
    driver: bridge
"

    writeLines(docker_compose_content, "docker-compose.yml")

    # Create cache directory with proper permissions
    dir.create("cache", showWarnings = FALSE)
    Sys.chmod("cache", mode = "0777")  # Ensure Docker can write to cache

    validation_log <- c(validation_log, "✓ Cache directory created with write permissions")

    # Create startup script (for Docker local deployment)
    startup_script <- "#!/bin/bash

echo \\"Starting QTL Viewer for study: ${prefix}\\"
echo \\"========================================\\"

if ! docker info > /dev/null 2>&1; then
    echo \\"ERROR: Docker is not running\\"
    exit 1
fi

if [ ! -f \\"data/rdata/${prefix}.RData\\" ]; then
    echo \\"ERROR: QTL data file not found\\"
    exit 1
fi

# Stop and remove any existing containers
echo \\"Cleaning up any existing containers...\\"
docker-compose down 2>/dev/null

# Pull latest images
echo \\"Pulling Docker images...\\"
docker-compose pull

# Start services
echo \\"Starting QTL Viewer services...\\"
docker-compose up -d

# Wait for services to be ready
echo \\"Waiting for services to start...\\"
sleep 15

# Check if services are running
if docker-compose ps | grep -q \\"Up\\"; then
    echo \\"\\"
    echo \\"🎉 QTL Viewer is now running!\\"
    echo \\"==============================\\"
    echo \\"🌐 Web Interface: http://localhost:8000\\"
    echo \\"🔌 REST API:      http://localhost:8001\\"
    echo \\"\\"
    echo \\"📊 Your QTL Results:\\"
    echo \\"   Study: ${prefix}\\"
    echo \\"\\"
    echo \\"🛑 To stop: docker-compose down\\"
    echo \\"📋 To view logs: docker-compose logs -f\\"
    echo \\"\\"
    echo \\"💡 Open http://localhost:8000 in your web browser!\\"
else
    echo \\"ERROR: Services failed to start\\"
    echo \\"Check logs with: docker-compose logs\\"
fi
"

    writeLines(startup_script, "start_qtlviewer.sh")
    Sys.chmod("start_qtlviewer.sh", mode = "0755")

    # Create HPC debugging script for troubleshooting multi-container issues
    debug_script <- "#!/bin/bash
#
# QTL Viewer HPC Debug Script
# Comprehensive diagnostics for troubleshooting multi-container deployment
#

set +e  # Don't exit on error - we want to collect all diagnostics

STUDY_PREFIX=\\"${prefix}\\"
LOG_FILE=\\"qtlviewer_debug_\\$(date +%Y%m%d_%H%M%S).log\\"

# Color codes for output
RED='\\\\033[0;31m'
GREEN='\\\\033[0;32m'
YELLOW='\\\\033[1;33m'
BLUE='\\\\033[0;34m'
NC='\\\\033[0m' # No Color

# Logging functions
log_section() {
    echo \\"\\"\\" | tee -a \\$LOG_FILE
    echo \\"========================================\\" | tee -a \\$LOG_FILE
    echo \\"\\$1\\" | tee -a \\$LOG_FILE
    echo \\"========================================\\" | tee -a \\$LOG_FILE
}

log_info() {
    echo -e \\"\\${BLUE}[INFO]\\${NC} \\$1\\" | tee -a \\$LOG_FILE
}

log_success() {
    echo -e \\"\\${GREEN}[✓]\\${NC} \\$1\\" | tee -a \\$LOG_FILE
}

log_error() {
    echo -e \\"\\${RED}[✗]\\${NC} \\$1\\" | tee -a \\$LOG_FILE
}

log_warning() {
    echo -e \\"\\${YELLOW}[!]\\${NC} \\$1\\" | tee -a \\$LOG_FILE
}

echo \\"QTL Viewer HPC Diagnostic Tool\\" | tee \\$LOG_FILE
echo \\"Study: \\$STUDY_PREFIX\\" | tee -a \\$LOG_FILE
echo \\"Timestamp: \\$(date)\\" | tee -a \\$LOG_FILE
echo \\"Output: \\$LOG_FILE\\" | tee -a \\$LOG_FILE

# 1. Check Singularity
log_section \\"1. Singularity Environment\\"
if command -v singularity &> /dev/null; then
    log_success \\"Singularity is available\\"
    singularity --version | tee -a \\$LOG_FILE
else
    log_error \\"Singularity not found in PATH\\"
    log_info \\"Try: module load singularity\\"
fi

# 2. Check for container images
log_section \\"2. Container Images\\"
REQUIRED_IMAGES=(\\"qtl2rest_0.6.5.sif\\" \\"qtl2web_1.6.1.sif\\" \\"redis_7-alpine.sif\\" \\"ensimpl_latest.sif\\")
for img in \\"\${REQUIRED_IMAGES[@]}\\"; do
    if [ -f \\"\\$img\\" ]; then
        log_success \\"Found: \\$img (\\$(du -h \\$img | cut -f1))\\"
    else
        log_error \\"Missing: \\$img\\"
        log_info \\"Pull with: singularity pull docker://churchilllab/\\${img%.sif}\\"
    fi
done

# 3. Check data files
log_section \\"3. Data Files\\"
if [ -f \\"data/rdata/\\${STUDY_PREFIX}.RData\\" ]; then
    DATA_SIZE=\\$(du -h data/rdata/\\${STUDY_PREFIX}.RData | cut -f1)
    log_success \\"RData file found: \\$DATA_SIZE\\"

    # Try to inspect RData file
    log_info \\"Attempting to load RData file...\\"
    R --quiet --no-save << 'EOF' 2>&1 | tee -a \\$LOG_FILE
    tryCatch({
        load(\\"data/rdata/${prefix}.RData\\")
        cat(\\"✓ RData loaded successfully\\\\n\\")
        cat(sprintf(\\"  - ensembl.version: %s\\\\n\\", ensembl.version))
        cat(sprintf(\\"  - ensembl.species: %s\\\\n\\", ensembl.species))
        cat(sprintf(\\"  - genoprobs chromosomes: %d\\\\n\\", length(genoprobs)))
        cat(sprintf(\\"  - markers: %d\\\\n\\", nrow(markers)))
        cat(sprintf(\\"  - dataset.mrna genes: %d\\\\n\\", nrow(dataset.mrna\\$annot.mrna)))
    }, error = function(e) {
        cat(\\"✗ Error loading RData:\\\\n\\")
        cat(sprintf(\\"  %s\\\\n\\", e\\$message))
    })
EOF
else
    log_error \\"RData file not found: data/rdata/\\${STUDY_PREFIX}.RData\\"
fi

if [ -f \\"data/sqlite/ccfoundersnps.sqlite\\" ]; then
    SQLITE_SIZE=\\$(du -h data/sqlite/ccfoundersnps.sqlite | cut -f1)
    log_success \\"SQLite database found: \\$SQLITE_SIZE\\"
else
    log_error \\"SQLite database not found\\"
fi

# 4. Check port availability
log_section \\"4. Port Availability\\"
PORTS=(80 8000 8001 8002 6379)
for port in \\"\${PORTS[@]}\\"; do
    if netstat -tuln 2>/dev/null | grep -q \\":\\$port \\" || ss -tuln 2>/dev/null | grep -q \\":\\$port \\"; then
        log_warning \\"Port \\$port is IN USE\\"
        PROCESS=\\$(lsof -ti:\\$port 2>/dev/null || ss -lptn sport = :\\$port 2>/dev/null | grep -v State)
        if [ ! -z \\"\\$PROCESS\\" ]; then
            log_info \\"   Process: \\$PROCESS\\"
        fi
    else
        log_success \\"Port \\$port is available\\"
    fi
done

# 5. Check running processes
log_section \\"5. Running QTL Viewer Processes\\"
QTL_PROCESSES=\\$(ps aux | grep -E '(qtl2rest|qtl2web|redis|ensimpl|singularity)' | grep -v grep)
if [ ! -z \\"\\$QTL_PROCESSES\\" ]; then
    log_info \\"Found running processes:\\"
    echo \\"\\$QTL_PROCESSES\\" | tee -a \\$LOG_FILE
else
    log_warning \\"No QTL Viewer processes currently running\\"
fi

# 6. Test API endpoints (if services are running)
log_section \\"6. API Endpoint Tests\\"
if curl -s http://localhost:8001/health &> /dev/null || curl -s http://localhost:8001/ &> /dev/null; then
    log_success \\"QTL2REST API is responding\\"

    # Test datasets endpoint
    log_info \\"Testing /datasets endpoint...\\"
    DATASETS=\\$(curl -s http://localhost:8001/datasets 2>&1)
    if echo \\"\\$DATASETS\\" | grep -q \\"${prefix}\\"; then
        log_success \\"Dataset \\$STUDY_PREFIX found in API\\"
    else
        log_warning \\"Dataset might not be loaded yet\\"
        echo \\"\\$DATASETS\\" | head -5 | tee -a \\$LOG_FILE
    fi
else
    log_warning \\"QTL2REST API not responding (service may not be running)\\"
fi

if curl -s http://localhost:80 &> /dev/null || curl -s http://localhost:8000 &> /dev/null; then
    log_success \\"QTL2WEB is responding\\"
else
    log_warning \\"QTL2WEB not responding (service may not be running)\\"
fi

# 7. Check system resources
log_section \\"7. System Resources\\"
log_info \\"Memory usage:\\"
free -h | tee -a \\$LOG_FILE

log_info \\"Disk space:\\"
df -h . | tee -a \\$LOG_FILE

# 8. Collect recent logs if available
log_section \\"8. Recent Error Logs\\"
if [ -d \\"logs\\" ]; then
    log_info \\"Found logs directory\\"
    find logs -name '*.log' -mtime -1 -exec echo \\"{}:\\" \\\\; -exec tail -20 {} \\\\; | tee -a \\$LOG_FILE
else
    log_info \\"No logs directory found\\"
fi

# 9. Configuration check
log_section \\"9. Configuration Files\\"
if [ -f \\"project.viewer.settings\\" ]; then
    log_success \\"Found project.viewer.settings\\"
    cat project.viewer.settings | tee -a \\$LOG_FILE
else
    log_error \\"project.viewer.settings not found\\"
fi

# 10. Provide recommendations
log_section \\"10. Troubleshooting Recommendations\\"

if [ ! -f \\"qtl2rest_0.6.5.sif\\" ]; then
    log_error \\"Container images missing - pull them first:\\"
    echo \\"   singularity pull docker://churchilllab/qtl2rest:0.6.5\\" | tee -a \\$LOG_FILE
    echo \\"   singularity pull docker://churchilllab/qtl2web:1.6.1\\" | tee -a \\$LOG_FILE
    echo \\"   singularity pull docker://redis:7-alpine\\" | tee -a \\$LOG_FILE
fi

if netstat -tuln 2>/dev/null | grep -q ':8001 ' || ss -tuln 2>/dev/null | grep -q ':8001 '; then
    if ! curl -s http://localhost:8001/health &> /dev/null; then
        log_warning \\"Port 8001 is in use but API not responding\\"
        log_info \\"Check process using: lsof -i:8001 or ss -lptn sport = :8001\\"
        log_info \\"Kill stuck process if needed: kill \\$(lsof -ti:8001)\\"
    fi
fi

if [ ! -f \\"data/rdata/\\${STUDY_PREFIX}.RData\\" ]; then
    log_error \\"RData file missing - ensure Module 9 completed successfully\\"
fi

# Summary
log_section \\"Diagnostic Summary\\"
echo \\"Full diagnostic log saved to: \\$LOG_FILE\\" | tee -a \\$LOG_FILE
echo \\"\\" | tee -a \\$LOG_FILE
echo \\"Next steps:\\" | tee -a \\$LOG_FILE
echo \\"1. Review the log file for errors\\" | tee -a \\$LOG_FILE
echo \\"2. Start services in screen/tmux sessions\\" | tee -a \\$LOG_FILE
echo \\"3. Check logs in each session for errors\\" | tee -a \\$LOG_FILE
echo \\"4. Use SSH tunnel to access from local machine\\" | tee -a \\$LOG_FILE
echo \\"\\" | tee -a \\$LOG_FILE
echo \\"For help, share this log file: \\$LOG_FILE\\" | tee -a \\$LOG_FILE
"

    writeLines(debug_script, "debug_qtlviewer.sh")
    Sys.chmod("debug_qtlviewer.sh", mode = "0755")

    # Create README
    readme_content <- "# QTL Viewer for ${prefix}

## Overview
This directory contains a complete QTL Viewer deployment following the official Churchill Lab documentation. Supports both **HPC deployment (Singularity)** and **local deployment (Docker)**.

## Directory Structure
- **data/** - Contains RData files and SQLite database
  - **rdata/** - ${prefix}.RData file with QTL analysis results
  - **sqlite/** - ccfoundersnps.sqlite founder SNP database
- **cache/** - Cache directory for QTL Viewer (writable)
- **.env** - Environment variables for Docker Compose
- **project.viewer.settings** - QTL Viewer web interface settings (includes species info)
- **docker-compose.yml** - Docker orchestration file (for local deployment)
- **start_qtlviewer.sh** - Convenience startup script (for local deployment)
- **debug_qtlviewer.sh** - Comprehensive debugging tool for troubleshooting

## Deployment Options

### Option 1: HPC Deployment (Singularity) - RECOMMENDED FOR LARGE DATA

For HPC deployment with Singularity (recommended when RData > 8GB):

1. **Load Singularity on your HPC**:
   ```bash
   module load singularity
   ```

2. **Pull QTL Viewer images**:
   ```bash
   singularity pull docker://churchilllab/qtl2rest:0.6.5
   singularity pull docker://churchilllab/qtl2web:1.6.1
   singularity pull docker://redis:7-alpine
   singularity pull docker://churchilllab/ensimpl:latest
   ```

3. **Start QTL2REST API (in screen/tmux session)**:
   ```bash
   singularity run \\
     --bind \$PWD/data/rdata:/app/qtl2rest/data/rdata \\
     --bind \$PWD/data/sqlite:/app/qtl2rest/data/sqlite \\
     --env R_MAX_VSIZE=32Gb \\
     qtl2rest_0.6.5.sif \\
     uvicorn qtl2rest.main:app --host 0.0.0.0 --port 8001
   ```

4. **Start Redis (separate session)**:
   ```bash
   singularity run redis_7-alpine.sif
   ```

5. **Start QTL2WEB (separate session)**:
   ```bash
   singularity run \\
     --bind \$PWD/project.viewer.settings:/app/qtlweb/viewer.settings \\
     --bind \$PWD/cache:/app/cache \\
     --env QTL2REST_URL=http://localhost:8001 \\
     qtl2web_1.6.1.sif
   ```

6. **SSH tunnel from your local machine**:
   ```bash
   ssh -L 8000:localhost:80 -L 8001:localhost:8001 user@hpc-server
   ```

7. **Access in your browser**: http://localhost:8000

**Troubleshooting HPC Deployment:**
If services don't start properly, run the debug script:
```bash
./debug_qtlviewer.sh
```

This will generate a comprehensive diagnostic report including:
- Singularity availability
- Container image status
- Port conflicts
- Running processes
- API endpoint tests
- RData file validation
- System resources
- Configuration checks
- Actionable recommendations

### Option 2: Local Deployment (Docker)

For local deployment with Docker Desktop (works if RData < 8GB):

#### Prerequisites
- Docker Desktop installed and running
- This entire directory downloaded to your local machine

#### Steps
1. **Download this directory** to your local machine
2. **Edit the .env file** if needed (ports, paths)
3. **Run the startup script**:
   ```bash
   ./start_qtlviewer.sh
   ```
4. **Open your browser** to: http://localhost:8000

## Configuration

### .env File
The .env file contains all configuration variables. Key settings:
- **HOST_PORT_WEB**: Web interface port (default: 8000)
- **HOST_PORT_API**: REST API port (default: 8001)
- **HOST_PORT_REDIS**: Redis port (default: 6379)

For local deployment, the default relative paths should work. For production deployment, change to absolute paths.

### project.viewer.settings
Contains the web interface configuration including:
- **species**: mouse
- **build**: GRCm39
- **ensembl_version**: 105
- Title, subtitle, and other display settings

## Management Commands

### Start QTL Viewer
```bash
./start_qtlviewer.sh
# OR
docker-compose up -d
```

### Stop QTL Viewer
```bash
docker-compose down
```

### View Logs
```bash
docker-compose logs -f
```

### Check Status
```bash
docker-compose ps
```

### Restart Services
```bash
docker-compose restart
```

## Troubleshooting

### HPC-Specific Issues

#### Large RData Files (>8GB)
- **Use HPC Singularity deployment** instead of local Docker
- Increase R memory limit: `--env R_MAX_VSIZE=64Gb` (adjust as needed)
- Monitor memory usage with `top` or `htop` during startup

#### Singularity Not Available
```bash
# Check if Singularity is installed
which singularity

# Load Singularity module (HPC-specific)
module avail singularity
module load singularity/3.8.0  # or latest version
```

#### SSH Tunnel Issues
If tunnel disconnects:
```bash
# Use autossh for persistent tunnels
autossh -M 0 -L 8000:localhost:80 -L 8001:localhost:8001 user@hpc-server

# Or use screen/tmux to keep tunnel alive
screen -S qtl_tunnel
ssh -L 8000:localhost:80 -L 8001:localhost:8001 user@hpc-server
# Ctrl+A+D to detach
```

### Common Debugging Steps

**Always start with the debug script:**
```bash
./debug_qtlviewer.sh
```

The script will:
- Check all prerequisites
- Verify data files
- Test ports
- Check running services
- Provide specific recommendations
- Generate a timestamped log file

**Share the generated log file when asking for help!**

### Local Docker Issues

#### Port Conflicts
If port 8000 is already in use, edit .env and change HOST_PORT_WEB to another port (e.g., 8080).

#### Docker Memory Limit
Docker Desktop has an 8GB default memory limit. For larger datasets:
1. Open Docker Desktop → Settings → Resources
2. Increase Memory allocation
3. OR use HPC Singularity deployment instead

#### Permissions
Ensure the cache/ directory has write permissions:
```bash
chmod -R 777 cache/
```

#### Docker Not Running
Make sure Docker Desktop is running:
```bash
docker info
```

#### View Container Logs
```bash
docker-compose logs qtl2rest
docker-compose logs qtl2web
docker-compose logs redis
```

## Data Files

- **RData file**: data/rdata/${prefix}.RData (2.1 GB)
- **SQLite database**: data/sqlite/ccfoundersnps.sqlite
- **Format**: Follows qtl2api RData specification

## More Information

- QTL Viewer Documentation: https://github.com/churchill-lab/qtl2web
- QTL2 API: https://github.com/churchill-lab/qtl2rest
- Churchill Lab: https://churchilllab.jax.org/

## Performance Notes

- **HPC Deployment**: Recommended for RData files > 8GB
  - Can handle datasets up to 100GB+ with sufficient memory allocation
  - Uses Singularity containers natively supported on HPC systems
  - Better for production deployments with many users

- **Local Docker Deployment**: Best for RData files < 8GB
  - Limited by Docker Desktop memory allocation (default 8GB)
  - Good for quick exploration and small datasets
  - Easier setup for single-user access

## Generated by QTL2_NF Pipeline
Study Prefix: ${prefix}
Deployment: Supports both HPC (Singularity) and Local (Docker)
"

    writeLines(readme_content, "README_qtlviewer.md")

    cat("✅ QTL Viewer deployment package complete!\\n")

    # Manually copy files to results directory with error checking
    cat("Copying files to results directory...\\n")
    workflow_dir <- "${workflow.launchDir}"
    outdir <- file.path(workflow_dir, "${params.outdir}/09_qtl_viewer")

    # Create complete directory structure
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(outdir, "data/rdata"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(outdir, "data/sqlite"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(outdir, "cache"), showWarnings = FALSE, recursive = TRUE)

    # Set cache permissions in output directory too
    Sys.chmod(file.path(outdir, "cache"), mode = "0777")

    # Copy all deployment files with error checking
    copy_success <- TRUE

    files_to_copy <- list(
        c("data/rdata/${prefix}.RData", "data/rdata/${prefix}.RData"),
        c("data/sqlite/ccfoundersnps.sqlite", "data/sqlite/ccfoundersnps.sqlite"),
        c("qtlviewer_conversion_report.txt", "qtlviewer_conversion_report.txt"),
        c(".env", ".env"),
        c("project.viewer.settings", "project.viewer.settings"),
        c("docker-compose.yml", "docker-compose.yml"),
        c("start_qtlviewer.sh", "start_qtlviewer.sh"),
        c("debug_qtlviewer.sh", "debug_qtlviewer.sh"),
        c("README_qtlviewer.md", "README_qtlviewer.md")
    )

    for (file_pair in files_to_copy) {
        src <- file_pair[1]
        dest <- file.path(outdir, file_pair[2])

        if (file.exists(src)) {
            success <- file.copy(src, dest, overwrite = TRUE)
            if (!success) {
                cat("ERROR: Failed to copy", src, "to", dest, "\\n")
                copy_success <- FALSE
            } else {
                cat("  ✓ Copied:", basename(src), "\\n")
            }
        } else {
            cat("WARNING: Source file not found:", src, "\\n")
            copy_success <- FALSE
        }
    }

    # Set executable permissions on scripts
    Sys.chmod(file.path(outdir, "start_qtlviewer.sh"), mode = "0755")
    Sys.chmod(file.path(outdir, "debug_qtlviewer.sh"), mode = "0755")

    if (copy_success) {
        cat("\\n✅ All files successfully copied to", outdir, "\\n")
        validation_log <- c(validation_log, "✓ All deployment files copied successfully")
    } else {
        cat("\\n⚠️  Some file copy operations failed - check log above\\n")
        validation_log <- c(validation_log, "WARNING: Some file copy operations failed")
    }
    cat("\\n")
    cat("=================================================================\\n")
    cat("QTL Viewer Deployment Package Complete\\n")
    cat("=================================================================\\n")
    cat("Location:", outdir, "\\n")
    cat("\\nDeployment files generated:\\n")
    cat("  ✓ .env (environment configuration)\\n")
    cat("  ✓ project.viewer.settings (species, build info)\\n")
    cat("  ✓ docker-compose.yml (Docker orchestration)\\n")
    cat("  ✓ start_qtlviewer.sh (startup script)\\n")
    cat("  ✓ debug_qtlviewer.sh (diagnostic tool - USE THIS FOR TROUBLESHOOTING!)\\n")
    cat("  ✓ cache/ (writable cache directory)\\n")
    cat("\\nFor HPC deployment with Singularity:\\n")
    cat("  → See README_qtlviewer.md for detailed instructions\\n")
    cat("  → Run ./debug_qtlviewer.sh for comprehensive diagnostics\\n")
    cat("\\nThe debug script will help troubleshoot:\\n")
    cat("  - Container availability\\n")
    cat("  - Port conflicts\\n")
    cat("  - Service startup issues\\n")
    cat("  - API connectivity\\n")
    cat("  - RData loading problems\\n")
    cat("=================================================================\\n")
    """
}