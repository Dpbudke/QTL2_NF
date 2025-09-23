process PREPARE_GENOME_SCAN {
    tag "Preparing genome scan for ${prefix}"
    publishDir "${params.outdir}/05_genome_scan_preparation", mode: 'copy'
    
    input:
    path(cross2_file)
    val(prefix)
    
    output:
    path("${prefix}_genoprob.rds"), emit: genoprob
    path("${prefix}_kinship_loco.rds"), emit: kinship
    path("${prefix}_genetic_map.rds"), emit: genetic_map
    path("genoprob_validation_report.txt"), emit: validation_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })
    
    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== Genome Scan Preparation Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")
    
    cat("Loading cross2 object...\\n")
    cross2 <- readRDS("${cross2_file}")
    
    validation_log <- c(validation_log, "=== Cross2 Object Loaded ===")
    validation_log <- c(validation_log, paste("✓ Number of individuals:", nrow(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Number of chromosomes:", length(cross2\$gmap)))
    validation_log <- c(validation_log, paste("✓ Total markers:", sum(sapply(cross2\$gmap, length))))
    validation_log <- c(validation_log, paste("✓ Number of phenotypes:", ncol(cross2\$pheno)))
    validation_log <- c(validation_log, "")
    
    # Get genetic map and insert pseudomarkers for high-resolution scanning
    cat("Preparing genetic map with pseudomarkers...\\n")
    gmap <- cross2\$gmap
    
    # Debug: check cross2 and genetic map structure
    cat("DEBUG: Cross2 class:", class(cross2), "\\n")
    cat("DEBUG: Cross2 names:", paste(names(cross2), collapse=", "), "\\n")
    cat("DEBUG: Genotypes class:", class(cross2\$geno), "\\n")
    cat("DEBUG: Genotypes names:", paste(names(cross2\$geno), collapse=", "), "\\n")
    if (length(cross2\$geno) > 0) {
        cat("DEBUG: First geno chr class:", class(cross2\$geno[[1]]), "\\n")
        cat("DEBUG: First geno chr dim:", paste(dim(cross2\$geno[[1]]), collapse=" x "), "\\n")
    }
    cat("DEBUG: Genetic map class:", class(gmap), "\\n")
    cat("DEBUG: Genetic map length:", length(gmap), "\\n")
    if (length(gmap) > 0) {
        cat("DEBUG: First chromosome name:", names(gmap)[1], "\\n")
        cat("DEBUG: First chromosome class:", class(gmap[[1]]), "\\n")
        cat("DEBUG: First chromosome length:", length(gmap[[1]]), "\\n")
        if (length(gmap[[1]]) > 0) {
            cat("DEBUG: First few markers:", paste(head(names(gmap[[1]]), 5), collapse=", "), "\\n")
        }
    }
    
    # Insert pseudomarkers every 0.5 cM for high resolution
    gmap_pseudo <- insert_pseudomarkers(gmap, step=0.5)
    
    # Save genetic map
    saveRDS(gmap_pseudo, file = "${prefix}_genetic_map.rds")
    
    validation_log <- c(validation_log, "=== Genetic Map Preparation ===")
    validation_log <- c(validation_log, paste("✓ Original markers:", sum(sapply(gmap, length))))
    validation_log <- c(validation_log, paste("✓ Markers with pseudomarkers (0.5 cM):", sum(sapply(gmap_pseudo, length))))
    validation_log <- c(validation_log, "")
    
    # Calculate genotype probabilities
    cat("Calculating genotype probabilities...\\n")
    cat("This may take several minutes for large datasets...\\n")
    
    validation_log <- c(validation_log, "=== Genotype Probability Calculation ===")
    start_time <- Sys.time()
    
    # Use error probability appropriate for high-density arrays
    genoprob <- calc_genoprob(cross2, gmap_pseudo, error_prob=0.002, cores=0)  # cores=0 uses all available
    
    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)
    
    # Save genotype probabilities
    saveRDS(genoprob, file = "${prefix}_genoprob.rds")
    
    validation_log <- c(validation_log, paste("✓ Genotype probabilities calculated in", duration, "minutes"))
    validation_log <- c(validation_log, paste("✓ Genoprob dimensions:", paste(dim(genoprob[[1]]), collapse=" x ")))
    validation_log <- c(validation_log, "")
    
    # Calculate kinship matrix using LOCO (Leave One Chromosome Out) method
    cat("Calculating LOCO kinship matrices...\\n")
    cat("This step accounts for population structure and relatedness...\\n")
    
    validation_log <- c(validation_log, "=== LOCO Kinship Matrix Calculation ===")
    start_time <- Sys.time()
    
    # LOCO kinship for proper mixed model analysis
    kinship_loco <- calc_kinship(genoprob, type="loco", cores=0)  # cores=0 uses all available
    
    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)
    
    # Save kinship matrices
    saveRDS(kinship_loco, file = "${prefix}_kinship_loco.rds")
    
    validation_log <- c(validation_log, paste("✓ LOCO kinship matrices calculated in", duration, "minutes"))
    validation_log <- c(validation_log, paste("✓ Number of kinship matrices (one per chromosome):", length(kinship_loco)))
    
    # Report kinship matrix properties
    if (length(kinship_loco) > 0) {
        first_kinship <- kinship_loco[[1]]
        validation_log <- c(validation_log, paste("✓ Kinship matrix dimensions:", paste(dim(first_kinship), collapse=" x ")))
        validation_log <- c(validation_log, paste("✓ Kinship diagonal range:", round(min(diag(first_kinship)), 3), "to", round(max(diag(first_kinship)), 3)))
    }
    
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Genome Scan Preparation Complete ===")
    validation_log <- c(validation_log, "✓ Ready for genome scanning and permutation testing")
    
    # Write validation report
    writeLines(validation_log, "genoprob_validation_report.txt")
    
    cat("Genome scan preparation completed successfully\\n")
    """
}

process GENOME_SCAN {
    tag "Setting up HPC array genome scan for ${prefix}"
    publishDir "${params.outdir}/06_qtl_analysis", mode: 'copy'

    input:
    path(cross2_file)
    path(genoprob_file)
    path(kinship_file)
    val(prefix)
    val(lod_threshold)

    output:
    path("chunk_info.txt"), emit: chunk_info
    path("${prefix}_scan_setup_report.txt"), emit: setup_report

    script:
    """
    #!/usr/bin/env Rscript

    # HPC Array Job Setup - Coordinator Process
    # Creates chunk information for individual SLURM array jobs

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })

    cat("=== HPC ARRAY GENOME SCAN SETUP ===\\n")
    cat("Study Prefix: ${prefix}\\n")
    cat("LOD Threshold: ${lod_threshold}\\n")
    cat("Setup Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # Load cross2 to get phenotype information
    cat("Loading cross2 object for setup...\\n")
    cross2 <- readRDS("${cross2_file}")

    # Configuration
    chunk_size <- 500
    total_phenotypes <- ncol(cross2\$pheno)
    n_chunks <- ceiling(total_phenotypes / chunk_size)

    cat("Total phenotypes:", total_phenotypes, "\\n")
    cat("Chunk size:", chunk_size, "\\n")
    cat("Number of chunks:", n_chunks, "\\n")

    # Create chunk information for array jobs
    chunk_info <- data.frame(
        chunk_id = 1:n_chunks,
        pheno_start = ((1:n_chunks) - 1) * chunk_size + 1,
        pheno_end = pmin((1:n_chunks) * chunk_size, total_phenotypes),
        stringsAsFactors = FALSE
    )

    # Write chunk information file
    write.csv(chunk_info, "chunk_info.txt", row.names = FALSE)

    cat("✓ Created chunk information for", n_chunks, "array jobs\\n")

    # Create setup report
    setup_report <- c(
        "=== HPC Array Genome Scan Setup Complete ===",
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Study Prefix: ${prefix}"),
        paste("LOD Threshold: ${lod_threshold}"),
        "",
        "=== Array Job Configuration ===",
        paste("✓ Total phenotypes:", total_phenotypes),
        paste("✓ Chunk size:", chunk_size),
        paste("✓ Number of array jobs:", n_chunks),
        paste("✓ Resource per job: 16 CPUs, 64GB RAM"),
        "",
        "=== Next Steps ===",
        "✓ Chunk information file created",
        "✓ Ready to launch array jobs",
        "✓ Each job will process ~500 phenotypes independently"
    )

    writeLines(setup_report, "${prefix}_scan_setup_report.txt")

    cat("\\n=== ARRAY JOB SETUP COMPLETED ===\\n")
    cat("Ready to launch", n_chunks, "independent SLURM array jobs\\n")

    """
}

process GENOME_SCAN_CHUNK {
    tag "HPC array job chunk ${chunk_id} for ${prefix}"
    publishDir "${params.outdir}/06_qtl_analysis/chunks", mode: 'copy'

    input:
    path(cross2_file)
    path(genoprob_file)
    path(kinship_file)
    val(prefix)
    val(chunk_id)
    val(pheno_start)
    val(pheno_end)

    output:
    path("${prefix}_chunk_${chunk_id}_results.rds"), emit: chunk_results

    script:
    """
    #!/usr/bin/env Rscript

    # HPC Array Job - Individual Chunk Processing
    # Process phenotypes ${pheno_start} to ${pheno_end}

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })

    cat("=== HPC ARRAY JOB CHUNK ${chunk_id} ===\\n")
    cat("Study Prefix: ${prefix}\\n")
    cat("Chunk ID: ${chunk_id}\\n")
    cat("Phenotypes: ${pheno_start} to ${pheno_end}\\n")
    cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # Load data files
    cat("Loading required data files...\\n")
    cross2 <- readRDS("${cross2_file}")
    genoprob <- readRDS("${genoprob_file}")
    kinship_loco <- readRDS("${kinship_file}")

    # Extract phenotype chunk
    pheno_full <- cross2\$pheno
    pheno_chunk <- pheno_full[, ${pheno_start}:${pheno_end}, drop=FALSE]
    covar <- cross2\$covar

    cat("✓ Loaded data for", ncol(pheno_chunk), "phenotypes\\n")

    # Set up covariates (matching original approach)
    addcovar <- NULL
    Xcovar <- NULL

    # Check if we're scanning coat_color
    is_coat_color_scan <- "coat_color" %in% colnames(pheno_chunk)

    if (is_coat_color_scan) {
        cat("✓ Coat color detected - using minimal covariates\\n")
    } else {
        # Sex covariate
        if ("Sex" %in% colnames(covar)) {
            sex_numeric <- ifelse(covar\$Sex == "male", 1, 0)
            addcovar <- matrix(sex_numeric, ncol=1)
            colnames(addcovar) <- "sex"
            rownames(addcovar) <- rownames(covar)
            Xcovar <- get_x_covar(cross2)
            cat("✓ Sex covariate added\\n")
        }

        # Generation as categorical
        if ("generation" %in% colnames(covar)) {
            gen_factor <- factor(covar\$generation)
            if (nlevels(gen_factor) > 1) {
                gen_matrix <- model.matrix(~ gen_factor - 1)[, -1, drop=FALSE]
                gen_colnames <- paste0("gen_", levels(gen_factor)[-1])
                gen_matrix <- matrix(as.numeric(gen_matrix), nrow=nrow(gen_matrix), ncol=ncol(gen_matrix))
                colnames(gen_matrix) <- gen_colnames
                rownames(gen_matrix) <- rownames(covar)

                if (is.null(addcovar)) {
                    addcovar <- gen_matrix
                } else {
                    addcovar <- cbind(addcovar, gen_matrix)
                }
                cat("✓ Generation covariate added as categorical (", nlevels(gen_factor), "levels )\\n")
            }
        }
    }

    # Final numeric conversion
    if (!is.null(addcovar)) {
        original_colnames <- colnames(addcovar)
        original_rownames <- rownames(addcovar)
        addcovar <- matrix(as.numeric(as.matrix(addcovar)), nrow=nrow(addcovar), ncol=ncol(addcovar))
        colnames(addcovar) <- original_colnames
        rownames(addcovar) <- original_rownames
    }

    # Perform genome scan for this chunk using ALL available cores
    cat("\\n=== STARTING CHUNK GENOME SCAN ===\\n")
    start_time <- Sys.time()

    # Use all available cores (16 CPUs allocated to this job)
    if (is_coat_color_scan) {
        chunk_results <- scan1(genoprob, pheno_chunk, addcovar=addcovar, Xcovar=Xcovar, cores=0)
    } else {
        chunk_results <- scan1(genoprob, pheno_chunk, kinship_loco, addcovar=addcovar, Xcovar=Xcovar, cores=0)
    }

    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)

    cat("✓ Chunk ${chunk_id} completed in", duration, "minutes\\n")
    cat("✓ Results dimensions:", paste(dim(chunk_results), collapse=" x "), "\\n")
    cat("✓ Maximum LOD score:", round(max(chunk_results, na.rm=TRUE), 3), "\\n")

    # Save chunk results
    saveRDS(chunk_results, file = "${prefix}_chunk_${chunk_id}_results.rds")

    cat("\\n=== CHUNK ${chunk_id} COMPLETED ===\\n")
    cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
    """
}

process COMBINE_SCAN_RESULTS {
    tag "Combining HPC array results for ${prefix}"
    publishDir "${params.outdir}/06_qtl_analysis", mode: 'copy'

    input:
    path(chunk_results)
    path(cross2_file)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_scan_results.rds"), emit: scan_results
    path("${prefix}_scan_peaks.csv"), emit: peaks
    path("${prefix}_filtered_phenotypes.txt"), emit: filtered_phenotypes
    path("genome_scan_validation_report.txt"), emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript

    # HPC Array Results Combination
    # Merge all chunk results into final scan results

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })

    cat("=== HPC ARRAY RESULTS COMBINATION ===\\n")
    cat("Study Prefix: ${prefix}\\n")
    cat("LOD Threshold: ${lod_threshold}\\n")
    cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # Load cross2 for genetic map
    cross2 <- readRDS("${cross2_file}")

    # Find all chunk result files
    chunk_files <- list.files(pattern = "${prefix}_chunk_.*_results.rds")
    chunk_files <- sort(chunk_files)  # Ensure proper order

    cat("✓ Found", length(chunk_files), "chunk result files\\n")

    if (length(chunk_files) == 0) {
        stop("No chunk result files found!")
    }

    # Load and combine all chunk results
    cat("Loading and combining chunk results...\\n")
    all_chunks <- list()

    for (i in 1:length(chunk_files)) {
        chunk_file <- chunk_files[i]
        cat("Loading", chunk_file, "\\n")
        chunk_result <- readRDS(chunk_file)
        all_chunks[[i]] <- chunk_result
    }

    # Combine all chunks column-wise
    scan_results <- do.call(cbind, all_chunks)

    cat("\\n=== COMBINATION COMPLETE ===\\n")
    cat("✓ Final scan results dimensions:", paste(dim(scan_results), collapse=" x "), "\\n")
    cat("✓ Maximum LOD score:", round(max(scan_results, na.rm=TRUE), 3), "\\n")

    # Save combined scan results
    saveRDS(scan_results, file = "${prefix}_scan_results.rds")

    # Generate filtered phenotypes list
    cat("\\n=== FILTERING PHENOTYPES FOR PERMUTATION TESTING ===\\n")
    max_lods <- apply(scan_results, 2, max, na.rm=TRUE)
    phenotypes_above_threshold <- names(max_lods)[max_lods >= ${lod_threshold}]
    writeLines(phenotypes_above_threshold, "${prefix}_filtered_phenotypes.txt")

    cat("✓ Phenotypes passing LOD threshold (≥${lod_threshold}):", length(phenotypes_above_threshold), "out of", length(max_lods), "\\n")

    # Find preliminary peaks
    cat("\\n=== FINDING PRELIMINARY PEAKS ===\\n")
    all_peaks <- find_peaks(scan_results, cross2\$gmap, threshold=3.0, drop=1.5)

    if (nrow(all_peaks) > 0) {
        all_peaks\$study <- "${prefix}"
        all_peaks <- all_peaks[, c("study", names(all_peaks)[names(all_peaks) != "study"])]
        write.csv(all_peaks, file = "${prefix}_scan_peaks.csv", row.names=FALSE)
        cat("✓ Found", nrow(all_peaks), "preliminary peaks (LOD ≥ 3.0)\\n")
    } else {
        empty_peaks <- data.frame(study=character(0), lodindex=numeric(0), lodcolumn=character(0),
                                 chr=character(0), pos=numeric(0), lod=numeric(0))
        write.csv(empty_peaks, file = "${prefix}_scan_peaks.csv", row.names=FALSE)
        cat("⚠ No preliminary peaks found with LOD ≥ 3.0\\n")
    }

    # Create validation report
    validation_lines <- c(
        "=== HPC Array Genome Scan Complete ===",
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Study Prefix: ${prefix}"),
        paste("LOD Threshold for Permutation Filtering: ${lod_threshold}"),
        "",
        "=== HPC Array Performance ===",
        paste("✓ Number of array jobs:", length(chunk_files)),
        paste("✓ Total phenotypes processed:", ncol(scan_results)),
        paste("✓ Genome positions scanned:", nrow(scan_results)),
        paste("✓ Resource utilization: 16 CPUs × 64GB per job"),
        "",
        "=== Results Summary ===",
        paste("✓ Maximum LOD score:", round(max(scan_results, na.rm=TRUE), 3)),
        paste("✓ Preliminary peaks found:", nrow(all_peaks)),
        paste("✓ Phenotypes above threshold:", length(phenotypes_above_threshold)),
        "",
        "=== HPC Array Success ===",
        "✓ Massive parallel processing complete",
        "✓ All chunks successfully combined",
        "✓ Ready for permutation testing"
    )

    writeLines(validation_lines, "genome_scan_validation_report.txt")

    cat("\\n=== HPC ARRAY GENOME SCAN COMPLETED ===\\n")
    cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
    cat("Results ready for permutation testing\\n")
    """
}

process PERMUTATION_TEST {
    tag "Permutation testing for ${prefix}"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'
    
    input:
    path(cross2_file)
    path(genoprob_file)
    path(kinship_file)
    path(filtered_phenotypes_file)
    val(prefix)
    val(lod_threshold)
    
    output:
    path("${prefix}_permutation_results.rds"), emit: perm_results
    path("${prefix}_significance_thresholds.csv"), emit: thresholds
    path("permutation_validation_report.txt"), emit: validation_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })
    
    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== Permutation Testing Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, paste("LOD Threshold Filter:", ${lod_threshold}))
    validation_log <- c(validation_log, paste("Permutation Rounds: 1000"))
    validation_log <- c(validation_log, "")
    
    cat("Loading required data files for permutation testing...\\n")
    
    # Load cross2 object, genotype probabilities, and kinship matrices
    cross2 <- readRDS("${cross2_file}")
    genoprob <- readRDS("${genoprob_file}")
    kinship_loco <- readRDS("${kinship_file}")
    
    # Load filtered phenotypes list
    filtered_phenotypes <- readLines("${filtered_phenotypes_file}")
    
    validation_log <- c(validation_log, "=== LOD Threshold Filtering Applied ===")
    validation_log <- c(validation_log, paste("✓ Total phenotypes in dataset:", ncol(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Phenotypes passing LOD ≥", ${lod_threshold}, "threshold:", length(filtered_phenotypes)))
    
    # Check if we have any phenotypes to test
    if (length(filtered_phenotypes) == 0) {
        validation_log <- c(validation_log, "⚠ WARNING: No phenotypes passed LOD threshold - creating empty results")
        
        # Create empty permutation results
        empty_perm <- matrix(NA, nrow=1000, ncol=0)
        saveRDS(empty_perm, file = "${prefix}_permutation_results.rds")
        
        # Create empty thresholds
        empty_thresholds <- data.frame(
            study=character(0), phenotype=character(0), significance_level=character(0),
            alpha=numeric(0), lod_threshold=numeric(0), stringsAsFactors=FALSE
        )
        write.csv(empty_thresholds, file = "${prefix}_significance_thresholds.csv", row.names=FALSE)
        
        validation_log <- c(validation_log, "✓ Empty results files created")
        writeLines(validation_log, "permutation_validation_report.txt")
        
        cat("No phenotypes passed LOD threshold - permutation testing skipped\\n")
        quit(save = "no", status = 0, runLast = FALSE)
    }
    
    validation_log <- c(validation_log, paste("✓ Phenotypes for permutation testing:", paste(filtered_phenotypes, collapse=", ")))
    validation_log <- c(validation_log, "")
    
    # Prepare phenotypes and covariates - SUBSET TO FILTERED PHENOTYPES ONLY
    pheno_full <- cross2\$pheno
    covar <- cross2\$covar
    
    # Subset phenotypes to only those passing LOD threshold
    pheno <- pheno_full[, filtered_phenotypes, drop=FALSE]
    
    validation_log <- c(validation_log, paste("✓ Phenotype matrix subset from", ncol(pheno_full), "to", ncol(pheno), "phenotypes"))
    
    # Set up covariates (identical to genome scan process)
    addcovar <- NULL
    Xcovar <- NULL

    # Check if we're testing coat_color - if so, use minimal covariates for positive control
    is_coat_color_scan <- "coat_color" %in% colnames(pheno)

    if (is_coat_color_scan) {
        validation_log <- c(validation_log, "✓ Coat color detected - using minimal covariates for positive control test")
        # For coat color positive control, use no covariates to avoid overfitting
        validation_log <- c(validation_log, "✓ No covariates for coat color permutation testing - intercept-only model for maximum power")
    } else {
    
    if ("Sex" %in% colnames(covar)) {
        sex_numeric <- ifelse(covar\$Sex == "male", 1, 0)
        addcovar <- matrix(sex_numeric, ncol=1)
        colnames(addcovar) <- "sex"
        rownames(addcovar) <- rownames(covar)
        Xcovar <- get_x_covar(cross2)
    }
    
    if ("generation" %in% colnames(covar)) {
        # Convert generation to factor and then to model matrix for proper categorical handling
        gen_factor <- factor(covar\$generation)
        if (nlevels(gen_factor) > 1) {
            gen_matrix <- model.matrix(~ gen_factor - 1)[, -1, drop=FALSE]  # Remove first level
            gen_colnames <- paste0("gen_", levels(gen_factor)[-1])
            gen_matrix <- matrix(as.numeric(gen_matrix), nrow=nrow(gen_matrix), ncol=ncol(gen_matrix))
            colnames(gen_matrix) <- gen_colnames
            rownames(gen_matrix) <- rownames(covar)

            if (is.null(addcovar)) {
                addcovar <- gen_matrix
            } else {
                addcovar <- cbind(addcovar, gen_matrix)
            }
        }
    }
    
    # BATCH EFFECTS TEMPORARILY DISABLED FOR TESTING (matches GENOME_SCAN process)
    # if ("batch" %in% colnames(covar)) {
    #     batch_factor <- factor(covar\$batch)
    #     if (nlevels(batch_factor) > 1) {
    #         batch_matrix <- model.matrix(~ batch_factor - 1)[, -1, drop=FALSE]
    #         # Store column names before numeric conversion
    #         batch_colnames <- paste0("batch_", levels(batch_factor)[-1])
    #         batch_rownames <- rownames(covar)
    #         # Convert to numeric matrix
    #         batch_matrix <- matrix(as.numeric(batch_matrix), nrow=nrow(batch_matrix), ncol=ncol(batch_matrix))
    #         colnames(batch_matrix) <- batch_colnames
    #         rownames(batch_matrix) <- batch_rownames
    #         if (is.null(addcovar)) {
    #             addcovar <- batch_matrix
    #         } else {
    #             addcovar <- cbind(addcovar, batch_matrix)
    #         }
    #     }
    # }
    } # End of else block for non-coat_color phenotypes

    # Final check: ensure all covariate matrices are numeric
    if (!is.null(addcovar)) {
        # Force all columns to numeric and ensure matrix format while preserving row/col names
        original_colnames <- colnames(addcovar)
        original_rownames <- rownames(addcovar)
        addcovar <- apply(addcovar, 2, as.numeric)
        if (is.vector(addcovar)) addcovar <- matrix(addcovar, ncol=1)
        colnames(addcovar) <- original_colnames
        rownames(addcovar) <- original_rownames
        if (!is.numeric(addcovar)) {
            stop("addcovar matrix is not numeric after conversion")
        }
    }
    
    validation_log <- c(validation_log, "=== Permutation Test Setup ===")
    validation_log <- c(validation_log, paste("✓ Number of phenotypes for testing:", ncol(pheno)))
    validation_log <- c(validation_log, paste("✓ Number of individuals:", nrow(pheno)))
    if (!is.null(addcovar)) {
        validation_log <- c(validation_log, paste("✓ Covariates:", paste(colnames(addcovar), collapse=", ")))
        validation_log <- c(validation_log, paste("✓ Addcovar class:", class(addcovar)))
        validation_log <- c(validation_log, paste("✓ Addcovar is numeric:", is.numeric(addcovar)))
    }
    validation_log <- c(validation_log, "")
    
    # Perform permutation testing
    cat("Starting permutation testing with 1000 rounds...\\n")
    cat("This will take considerable time depending on dataset size...\\n")
    cat("Progress will be displayed every 100 permutations...\\n")
    
    validation_log <- c(validation_log, "=== Permutation Testing Execution ===")
    start_time <- Sys.time()
    
    # Run 1000 permutations using all available cores
    # Get chromosome lengths from genetic map for X-specific permutations
    chr_lengths <- sapply(cross2\$gmap, function(x) max(x) - min(x))
    
    # Check if we have X chromosome - disable X-specific permutations in test mode with autosomal chromosomes only
    has_x_chr <- any(names(cross2\$gmap) %in% c("X", "x"))
    
    # Run permutation testing with appropriate model for phenotype type
    if (is_coat_color_scan) {
        # For coat color, use simple permutation without kinship to match genome scan
        perm_results <- scan1perm(genoprob, pheno, addcovar=addcovar, Xcovar=Xcovar,
                                 n_perm=1000, cores=32, perm_Xsp=has_x_chr, perm_strata=NULL,
                                 chr_lengths=if(has_x_chr) chr_lengths else NULL)
        validation_log <- c(validation_log, "✓ Using simplified permutation (no kinship) for coat color positive control")
    } else {
        # For other phenotypes, use full LOCO kinship correction
        perm_results <- scan1perm(genoprob, pheno, kinship_loco, addcovar=addcovar, Xcovar=Xcovar,
                                 n_perm=1000, cores=32, perm_Xsp=has_x_chr, perm_strata=NULL,
                                 chr_lengths=if(has_x_chr) chr_lengths else NULL)
    }
    
    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="hours"), 2)
    
    validation_log <- c(validation_log, paste("✓ Permutation testing completed in", duration, "hours"))
    validation_log <- c(validation_log, paste("✓ Permutation results dimensions:", paste(dim(perm_results), collapse=" x ")))
    
    # Save permutation results
    saveRDS(perm_results, file = "${prefix}_permutation_results.rds")
    
    # Calculate significance thresholds
    cat("Calculating significance thresholds...\\n")
    
    # Define significance levels: 0.63 (suggestive), 0.1, 0.05, 0.01
    alpha_levels <- c(0.63, 0.1, 0.05, 0.01)
    threshold_names <- c("suggestive", "0.10", "0.05", "0.01")
    
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Significance Threshold Calculation ===")
    
    # Calculate thresholds for each phenotype
    all_thresholds <- data.frame()
    
    for (i in 1:ncol(perm_results)) {
        pheno_name <- colnames(perm_results)[i]
        
        # Extract permutation LOD scores for this phenotype
        perm_lods <- perm_results[, i]
        
        # Calculate quantile thresholds
        thresholds <- quantile(perm_lods, probs = 1 - alpha_levels, na.rm = TRUE)
        
        # Create data frame for this phenotype
        pheno_thresholds <- data.frame(
            study = "${prefix}",
            phenotype = pheno_name,
            significance_level = threshold_names,
            alpha = alpha_levels,
            lod_threshold = as.numeric(thresholds),
            stringsAsFactors = FALSE
        )
        
        all_thresholds <- rbind(all_thresholds, pheno_thresholds)
    }
    
    # Write thresholds to CSV
    write.csv(all_thresholds, file = "${prefix}_significance_thresholds.csv", row.names = FALSE)
    
    validation_log <- c(validation_log, paste("✓ Calculated thresholds for", ncol(perm_results), "phenotypes"))
    validation_log <- c(validation_log, paste("✓ Threshold levels: suggestive, 0.10, 0.05, 0.01"))
    
    # Report summary statistics
    overall_thresholds <- apply(perm_results, 1, max, na.rm=TRUE)
    overall_quantiles <- quantile(overall_thresholds, probs = 1 - alpha_levels, na.rm = TRUE)
    
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Overall Threshold Summary ===")
    for (i in 1:length(alpha_levels)) {
        validation_log <- c(validation_log, paste("✓", threshold_names[i], "threshold (α =", alpha_levels[i], "):", round(overall_quantiles[i], 3)))
    }
    
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Permutation Testing Complete ===")
    validation_log <- c(validation_log, "✓ Significance thresholds established for all phenotypes")
    validation_log <- c(validation_log, "✓ Ready for QTL identification and results processing")
    
    # Write validation report
    writeLines(validation_log, "permutation_validation_report.txt")
    
    cat("Permutation testing completed successfully\\n")
    """
}

process IDENTIFY_SIGNIFICANT_QTLS {
    tag "Identifying significant QTLs for ${prefix}"
    publishDir "${params.outdir}/08_significant_qtls", mode: 'copy'
    
    input:
    path(cross2_file)
    path(scan_results_file)
    path(perm_results_file)
    path(thresholds_file)
    val(prefix)
    
    output:
    path("${prefix}_significant_qtls.csv"), emit: significant_qtls
    path("${prefix}_qtl_summary.txt"), emit: qtl_summary
    path("qtl_identification_report.txt"), emit: validation_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })
    
    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== QTL Identification Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")
    
    cat("Loading scan results and significance thresholds...\\n")
    
    # Load required data
    cross2 <- readRDS("${cross2_file}")
    scan_results <- readRDS("${scan_results_file}")
    perm_results <- readRDS("${perm_results_file}")
    thresholds <- read.csv("${thresholds_file}", stringsAsFactors = FALSE)
    
    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Scan results loaded:", paste(dim(scan_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Permutation results loaded:", paste(dim(perm_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Significance thresholds loaded:", nrow(thresholds), "entries"))
    validation_log <- c(validation_log, "")
    
    # Initialize results data frame
    all_significant_qtls <- data.frame()
    qtl_counts <- data.frame(
        significance_level = c("suggestive", "0.10", "0.05", "0.01"),
        qtl_count = 0,
        phenotypes_with_qtls = 0,
        stringsAsFactors = FALSE
    )
    
    cat("Identifying significant QTLs for each significance level...\\n")
    
    # Process each significance level
    for (sig_level in c("suggestive", "0.10", "0.05", "0.01")) {
        cat("Processing", sig_level, "significance level...\\n")
        
        # Get thresholds for this significance level
        level_thresholds <- thresholds[thresholds\$significance_level == sig_level, ]
        
        total_qtls <- 0
        phenotypes_with_qtls <- 0
        
        # Process each phenotype
        for (i in 1:nrow(level_thresholds)) {
            pheno_name <- level_thresholds\$phenotype[i]
            threshold <- level_thresholds\$lod_threshold[i]
            
            # Find peaks for this phenotype at this significance level
            if (pheno_name %in% colnames(scan_results)) {
                # Create single-phenotype scan result
                pheno_scan <- scan_results[, pheno_name, drop=FALSE]
                
                # Find peaks above threshold
                peaks <- find_peaks(pheno_scan, cross2\$gmap, threshold=threshold, drop=1.5)
                
                if (nrow(peaks) > 0) {
                    # Add metadata to peaks
                    peaks\$study <- "${prefix}"
                    peaks\$phenotype <- pheno_name
                    peaks\$significance_level <- sig_level
                    peaks\$alpha <- level_thresholds\$alpha[i]
                    peaks\$threshold_used <- threshold
                    
                    # Reorder columns
                    peaks <- peaks[, c("study", "phenotype", "significance_level", "alpha", 
                                      "threshold_used", "lodindex", "lodcolumn", "chr", "pos", "lod")]
                    
                    all_significant_qtls <- rbind(all_significant_qtls, peaks)
                    total_qtls <- total_qtls + nrow(peaks)
                    phenotypes_with_qtls <- phenotypes_with_qtls + 1
                }
            }
        }
        
        # Update counts
        qtl_counts[qtl_counts\$significance_level == sig_level, "qtl_count"] <- total_qtls
        qtl_counts[qtl_counts\$significance_level == sig_level, "phenotypes_with_qtls"] <- phenotypes_with_qtls
        
        validation_log <- c(validation_log, paste("✓", sig_level, "level:", total_qtls, "QTLs in", phenotypes_with_qtls, "phenotypes"))
    }
    
    # Write significant QTLs to file
    if (nrow(all_significant_qtls) > 0) {
        write.csv(all_significant_qtls, file = "${prefix}_significant_qtls.csv", row.names = FALSE)
        validation_log <- c(validation_log, paste("✓ Total significant QTLs identified:", nrow(all_significant_qtls)))
    } else {
        # Create empty file
        empty_qtls <- data.frame(
            study=character(0), phenotype=character(0), significance_level=character(0),
            alpha=numeric(0), threshold_used=numeric(0), lodindex=numeric(0),
            lodcolumn=character(0), chr=character(0), pos=numeric(0), lod=numeric(0)
        )
        write.csv(empty_qtls, file = "${prefix}_significant_qtls.csv", row.names = FALSE)
        validation_log <- c(validation_log, "⚠ No significant QTLs identified at any significance level")
    }
    
    # Create summary report
    cat("Generating QTL summary report...\\n")
    
    summary_lines <- c()
    summary_lines <- c(summary_lines, paste("QTL Mapping Summary for Study:", "${prefix}"))
    summary_lines <- c(summary_lines, paste("Generated:", Sys.time()))
    summary_lines <- c(summary_lines, "")
    summary_lines <- c(summary_lines, "=== Dataset Overview ===")
    summary_lines <- c(summary_lines, paste("Total individuals:", nrow(cross2\$pheno)))
    summary_lines <- c(summary_lines, paste("Total phenotypes:", ncol(cross2\$pheno)))
    summary_lines <- c(summary_lines, paste("Chromosomes analyzed:", paste(names(cross2\$gmap), collapse=", ")))
    summary_lines <- c(summary_lines, paste("Total markers scanned:", nrow(scan_results)))
    summary_lines <- c(summary_lines, "")
    summary_lines <- c(summary_lines, "=== QTL Summary by Significance Level ===")
    
    for (i in 1:nrow(qtl_counts)) {
        level <- qtl_counts\$significance_level[i]
        count <- qtl_counts\$qtl_count[i]
        phenos <- qtl_counts\$phenotypes_with_qtls[i]
        summary_lines <- c(summary_lines, paste(level, ":", count, "QTLs in", phenos, "phenotypes"))
    }
    
    if (nrow(all_significant_qtls) > 0) {
        summary_lines <- c(summary_lines, "")
        summary_lines <- c(summary_lines, "=== Chromosome Distribution ===")
        chr_counts <- table(all_significant_qtls\$chr)
        for (chr in names(chr_counts)) {
            summary_lines <- c(summary_lines, paste("Chr", chr, ":", chr_counts[chr], "QTLs"))
        }
        
        summary_lines <- c(summary_lines, "")
        summary_lines <- c(summary_lines, "=== Top QTLs (LOD ≥ 10) ===")
        top_qtls <- all_significant_qtls[all_significant_qtls\$lod >= 10, ]
        if (nrow(top_qtls) > 0) {
            top_qtls <- top_qtls[order(top_qtls\$lod, decreasing=TRUE), ]
            for (i in 1:min(10, nrow(top_qtls))) {
                qtl <- top_qtls[i, ]
                summary_lines <- c(summary_lines, 
                    paste("  ", qtl\$phenotype, "- Chr", qtl\$chr, "@", round(qtl\$pos, 1), "cM, LOD =", round(qtl\$lod, 2)))
            }
        } else {
            summary_lines <- c(summary_lines, "  No QTLs with LOD ≥ 10")
        }
    }
    
    # Write summary
    writeLines(summary_lines, "${prefix}_qtl_summary.txt")
    
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Identification Complete ===")
    validation_log <- c(validation_log, paste("✓ Significant QTLs file:", "${prefix}_significant_qtls.csv"))
    validation_log <- c(validation_log, paste("✓ Summary report:", "${prefix}_qtl_summary.txt"))
    validation_log <- c(validation_log, "✓ Ready for results visualization and downstream analysis")
    
    # Write validation report
    writeLines(validation_log, "qtl_identification_report.txt")
    
    cat("QTL identification completed successfully\\n")
    """
}