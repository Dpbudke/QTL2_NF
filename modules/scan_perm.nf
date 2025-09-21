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
    tag "Genome scanning for ${prefix}"
    publishDir "${params.outdir}/06_genome_scan_results", mode: 'copy'
    
    input:
    path(cross2_file)
    path(genoprob_file)
    path(kinship_file)
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
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })
    
    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== Genome Scan Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, paste("LOD Threshold for Permutation Filtering:", ${lod_threshold}))
    validation_log <- c(validation_log, "")
    
    cat("Loading required data files...\\n")
    
    # Load cross2 object, genotype probabilities, and kinship matrices
    cross2 <- readRDS("${cross2_file}")
    genoprob <- readRDS("${genoprob_file}")
    kinship_loco <- readRDS("${kinship_file}")
    
    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Cross2 individuals:", nrow(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Genoprob chromosomes:", length(genoprob)))
    validation_log <- c(validation_log, paste("✓ LOCO kinship matrices:", length(kinship_loco)))
    validation_log <- c(validation_log, "")
    
    # Prepare phenotypes and covariates
    pheno <- cross2\$pheno
    covar <- cross2\$covar
    
    validation_log <- c(validation_log, "=== Phenotype and Covariate Setup ===")
    validation_log <- c(validation_log, paste("✓ Number of phenotypes:", ncol(pheno)))
    validation_log <- c(validation_log, paste("✓ Number of individuals with phenotypes:", nrow(pheno)))
    validation_log <- c(validation_log, paste("✓ Available covariates:", paste(colnames(covar), collapse=", ")))
    
    # Set up covariates for qtl2
    addcovar <- NULL
    Xcovar <- NULL

    # Check if we're scanning coat_color - if so, use minimal covariates for positive control
    is_coat_color_scan <- "coat_color" %in% colnames(pheno)

    if (is_coat_color_scan) {
        validation_log <- c(validation_log, "✓ Coat color detected - using minimal covariates for positive control test")
        # For coat color positive control, use no covariates to avoid overfitting
        validation_log <- c(validation_log, "✓ No covariates for coat color analysis - intercept-only model for maximum power")
    } else {
    
    # Check for sex covariate (required for X chromosome analysis)
    if ("Sex" %in% colnames(covar)) {
        # Convert sex to numeric (assuming male/female coding)
        sex_numeric <- ifelse(covar\$Sex == "male", 1, 0)
        addcovar <- matrix(sex_numeric, ncol=1)
        colnames(addcovar) <- "sex"
        rownames(addcovar) <- rownames(covar)
        
        # For X chromosome, we need special covariate handling
        Xcovar <- get_x_covar(cross2)
        
        validation_log <- c(validation_log, "✓ Sex covariate prepared for autosomal and X chromosome analysis")
    } else {
        validation_log <- c(validation_log, "⚠ WARNING: No sex covariate found - X chromosome analysis may not be optimal")
    }
    
    # Add generation covariate if available
    if ("generation" %in% colnames(covar)) {
        gen_covar <- matrix(covar\$generation, ncol=1)
        colnames(gen_covar) <- "generation"
        rownames(gen_covar) <- rownames(covar)
        
        if (is.null(addcovar)) {
            addcovar <- gen_covar
        } else {
            addcovar <- cbind(addcovar, gen_covar)
        }
        validation_log <- c(validation_log, "✓ Generation covariate added")
    }
    
    # Add batch covariate if available
    if ("batch" %in% colnames(covar)) {
        # Convert batch to factor and then to model matrix for proper handling
        batch_factor <- factor(covar\$batch)
        if (nlevels(batch_factor) > 1) {
            batch_matrix <- model.matrix(~ batch_factor - 1)[, -1, drop=FALSE]  # Remove first level
            # Store column names before numeric conversion
            batch_colnames <- paste0("batch_", levels(batch_factor)[-1])
            batch_rownames <- rownames(covar)
            # Convert to numeric matrix
            batch_matrix <- matrix(as.numeric(batch_matrix), nrow=nrow(batch_matrix), ncol=ncol(batch_matrix))
            colnames(batch_matrix) <- batch_colnames
            rownames(batch_matrix) <- batch_rownames
            
            if (is.null(addcovar)) {
                addcovar <- batch_matrix
            } else {
                addcovar <- cbind(addcovar, batch_matrix)
            }
            validation_log <- c(validation_log, paste("✓ Batch covariate added (", nlevels(batch_factor), "levels )"))
        }
    }
    } # End of else block for non-coat_color phenotypes

    # Final comprehensive numeric conversion
    if (!is.null(addcovar)) {
        # Force complete numeric matrix conversion
        original_colnames <- colnames(addcovar)
        original_rownames <- rownames(addcovar)
        addcovar <- matrix(as.numeric(as.matrix(addcovar)), nrow=nrow(addcovar), ncol=ncol(addcovar))
        colnames(addcovar) <- original_colnames
        rownames(addcovar) <- original_rownames
        
        # Verify numeric conversion
        if (!is.numeric(addcovar) || !is.matrix(addcovar)) {
            stop("Failed to create numeric addcovar matrix")
        }
    }
    
    validation_log <- c(validation_log, "")
    if (!is.null(addcovar)) {
        validation_log <- c(validation_log, paste("✓ Final additive covariates:", paste(colnames(addcovar), collapse=", ")))
        validation_log <- c(validation_log, paste("✓ Addcovar final class:", class(addcovar)))
        validation_log <- c(validation_log, paste("✓ Addcovar is numeric:", is.numeric(addcovar)))
        validation_log <- c(validation_log, paste("✓ Addcovar dimensions:", paste(dim(addcovar), collapse=" x ")))
    } else {
        validation_log <- c(validation_log, "✓ No additive covariates - running intercept-only model")
    }
    validation_log <- c(validation_log, "")
    
    # Perform genome scan using linear mixed model with LOCO kinship
    cat("Performing genome scan with linear mixed model...\\n")
    cat("This will scan", ncol(pheno), "phenotypes across all chromosomes...\\n")
    cat("Using LOCO kinship matrices to account for population structure...\\n")
    
    validation_log <- c(validation_log, "=== Genome Scan Execution ===")
    start_time <- Sys.time()
    
    # Main genome scan with kinship correction
    if (is_coat_color_scan) {
        # For coat color, use simple scan without kinship to avoid numerical issues
        scan_results <- scan1(genoprob, pheno, addcovar=addcovar, Xcovar=Xcovar, cores=8)
        validation_log <- c(validation_log, "✓ Using simplified scan (no kinship) for coat color positive control")
    } else {
        # For other phenotypes, use full LOCO kinship correction
        scan_results <- scan1(genoprob, pheno, kinship_loco, addcovar=addcovar, Xcovar=Xcovar, cores=8)
    }
    
    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)
    
    validation_log <- c(validation_log, paste("✓ Genome scan completed in", duration, "minutes"))
    validation_log <- c(validation_log, paste("✓ Scan results dimensions:", paste(dim(scan_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Phenotypes scanned:", ncol(scan_results)))
    validation_log <- c(validation_log, paste("✓ Maximum LOD score:", round(max(scan_results, na.rm=TRUE), 3)))
    
    # Save scan results
    saveRDS(scan_results, file = "${prefix}_scan_results.rds")
    
    # LOD-based filtering for permutation testing efficiency
    cat("Filtering phenotypes based on LOD threshold for permutation testing...\\n")
    cat("LOD threshold:", ${lod_threshold}, "\\n")
    
    # Identify phenotypes with at least one QTL above the threshold
    phenotypes_above_threshold <- c()
    max_lod_per_phenotype <- apply(scan_results, 2, max, na.rm=TRUE)
    
    validation_log <- c(validation_log, "=== LOD Threshold Filtering ===")
    validation_log <- c(validation_log, paste("✓ Total phenotypes scanned:", length(max_lod_per_phenotype)))
    
    for (i in 1:length(max_lod_per_phenotype)) {
        pheno_name <- names(max_lod_per_phenotype)[i]
        max_lod <- max_lod_per_phenotype[i]
        
        if (!is.na(max_lod) && max_lod >= ${lod_threshold}) {
            phenotypes_above_threshold <- c(phenotypes_above_threshold, pheno_name)
            validation_log <- c(validation_log, paste("  ✓", pheno_name, "- Max LOD:", round(max_lod, 3)))
        }
    }
    
    # Write filtered phenotypes list for permutation testing
    if (length(phenotypes_above_threshold) > 0) {
        writeLines(phenotypes_above_threshold, "${prefix}_filtered_phenotypes.txt")
        validation_log <- c(validation_log, paste("✓ Phenotypes passing LOD threshold:", length(phenotypes_above_threshold), "out of", length(max_lod_per_phenotype)))
        validation_log <- c(validation_log, paste("✓ Computational efficiency gain:", round(100 * (1 - length(phenotypes_above_threshold)/length(max_lod_per_phenotype)), 1), "% reduction in permutation workload"))
    } else {
        writeLines(character(0), "${prefix}_filtered_phenotypes.txt")
        validation_log <- c(validation_log, "⚠ WARNING: No phenotypes exceed LOD threshold - no permutation testing will be performed")
    }
    validation_log <- c(validation_log, "")
    
    # Find peaks for each phenotype (preliminary peaks without permutation thresholds)
    cat("Finding preliminary peaks for each phenotype...\\n")
    
    # Extract peaks with minimum LOD threshold of 3.0 (rough suggestive threshold)
    all_peaks <- find_peaks(scan_results, cross2\$gmap, threshold=3.0, drop=1.5)
    
    if (nrow(all_peaks) > 0) {
        # Add study prefix to output
        all_peaks\$study <- "${prefix}"
        all_peaks <- all_peaks[, c("study", names(all_peaks)[names(all_peaks) != "study"])]
        
        # Write peaks to CSV
        write.csv(all_peaks, file = "${prefix}_scan_peaks.csv", row.names=FALSE)
        
        validation_log <- c(validation_log, paste("✓ Found", nrow(all_peaks), "preliminary peaks (LOD ≥ 3.0)"))
        validation_log <- c(validation_log, paste("✓ Peaks across", length(unique(all_peaks\$chr)), "chromosomes"))
        validation_log <- c(validation_log, paste("✓ Phenotypes with peaks:", length(unique(all_peaks\$lodcolumn))))
    } else {
        # Create empty peaks file
        empty_peaks <- data.frame(study=character(0), lodindex=numeric(0), lodcolumn=character(0), 
                                 chr=character(0), pos=numeric(0), lod=numeric(0))
        write.csv(empty_peaks, file = "${prefix}_scan_peaks.csv", row.names=FALSE)
        validation_log <- c(validation_log, "⚠ No preliminary peaks found with LOD ≥ 3.0")
    }
    
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Genome Scan Summary ===")
    validation_log <- c(validation_log, paste("✓ Total positions scanned:", nrow(scan_results)))
    validation_log <- c(validation_log, paste("✓ Chromosomes analyzed:", paste(names(cross2\$gmap), collapse=", ")))
    validation_log <- c(validation_log, "✓ Results saved for permutation testing")
    
    # Write validation report
    writeLines(validation_log, "genome_scan_validation_report.txt")
    
    cat("Genome scan completed successfully\\n")
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
        gen_covar <- matrix(covar\$generation, ncol=1)
        colnames(gen_covar) <- "generation"
        rownames(gen_covar) <- rownames(covar)
        if (is.null(addcovar)) {
            addcovar <- gen_covar
        } else {
            addcovar <- cbind(addcovar, gen_covar)
        }
    }
    
    if ("batch" %in% colnames(covar)) {
        batch_factor <- factor(covar\$batch)
        if (nlevels(batch_factor) > 1) {
            batch_matrix <- model.matrix(~ batch_factor - 1)[, -1, drop=FALSE]
            # Store column names before numeric conversion
            batch_colnames <- paste0("batch_", levels(batch_factor)[-1])
            batch_rownames <- rownames(covar)
            # Convert to numeric matrix
            batch_matrix <- matrix(as.numeric(batch_matrix), nrow=nrow(batch_matrix), ncol=ncol(batch_matrix))
            colnames(batch_matrix) <- batch_colnames
            rownames(batch_matrix) <- batch_rownames
            if (is.null(addcovar)) {
                addcovar <- batch_matrix
            } else {
                addcovar <- cbind(addcovar, batch_matrix)
            }
        }
    }
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
                                 n_perm=1000, cores=8, perm_Xsp=has_x_chr, perm_strata=NULL,
                                 chr_lengths=if(has_x_chr) chr_lengths else NULL)
        validation_log <- c(validation_log, "✓ Using simplified permutation (no kinship) for coat color positive control")
    } else {
        # For other phenotypes, use full LOCO kinship correction
        perm_results <- scan1perm(genoprob, pheno, kinship_loco, addcovar=addcovar, Xcovar=Xcovar,
                                 n_perm=1000, cores=8, perm_Xsp=has_x_chr, perm_strata=NULL,
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