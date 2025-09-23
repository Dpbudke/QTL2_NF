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