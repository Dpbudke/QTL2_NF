process PHENOTYPE_PROCESS {
    tag "Processing phenotypes for ${prefix}"
    publishDir "${params.outdir}/01_phenotype_processing", mode: 'copy'
    
    input:
    path(phenotype_file)
    val(prefix)
    val(sample_filter_json)
    
    output:
    path("${prefix}_covar.csv"), emit: covar
    path("${prefix}_pheno.csv"), emit: pheno
    path("${prefix}_valid_samples.txt"), emit: valid_samples
    path("diagnostic_plots/*"), emit: diagnostics, optional: true
    path("validation_report.txt"), emit: validation_report
    path("${prefix}_sample_filter_report.txt"), emit: filter_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(dplyr)
        library(gplots)
        library(qtl2)
        library(jsonlite)  # For parsing sample filter JSON
    })
    
    # Source custom functions (now from /opt/bin in container)
    source("/opt/bin/robustZmat.R") 
    source("/opt/bin/covCheck.R")
    
    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== Phenotype Processing Validation Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, paste("Input File:", "${phenotype_file}"))
    validation_log <- c(validation_log, "")
    
    # Read the phenotype file using params for na.strings
    phenPath <- "${phenotype_file}"
    
    # Validate file exists and is readable
    if (!file.exists(phenPath)) {
        stop("Phenotype file not found: ", phenPath)
    }
    
    validation_log <- c(validation_log, "✓ Input file exists and is readable")
    
    # Find the delineation between covariates and phenotypes
    # This handles the DO_pipe format with labels in first row
    firstRow <- tryCatch({
        read.csv(phenPath, nrows=1, fill = TRUE, header = FALSE, stringsAsFactors = FALSE)
    }, error = function(e) {
        stop("Error reading first row of phenotype file: ", e\$message)
    })
    
    # Find where "phenotype" label appears (case insensitive)
    pheno_position <- which(tolower(as.character(firstRow[1,])) == "phenotype")
    
    if (length(pheno_position) == 0) {
        stop("Could not find 'phenotype' label in first row - required for DO_pipe format")
    }
    
    # Since we skip row 1 when reading, the phenotype column position is pheno_position - 1
    # (because the first column becomes row names)
    pheno_start_col <- pheno_position - 1
    
    validation_log <- c(validation_log, paste("✓ Found 'phenotype' label at position", pheno_position))
    validation_log <- c(validation_log, paste("✓ Phenotypes will start at column", pheno_start_col, "after reading data"))
    
    # Pull in phenotype file
    mPhen <- tryCatch({
        read.csv(phenPath, 
                 header = TRUE,
                 row.names = 1,
                 na.strings = c("na", "NA", "N/A", ""),
                 skip = 1)
    }, error = function(e) {
        stop("Error reading phenotype data: ", e\$message)
    })
    
    validation_log <- c(validation_log, paste("✓ Successfully read phenotype file"))
    validation_log <- c(validation_log, paste("  - Initial dimensions:", nrow(mPhen), "individuals x", ncol(mPhen), "variables"))
    
    # =============================================================================
    # SAMPLE FILTERING BY COVARIATES
    # =============================================================================
    
    filter_log <- c()
    filter_log <- c(filter_log, paste("=== Sample Filtering Report ==="))
    filter_log <- c(filter_log, paste("Timestamp:", Sys.time()))
    filter_log <- c(filter_log, paste("Study Prefix:", "${prefix}"))
    filter_log <- c(filter_log, "")
    
    # Parse sample filter JSON
    sample_filter_str <- '${sample_filter_json}'

    # In test mode, ignore sample filtering to use all samples for maximum power
    if ("${params.test_mode}" == "true") {
        cat("TEST MODE: Ignoring sample filtering to use all available samples for positive control...\\n")
        filter_log <- c(filter_log, "=== TEST MODE: Sample Filtering Disabled ===")
        filter_log <- c(filter_log, "✓ Using ALL samples for coat color positive control test")
        filter_log <- c(filter_log, paste("✓ Total samples:", nrow(mPhen)))
        validation_log <- c(validation_log, "✓ TEST MODE: Sample filtering disabled - using all samples")
    } else if (sample_filter_str != "null" && nchar(sample_filter_str) > 0) {
        cat("Applying sample filtering based on covariates...\\n")
        
        # Parse JSON filter
        filter_criteria <- tryCatch({
            jsonlite::fromJSON(sample_filter_str)
        }, error = function(e) {
            stop("Error parsing sample_filter JSON: ", e\$message, "\\nProvided: ", sample_filter_str)
        })
        
        filter_log <- c(filter_log, "=== Filter Criteria Applied ===")
        for (covar_name in names(filter_criteria)) {
            allowed_values <- filter_criteria[[covar_name]]
            filter_log <- c(filter_log, paste("✓", covar_name, "in", paste(allowed_values, collapse=", ")))
        }
        filter_log <- c(filter_log, "")
        
        # Apply filters progressively
        original_n <- nrow(mPhen)
        filtered_samples <- rownames(mPhen)
        
        for (covar_name in names(filter_criteria)) {
            if (covar_name %in% colnames(mPhen)) {
                allowed_values <- filter_criteria[[covar_name]]
                
                # Get current covariate values for remaining samples
                current_values <- mPhen[filtered_samples, covar_name]
                
                # Filter samples
                keep_samples <- filtered_samples[current_values %in% allowed_values]
                n_before <- length(filtered_samples)
                n_after <- length(keep_samples)
                
                filter_log <- c(filter_log, paste("Filter by", covar_name, ":"))
                filter_log <- c(filter_log, paste("  Before:", n_before, "samples"))
                filter_log <- c(filter_log, paste("  After:", n_after, "samples"))
                filter_log <- c(filter_log, paste("  Removed:", n_before - n_after, "samples"))
                
                filtered_samples <- keep_samples
            } else {
                warning("Covariate '", covar_name, "' not found in data. Available: ", paste(colnames(mPhen), collapse=", "))
                filter_log <- c(filter_log, paste("⚠ WARNING: Covariate", covar_name, "not found in data"))
            }
        }
        
        # Apply final filtering to data
        if (length(filtered_samples) == 0) {
            stop("No samples remain after filtering! Please check filter criteria.")
        }
        
        mPhen <- mPhen[filtered_samples, ]
        
        filter_log <- c(filter_log, "")
        filter_log <- c(filter_log, "=== Final Filter Results ===")
        filter_log <- c(filter_log, paste("✓ Original samples:", original_n))
        filter_log <- c(filter_log, paste("✓ Filtered samples:", nrow(mPhen)))
        filter_log <- c(filter_log, paste("✓ Reduction:", round(100 * (1 - nrow(mPhen)/original_n), 1), "%"))
        
        # Report final covariate distributions
        for (covar_name in names(filter_criteria)) {
            if (covar_name %in% colnames(mPhen)) {
                covar_table <- table(mPhen[, covar_name])
                filter_log <- c(filter_log, paste("✓", covar_name, "distribution:", paste(names(covar_table), "=", covar_table, collapse=", ")))
            }
        }
        
        validation_log <- c(validation_log, paste("✓ Sample filtering applied - reduced to", nrow(mPhen), "samples"))
        
    } else {
        filter_log <- c(filter_log, "=== No Filtering Applied ===")
        filter_log <- c(filter_log, "✓ Using all available samples (sample_filter = null)")
        filter_log <- c(filter_log, paste("✓ Total samples:", nrow(mPhen)))
        validation_log <- c(validation_log, "✓ No sample filtering requested - using all samples")
    }
    
    validation_log <- c(validation_log, paste("  - Final dimensions:", nrow(mPhen), "individuals x", ncol(mPhen), "variables"))
    
    # Handle row name transformations (no changeRN.R needed)
    original_rownames <- rownames(mPhen)
    
    # Universal sample ID standardization for cross-project compatibility
    # Genotype files (FinalReport) use numeric IDs, so standardize all IDs to numeric format
    original_ids <- rownames(mPhen)
    
    # Detect and standardize sample ID format
    # Check if sample IDs follow pattern: ProjectName_Number (e.g., AFRI_2, DOchln_1, DOF1_15, etc.)
    # Look for underscore followed by numbers
    has_project_prefix <- any(grepl("_[0-9]+", original_ids))
    
    if (has_project_prefix) {
        # Extract numeric suffix from project-prefixed sample IDs  
        numeric_ids <- sub(".*_([0-9]+)", "\\\\1", original_ids)
        rownames(mPhen) <- numeric_ids
        
        # Get example for validation log
        example_original <- original_ids[1]
        example_new <- numeric_ids[1]
        
        validation_log <- c(validation_log, "✓ Standardized sample IDs to numeric format for universal compatibility")
        validation_log <- c(validation_log, paste("  Example:", example_original, "->", example_new))
        validation_log <- c(validation_log, paste("  Total IDs converted:", length(original_ids)))
    } else if ("${params.auto_prefix_samples}" == "true") {
        # Only apply prefix if explicitly requested and IDs are already numeric
        rownames(mPhen) <- paste0("${prefix}_", rownames(mPhen))
        validation_log <- c(validation_log, "✓ Applied automatic sample ID prefixing")
    } else {
        validation_log <- c(validation_log, "✓ Sample IDs already in compatible format (no changes needed)")
    }
    
    # Extract covariate and phenotype data based on phenotype label position
    covar_cols <- 1:(pheno_start_col - 1)
    pheno_cols <- pheno_start_col:ncol(mPhen)
    
    # Validate covariate section
    covar_data <- mPhen[, covar_cols, drop = FALSE]
    pheno_data <- mPhen[, pheno_cols, drop = FALSE]
    
    validation_log <- c(validation_log, paste("✓ Extracted covariates:", ncol(covar_data), "columns"))
    validation_log <- c(validation_log, paste("  Covariate names:", paste(colnames(covar_data), collapse = ", ")))
    validation_log <- c(validation_log, paste("✓ Extracted phenotypes:", ncol(pheno_data), "columns"))

    # =============================================================================
    # TEST MODE: COAT COLOR POSITIVE CONTROL
    # =============================================================================

    if ("${params.test_mode}" == "true") {
        cat("TEST MODE: Converting coat_color to phenotype for positive control test...\\n")

        if ("coat_color" %in% colnames(covar_data)) {
            # Convert coat_color from covariate to numeric phenotype
            coat_color_values <- covar_data\$coat_color

            # Create numeric phenotype coding for all coat colors:
            # agouti = 3 (wild-type, highest value - dominant/complex genetics)
            # black = 2 (common recessive)
            # white = 1 (albino-related)
            # tan = 0 (rare variant)
            coat_color_numeric <- ifelse(coat_color_values == "agouti", 3,
                                       ifelse(coat_color_values == "black", 2,
                                            ifelse(coat_color_values == "white", 1,
                                                 ifelse(coat_color_values == "tan", 0, NA))))

            # Replace all phenotypes with just coat_color
            pheno_data <- data.frame(coat_color = coat_color_numeric, row.names = rownames(pheno_data))

            # Remove coat_color from covariates
            covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]

            validation_log <- c(validation_log, "")
            validation_log <- c(validation_log, "=== TEST MODE: Coat Color Positive Control ===")
            validation_log <- c(validation_log, "✓ Converted coat_color from covariate to phenotype")
            validation_log <- c(validation_log, paste("✓ Agouti samples (coded as 3):", sum(coat_color_numeric == 3, na.rm = TRUE)))
            validation_log <- c(validation_log, paste("✓ Black samples (coded as 2):", sum(coat_color_numeric == 2, na.rm = TRUE)))
            validation_log <- c(validation_log, paste("✓ White samples (coded as 1):", sum(coat_color_numeric == 1, na.rm = TRUE)))
            validation_log <- c(validation_log, paste("✓ Tan samples (coded as 0):", sum(coat_color_numeric == 0, na.rm = TRUE)))
            validation_log <- c(validation_log, paste("✓ Missing coat_color:", sum(is.na(coat_color_numeric))))
            validation_log <- c(validation_log, "✓ Phenotype matrix now contains ONLY coat_color for chr 2 positive control")
            validation_log <- c(validation_log, paste("✓ Updated covariates:", ncol(covar_data), "columns"))
            validation_log <- c(validation_log, paste("✓ Updated phenotypes:", ncol(pheno_data), "columns"))

        } else {
            stop("TEST MODE ERROR: coat_color not found in covariates. Available covariates: ",
                 paste(colnames(covar_data), collapse = ", "))
        }
    }
    # ALWAYS remove coat_color from covariates (it's a genetic phenotype, not a covariate)
    if ("coat_color" %in% colnames(covar_data)) {
        covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
        validation_log <- c(validation_log, "")
        validation_log <- c(validation_log, "=== Covariate Filtering ===")
        validation_log <- c(validation_log, "✓ Removed coat_color from covariates (genetic phenotype, not environmental)")
        validation_log <- c(validation_log, paste("✓ Remaining covariates:", paste(colnames(covar_data), collapse = ", ")))
    }

    validation_log <- c(validation_log, paste("  First 5 phenotype names:", paste(head(colnames(pheno_data), 5), collapse = ", ")))

    # Validate required covariates for r/qtl2
    required_covars <- c("sex", "ngen")
    # Also check for alternative column names
    covar_names_lower <- tolower(colnames(covar_data))
    
    # Check for sex column (case insensitive)
    has_sex <- any(c("sex") %in% covar_names_lower)
    
    # Check for generation column (ngen or generation, case insensitive)
    has_ngen <- any(c("ngen", "generation") %in% covar_names_lower)
    
    missing_covars <- c()
    if (!has_sex) missing_covars <- c(missing_covars, "sex")
    if (!has_ngen) missing_covars <- c(missing_covars, "ngen/generation")
    
    if (length(missing_covars) > 0) {
        validation_log <- c(validation_log, paste("⚠ WARNING: Missing required covariates:", paste(missing_covars, collapse=", ")))
        validation_log <- c(validation_log, "  These are required for proper r/qtl2 cross object creation")
        validation_log <- c(validation_log, paste("  Available covariates:", paste(colnames(covar_data), collapse=", ")))
    } else {
        validation_log <- c(validation_log, "✓ All required covariates present (sex, ngen/generation)")
    }
    
    # Check phenotype data is strictly numeric
    non_numeric_phenos <- !sapply(pheno_data, is.numeric)
    if (any(non_numeric_phenos)) {
        problematic_cols <- colnames(pheno_data)[non_numeric_phenos]
        stop("Non-numeric data found in phenotype columns: ", paste(problematic_cols, collapse=", "),
             "\\nPhenotype matrix must be strictly numeric for r/qtl2 compatibility")
    }
    
    validation_log <- c(validation_log, "✓ All phenotype columns are numeric")
    
    # Check for sample ID consistency
    sample_ids <- rownames(mPhen)
    duplicate_ids <- sample_ids[duplicated(sample_ids)]
    if (length(duplicate_ids) > 0) {
        stop("Duplicate sample IDs found: ", paste(duplicate_ids, collapse=", "))
    }
    
    validation_log <- c(validation_log, "✓ No duplicate sample IDs")
    validation_log <- c(validation_log, paste("✓ Sample ID format example:", head(sample_ids, 3)))
    
    # Write out the covar file
    write.csv(covar_data, 
              file = paste0("${prefix}_covar.csv"), 
              row.names = TRUE)
    
    # Write out the pheno file  
    write.csv(pheno_data, 
              file = paste0("${prefix}_pheno.csv"), 
              row.names = TRUE)
    
    validation_log <- c(validation_log, "✓ Successfully wrote covariate file")
    validation_log <- c(validation_log, "✓ Successfully wrote phenotype file")
    
    # Generate diagnostic plots if data is suitable
    tryCatch({
        # Create diagnostic directory
        dir.create("diagnostic_plots", recursive = TRUE)
        
        # Check if we have enough data for diagnostics
        if (ncol(pheno_data) < 2 || nrow(pheno_data) < 10) {
            validation_log <<- c(validation_log, "⚠ WARNING: Insufficient data for diagnostic plots (need ≥2 phenotypes, ≥10 samples)")
        } else if (ncol(pheno_data) > 100) {
            validation_log <<- c(validation_log, "⚠ INFO: Skipping diagnostic plots for large phenotype dataset (>100 phenotypes)")
            validation_log <<- c(validation_log, "  This is typical for eQTL datasets to avoid memory/time constraints")
            validation_log <<- c(validation_log, "  Data validation and file processing completed successfully")
        } else {
            
            # Plot 1 & 2 - Robust Z-score matrices
            tryCatch({
                # robustZmat function already loaded from /opt/bin
                
                # Full matrix plot
                robustZmat(pheno_data,
                           prefix = "${prefix}",
                           path = './diagnostic_plots/',
                           rowFont = 1,
                           colFont = 10,
                           margins = c(150,45),
                           pdfWid = 150,
                           pdfHei = 150)
                
                # Zoomed outlier plot
                robustZmat(pheno_data,
                           prefix = "${prefix}",
                           path = './diagnostic_plots/',
                           rowFont = 10,
                           colFont = 20,
                           zoom = TRUE,
                           valSize = 20,
                           margins = c(100,60))
                
                validation_log <<- c(validation_log, "✓ Generated robust Z-score diagnostic plots")
                
            }, error = function(e) {
                if (grepl("no outliers", e\$message)) {
                    validation_log <<- c(validation_log, "✓ No outliers detected in data (robust Z-score plots skipped)")
                } else {
                    validation_log <<- c(validation_log, paste("⚠ WARNING: Robust Z-score plots failed:", e\$message))
                }
            })
            
            # Plot 3 - Batch effect diagnostic  
            tryCatch({
                # covCheck function already loaded from /opt/bin
                
                covCheck(covar_data, pheno_data,
                         prefix = "${prefix}",
                         path = './diagnostic_plots/',
                         pdfWid = 100,
                         pdfHei = 100,
                         margins = c(100,80),
                         rowFont = 15, 
                         colFont = 10)
                
                validation_log <<- c(validation_log, "✓ Generated batch effect diagnostic plots")
                
            }, error = function(e) {
                validation_log <<- c(validation_log, paste("⚠ WARNING: Batch effect plots failed:", e\$message))
            })
        }
        
    }, error = function(e) {
        validation_log <<- c(validation_log, paste("⚠ WARNING: Diagnostic plot generation failed:", e\$message))
    })
    
    # Final validation summary
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Summary ===")
    validation_log <- c(validation_log, paste("Total individuals:", nrow(mPhen)))
    validation_log <- c(validation_log, paste("Total covariates:", ncol(covar_data)))
    validation_log <- c(validation_log, paste("Total phenotypes:", ncol(pheno_data)))
    validation_log <- c(validation_log, "Files ready for r/qtl2 processing")
    
    # Write validation report
    writeLines(validation_log, "validation_report.txt")
    
    # Write sample filter report
    writeLines(filter_log, "${prefix}_sample_filter_report.txt")

    # Write valid sample list for Module 2 to use for genotype filtering
    valid_sample_ids <- rownames(mPhen)
    writeLines(valid_sample_ids, "${prefix}_valid_samples.txt")
    validation_log <- c(validation_log, paste("✓ Wrote valid sample list:", length(valid_sample_ids), "samples"))

    cat("Phenotype processing completed successfully\\n")
    cat("Check validation_report.txt for detailed information\\n")
    """
}