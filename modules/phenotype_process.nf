process PHENOTYPE_PROCESS {
    tag "Processing phenotypes for ${prefix}"
    publishDir "${params.outdir}/01_phenotype_processing", mode: 'copy'
    
    input:
    path(phenotype_file)
    val(prefix)
    
    output:
    path("${prefix}_covar.csv"), emit: covar
    path("${prefix}_pheno.csv"), emit: pheno
    path("diagnostic_plots/*"), emit: diagnostics, optional: true
    path("validation_report.txt"), emit: validation_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(dplyr)
        library(gplots)
        library(qtl2)
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
        read.csv(phenPath, nrows=1, fill = TRUE, header = FALSE)
    }, error = function(e) {
        stop("Error reading first row of phenotype file: ", e\$message)
    })
    
    phenStart <- which(firstRow != "")
    
    if (length(phenStart) < 2) {
        stop("Could not find proper covariate/phenotype delineation in first row")
    }
    
    validation_log <- c(validation_log, paste("✓ Found covariate/phenotype boundaries at positions:", paste(phenStart, collapse=", ")))
    
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
    validation_log <- c(validation_log, paste("  - Dimensions:", nrow(mPhen), "individuals x", ncol(mPhen), "variables"))
    
    # Handle row name transformations (no changeRN.R needed)
    original_rownames <- rownames(mPhen)
    
    # Optional: Apply study prefix to sample IDs automatically
    if ("${params.auto_prefix_samples}" == "true") {
        rownames(mPhen) <- paste0("${prefix}_", rownames(mPhen))
        validation_log <- c(validation_log, "✓ Applied automatic sample ID prefixing")
    } else {
        validation_log <- c(validation_log, "✓ Kept original sample IDs (no prefix applied)")
    }
    
    # Extract covariate and phenotype data
    covar_cols <- 1:(phenStart[2]-1)
    pheno_cols <- phenStart[2]:ncol(mPhen)
    
    # Validate covariate section
    covar_data <- mPhen[, covar_cols, drop = FALSE]
    pheno_data <- mPhen[, pheno_cols, drop = FALSE]
    
    validation_log <- c(validation_log, paste("✓ Extracted covariates:", ncol(covar_data), "columns"))
    validation_log <- c(validation_log, paste("✓ Extracted phenotypes:", ncol(pheno_data), "columns"))
    
    # Validate required covariates for r/qtl2
    required_covars <- c("sex", "ngen")
    # Also check for alternative column names
    covar_names_lower <- tolower(colnames(covar_data))
    
    # Check for sex column
    has_sex <- any(c("sex") %in% covar_names_lower)
    
    # Check for generation column (ngen or generation)
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
    
    cat("Phenotype processing completed successfully\\n")
    cat("Check validation_report.txt for detailed information\\n")
    """
}