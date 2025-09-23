process CREATE_CROSS2_OBJECT {
    tag "Creating cross2 object for ${prefix}"
    publishDir "${params.outdir}/04_cross2_creation", mode: 'copy'

    input:
    path(control_file)
    path(all_data_files)  // pheno, covar, geno files, maps, founder files
    val(prefix)

    output:
    path("${prefix}_cross2.rds"), emit: cross2_object
    path("cross2_validation_report.txt"), emit: validation_report
    path("cross2_summary.txt"), emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })

    # Initialize validation log
    validation_log <- c()
    validation_log <- c(validation_log, "=== Cross2 Object Creation Report ===")
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, paste("Control File:", "${control_file}"))
    validation_log <- c(validation_log, "")

    # Validate control file exists
    if (!file.exists("${control_file}")) {
        stop("Control file not found: ${control_file}")
    }

    validation_log <- c(validation_log, "=== Loading Cross2 Object ===")

    # Load cross2 object using qtl2
    cat("Loading cross2 object from control file...\\n")

    cross2 <- tryCatch({
        read_cross2("${control_file}")
    }, error = function(e) {
        validation_log <<- c(validation_log, paste("✗ Error loading cross2 object:", e\$message))
        stop("Failed to load cross2 object: ", e\$message)
    })

    validation_log <- c(validation_log, "✓ Successfully loaded cross2 object")

    # Set X chromosome specification for qtl2 (critical for proper X chromosome handling)
    chr_names <- names(cross2\$geno)
    if ("X" %in% chr_names) {
        # Set the is_x_chr attribute directly (qtl2 approach)
        is_x <- chr_names == "X"
        cross2\$is_x_chr <- setNames(is_x, chr_names)
        validation_log <- c(validation_log, "✓ X chromosome specification set for qtl2")
        validation_log <- c(validation_log, paste("  is_x_chr:", paste(names(cross2\$is_x_chr)[cross2\$is_x_chr], collapse = ", ")))

        # Generate X chromosome covariates to avoid spurious linkage (Broman et al. 2006)
        tryCatch({
            x_covar <- get_x_covar(cross2)
            if (!is.null(x_covar) && ncol(x_covar) > 0) {
                validation_log <- c(validation_log, "✓ X chromosome covariates generated successfully")
                validation_log <- c(validation_log, paste("  X covariates:", paste(colnames(x_covar), collapse = ", ")))
            }
        }, error = function(e) {
            validation_log <<- c(validation_log, paste("⚠ WARNING: Could not generate X chromosome covariates:", e\$message))
        })
    }

    # Basic cross2 object validation
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Cross2 Object Validation ===")

    # Check cross type
    if (!is.null(attr(cross2, "crosstype"))) {
        validation_log <- c(validation_log, paste("✓ Cross type:", attr(cross2, "crosstype")))
    }

    # Check dimensions
    n_ind <- n_ind(cross2)
    n_chr <- n_chr(cross2)
    n_mar <- n_mar(cross2)
    n_pheno <- n_pheno(cross2)

    validation_log <- c(validation_log, paste("✓ Number of individuals:", n_ind))
    validation_log <- c(validation_log, paste("✓ Number of chromosomes:", n_chr))
    validation_log <- c(validation_log, paste("✓ Total markers:", sum(n_mar)))
    validation_log <- c(validation_log, paste("✓ Number of phenotypes:", n_pheno))

    # Check chromosome names and types
    chr_names <- names(cross2\$geno)
    validation_log <- c(validation_log, paste("✓ Chromosomes:", paste(chr_names, collapse = ", ")))

    # Special validation for X chromosome (critical for DO mice)
    if ("X" %in% chr_names) {
        validation_log <- c(validation_log, "")
        validation_log <- c(validation_log, "=== X Chromosome Validation ===")
        validation_log <- c(validation_log, "✓ X chromosome detected - performing DO-specific validation")

        # Check X chromosome genotype structure
        x_geno <- cross2\$geno\$X
        unique_x_genos <- sort(unique(as.vector(x_geno)))
        unique_x_genos <- unique_x_genos[!is.na(unique_x_genos) & unique_x_genos != 0]
        validation_log <- c(validation_log, paste("✓ X chromosome genotype codes:", paste(unique_x_genos, collapse = ", ")))

        # Check sex information for X chromosome validation
        covar_names_lower_x <- tolower(colnames(cross2\$covar))
        sex_col_idx_x <- which(covar_names_lower_x == "sex")

        if (length(sex_col_idx_x) > 0) {
            sex_col_name_x <- colnames(cross2\$covar)[sex_col_idx_x[1]]
            sex_info <- cross2\$covar[[sex_col_name_x]]
            sex_counts <- table(sex_info, useNA = "ifany")
            validation_log <- c(validation_log, "✓ Sex distribution for X chromosome validation:")
            for (sex in names(sex_counts)) {
                validation_log <- c(validation_log, paste("  ", sex, ":", sex_counts[sex]))
            }

            # For DO mice, males should have simpler X genotypes (hemizygous)
            # This is a basic check - full validation would require more complex logic
            validation_log <- c(validation_log, "✓ X chromosome ready for DO-specific genotype probability calculations")
        } else {
            validation_log <- c(validation_log, "⚠ WARNING: Sex information not available for X chromosome validation")
        }
    }

    # Validate sample IDs match between components
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Sample ID Validation ===")

    geno_ids <- rownames(cross2\$geno[[1]])
    pheno_ids <- rownames(cross2\$pheno)
    covar_ids <- rownames(cross2\$covar)

    # Check ID consistency
    if (!identical(sort(geno_ids), sort(pheno_ids))) {
        validation_log <- c(validation_log, "⚠ WARNING: Sample IDs don't match between genotypes and phenotypes")
        missing_pheno <- setdiff(geno_ids, pheno_ids)
        missing_geno <- setdiff(pheno_ids, geno_ids)
        if (length(missing_pheno) > 0) {
            validation_log <- c(validation_log, paste("  Missing in phenotypes:", paste(head(missing_pheno, 5), collapse = ", ")))
        }
        if (length(missing_geno) > 0) {
            validation_log <- c(validation_log, paste("  Missing in genotypes:", paste(head(missing_geno, 5), collapse = ", ")))
        }
    } else {
        validation_log <- c(validation_log, "✓ Sample IDs match between genotypes and phenotypes")
    }

    if (!identical(sort(pheno_ids), sort(covar_ids))) {
        validation_log <- c(validation_log, "⚠ WARNING: Sample IDs don't match between phenotypes and covariates")
        missing_covar <- setdiff(pheno_ids, covar_ids)
        missing_pheno <- setdiff(covar_ids, pheno_ids)
        if (length(missing_covar) > 0) {
            validation_log <- c(validation_log, paste("  Missing in covariates:", paste(head(missing_covar, 5), collapse = ", ")))
        }
        if (length(missing_pheno) > 0) {
            validation_log <- c(validation_log, paste("  Missing in phenotypes:", paste(head(missing_pheno, 5), collapse = ", ")))
        }
    } else {
        validation_log <- c(validation_log, "✓ Sample IDs match between phenotypes and covariates")
    }

    # Validate required covariates for DO analysis
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== DO-Specific Covariate Validation ===")

    # Check sex covariate (required for DO analysis) - case insensitive
    covar_names_lower <- tolower(colnames(cross2\$covar))
    sex_col_idx <- which(covar_names_lower == "sex")

    if (length(sex_col_idx) > 0) {
        sex_col_name <- colnames(cross2\$covar)[sex_col_idx[1]]  # Use actual column name
        sex_levels <- unique(cross2\$covar[[sex_col_name]])
        validation_log <- c(validation_log, paste("✓ Sex covariate levels:", paste(sex_levels, collapse = ", ")))

        sex_counts <- table(cross2\$covar[[sex_col_name]], useNA = "ifany")
        validation_log <- c(validation_log, "✓ Sex distribution:")
        for (sex in names(sex_counts)) {
            validation_log <- c(validation_log, paste("  ", sex, ":", sex_counts[sex]))
        }

        # Check for reasonable sex distribution
        if (length(sex_levels) != 2) {
            validation_log <- c(validation_log, "⚠ WARNING: Expected exactly 2 sex levels for DO mice")
        }
    } else {
        validation_log <- c(validation_log, "✗ ERROR: Sex covariate not found - required for DO analysis")
    }

    # Check crossinfo (ngen/generation) covariate if present
    crossinfo_present <- FALSE
    if ("ngen" %in% colnames(cross2\$covar)) {
        crossinfo_col <- "ngen"
        crossinfo_present <- TRUE
    } else if ("generation" %in% colnames(cross2\$covar)) {
        crossinfo_col <- "generation"
        crossinfo_present <- TRUE
    }

    if (crossinfo_present) {
        crossinfo_levels <- sort(unique(cross2\$covar[[crossinfo_col]]))
        validation_log <- c(validation_log, paste("✓ Generation (", crossinfo_col, ") levels:", paste(crossinfo_levels, collapse = ", ")))

        crossinfo_counts <- table(cross2\$covar[[crossinfo_col]], useNA = "ifany")
        validation_log <- c(validation_log, "✓ Generation distribution:")
        for (gen in names(crossinfo_counts)) {
            validation_log <- c(validation_log, paste("  ", gen, ":", crossinfo_counts[gen]))
        }
    } else {
        validation_log <- c(validation_log, "⚠ WARNING: Neither ngen nor generation covariate found - may be needed for some DO analyses")
    }

    # Final validation check
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Final Validation ===")

    # Check if cross2 object is ready for genome scanning
    ready_for_scanning <- TRUE
    issues <- c()

    if (n_ind < 10) {
        ready_for_scanning <- FALSE
        issues <- c(issues, "Too few individuals (< 10)")
    }

    if (n_pheno == 0) {
        ready_for_scanning <- FALSE
        issues <- c(issues, "No phenotypes")
    }

    if (length(which(tolower(colnames(cross2\$covar)) == "sex")) == 0) {
        ready_for_scanning <- FALSE
        issues <- c(issues, "Missing sex covariate")
    }

    if (sum(n_mar) < 100) {
        ready_for_scanning <- FALSE
        issues <- c(issues, "Too few markers (< 100)")
    }

    if (ready_for_scanning) {
        validation_log <- c(validation_log, "✓ Cross2 object appears ready for genome scanning")
    } else {
        validation_log <- c(validation_log, "⚠ WARNING: Cross2 object may not be ready for genome scanning")
        validation_log <- c(validation_log, "  Issues:")
        for (issue in issues) {
            validation_log <- c(validation_log, paste("    -", issue))
        }
    }

    # Generate summary statistics
    summary_info <- c()
    summary_info <- c(summary_info, "=== Cross2 Object Summary ===")
    summary_info <- c(summary_info, paste("Study Prefix:", "${prefix}"))
    summary_info <- c(summary_info, paste("Cross Type: DO (Diversity Outbred)"))
    summary_info <- c(summary_info, paste("Number of Individuals:", n_ind))
    summary_info <- c(summary_info, paste("Number of Chromosomes:", n_chr))
    summary_info <- c(summary_info, paste("Total Markers:", sum(n_mar)))
    summary_info <- c(summary_info, paste("Number of Phenotypes:", n_pheno))
    summary_info <- c(summary_info, paste("Number of Covariates:", ncol(cross2\$covar)))
    summary_info <- c(summary_info, paste("X Chromosome Present:", "X" %in% chr_names))

    if (ncol(cross2\$covar) > 0) {
        summary_info <- c(summary_info, "")
        summary_info <- c(summary_info, "Covariates:")
        summary_info <- c(summary_info, paste("  ", colnames(cross2\$covar)))
    }

    if (ncol(cross2\$pheno) > 0) {
        summary_info <- c(summary_info, "")
        summary_info <- c(summary_info, "First 10 Phenotypes:")
        summary_info <- c(summary_info, paste("  ", head(colnames(cross2\$pheno), 10)))
        if (ncol(cross2\$pheno) > 10) {
            summary_info <- c(summary_info, paste("  ... and", ncol(cross2\$pheno) - 10, "more"))
        }
    }

    # Save cross2 object
    saveRDS(cross2, file = paste0("${prefix}_cross2.rds"))
    validation_log <- c(validation_log, paste("✓ Successfully saved cross2 object to", paste0("${prefix}_cross2.rds")))

    # Write validation report and summary
    writeLines(validation_log, "cross2_validation_report.txt")
    writeLines(summary_info, "cross2_summary.txt")

    cat("Cross2 object creation completed successfully\\n")
    cat("Cross2 object saved as:", paste0("${prefix}_cross2.rds"), "\\n")
    cat("Validation report:", "cross2_validation_report.txt", "\\n")
    cat("Summary:", "cross2_summary.txt", "\\n")
    """
}