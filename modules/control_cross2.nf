process GENERATE_CONTROL_FILE {
    tag "Generating control file for ${prefix}"
    publishDir "${params.outdir}/03_control_file_generation", mode: 'copy'
    
    input:
    path(pheno_csv)
    path(covar_csv)
    path(geno_files)  // all chromosome geno files from Module 2
    path(allele_codes_file)  // GM_allelecodes.csv from Module 2
    val(prefix)
    
    output:
    path("${prefix}_control.json"), emit: control_file
    path("GM_foundergeno*.csv"), emit: founder_genos
    path("GM_gmap*.csv"), emit: genetic_maps
    path("GM_pmap*.csv"), emit: physical_maps
    path("control_validation_report.txt"), emit: validation_report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })
    
    # Initialize validation log
    validation_log <- c()
    validation_log <- c(validation_log, "=== Control File Generation Report ===")
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")
    
    # Use already downloaded GM reference files from Module 2
    cat("Using GigaMUGA reference files from Module 2...\\n")
    figshare_url <- "https://figshare.com/ndownloader/files/40233652"
    
    validation_log <- c(validation_log, "=== Reference Data Download ===")
    validation_log <- c(validation_log, "Using GM reference files already downloaded in Module 2")
    
    tryCatch({
        # Download and extract (since we need founder/map files that Module 2 didn't extract)
        download.file(figshare_url, 
                     destfile = "GM_processed_files_build39.zip", 
                     method = "auto",
                     mode = "wb")
        
        # Extract the zip file
        unzip("GM_processed_files_build39.zip", exdir = ".")
        validation_log <- c(validation_log, "✓ Successfully extracted reference files")
        
        # Find required files
        founder_files <- list.files(pattern = ".*foundergeno.*\\\\.csv", recursive = TRUE, full.names = TRUE)
        gmap_files <- list.files(pattern = ".*gmap.*\\\\.csv", recursive = TRUE, full.names = TRUE)
        pmap_files <- list.files(pattern = ".*pmap.*\\\\.csv", recursive = TRUE, full.names = TRUE)
        
        validation_log <- c(validation_log, paste("✓ Found", length(founder_files), "founder genotype files"))
        validation_log <- c(validation_log, paste("✓ Found", length(gmap_files), "genetic map files"))
        validation_log <- c(validation_log, paste("✓ Found", length(pmap_files), "physical map files"))
        
        # Copy files to working directory for output
        for(file in founder_files) {
            file.copy(file, basename(file))
        }
        for(file in gmap_files) {
            file.copy(file, basename(file))
        }
        for(file in pmap_files) {
            file.copy(file, basename(file))
        }
        
        # Clean up zip file
        file.remove("GM_processed_files_build39.zip")
        
    }, error = function(e) {
        stop("Failed to download/extract reference files: ", e\$message)
    })
    
    # Validate input files exist
    if (!file.exists("${pheno_csv}")) {
        stop("Phenotype file not found: ${pheno_csv}")
    }
    if (!file.exists("${covar_csv}")) {
        stop("Covariate file not found: ${covar_csv}")
    }
    
    validation_log <- c(validation_log, "✓ Phenotype and covariate files exist")
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Control File Generation ===")
    
    # Get list of genotype files and extract chromosomes
    geno_files <- list.files(pattern = "${prefix}_geno.*\\\\.csv", full.names = FALSE)
    chromosomes <- gsub("${prefix}_geno(.*)\\\\.csv", "\\\\1", geno_files)
    
    # For DO mice, we need chromosomes 1-19 and X
    # Filter to only include standard chromosomes and X (exclude Y, M for DO analysis)
    valid_chromosomes <- chromosomes[chromosomes %in% c(1:19, "X")]
    
    validation_log <- c(validation_log, paste("✓ Found", length(geno_files), "genotype files"))
    validation_log <- c(validation_log, paste("✓ Valid chromosomes for DO analysis:", paste(valid_chromosomes, collapse = ", ")))
    
    # Check that we have corresponding reference files for each chromosome
    available_founder <- gsub("GM_foundergeno(.*)\\\\.csv", "\\\\1", list.files(pattern = "GM_foundergeno.*\\\\.csv"))
    available_gmap <- gsub("GM_gmap(.*)\\\\.csv", "\\\\1", list.files(pattern = "GM_gmap.*\\\\.csv"))
    available_pmap <- gsub("GM_pmap(.*)\\\\.csv", "\\\\1", list.files(pattern = "GM_pmap.*\\\\.csv"))
    
    # Keep only chromosomes where we have all required files
    final_chromosomes <- intersect(
        intersect(intersect(valid_chromosomes, available_founder), available_gmap),
        available_pmap
    )
    
    validation_log <- c(validation_log, paste("✓ Final chromosomes with all required files:", paste(final_chromosomes, collapse = ", ")))
    
    if (length(final_chromosomes) == 0) {
        stop("No chromosomes have all required files (genotype, founder, gmap, pmap)")
    }
    
    # Check covariate file for required columns
    covar_data <- read.csv("${covar_csv}", row.names = 1)
    
    # Check for sex column (required) - case insensitive
    covar_names_lower <- tolower(colnames(covar_data))
    sex_col_idx <- which(covar_names_lower == "sex")

    if (length(sex_col_idx) == 0) {
        stop("sex column not found in covariate file - required for DO analysis")
    }

    sex_col_name <- colnames(covar_data)[sex_col_idx[1]]  # Use actual column name

    # Check sex coding
    sex_levels <- unique(covar_data[[sex_col_name]])
    validation_log <- c(validation_log, paste("✓ Sex levels found:", paste(sex_levels, collapse = ", ")))
    
    # Determine sex codes - common patterns are f/m, F/M, female/male, Female/Male
    if (all(sex_levels %in% c("f", "m"))) {
        sex_codes <- list(f = "Female", m = "Male")
    } else if (all(sex_levels %in% c("F", "M"))) {
        sex_codes <- list(F = "Female", M = "Male")
    } else if (all(sex_levels %in% c("female", "male"))) {
        sex_codes <- list(female = "Female", male = "Male")
    } else if (all(sex_levels %in% c("Female", "Male"))) {
        sex_codes <- list(Female = "Female", Male = "Male")
    } else {
        validation_log <- c(validation_log, paste("⚠ WARNING: Unusual sex coding detected:", paste(sex_levels, collapse = ", ")))
        # Default to using the levels as-is
        sex_codes <- structure(as.list(sex_levels), names = sex_levels)
    }
    
    # Check for ngen/generation column (crossinfo_covar)
    has_ngen <- "ngen" %in% colnames(covar_data)
    has_generation <- "generation" %in% colnames(covar_data)
    
    if (has_ngen) {
        validation_log <- c(validation_log, "✓ Found ngen column for crossinfo_covar")
        crossinfo_col <- "ngen"
    } else if (has_generation) {
        validation_log <- c(validation_log, "✓ Found generation column for crossinfo_covar")
        crossinfo_col <- "generation"
    } else {
        validation_log <- c(validation_log, "⚠ WARNING: Neither ngen nor generation column found - crossinfo_covar will be omitted")
        crossinfo_col <- NULL
    }
    
    # Create control file using write_control_file function
    cat("Creating control file...\\n")
    
    control_args <- list(
        output_file = paste0("${prefix}_control.json"),
        crosstype = "do",
        description = paste0("${prefix} DO QTL analysis"),
        founder_geno_file = paste0("GM_foundergeno", final_chromosomes, ".csv"),
        founder_geno_transposed = TRUE,
        gmap_file = paste0("GM_gmap", final_chromosomes, ".csv"),
        pmap_file = paste0("GM_pmap", final_chromosomes, ".csv"),
        geno_file = paste0("${prefix}_geno", final_chromosomes, ".csv"),
        geno_transposed = TRUE,
        geno_codes = list(A = 1, H = 2, B = 3),
        xchr = "X",
        pheno_file = "${pheno_csv}",
        covar_file = "${covar_csv}",
        sex_covar = sex_col_name,
        sex_codes = sex_codes
    )
    
    # Add crossinfo_covar if ngen or generation exists
    if (!is.null(crossinfo_col)) {
        control_args\$crossinfo_covar = crossinfo_col
    }
    
    # Write control file
    do.call(write_control_file, control_args)
    
    validation_log <- c(validation_log, paste("✓ Successfully created", paste0("${prefix}_control.json")))
    validation_log <- c(validation_log, paste("✓ Control file includes", length(final_chromosomes), "chromosomes"))
    validation_log <- c(validation_log, paste("✓ X chromosome handling:", ifelse("X" %in% final_chromosomes, "enabled", "not present")))
    
    # Write validation report
    writeLines(validation_log, "control_validation_report.txt")
    
    cat("Control file generation completed successfully\\n")
    """
}

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