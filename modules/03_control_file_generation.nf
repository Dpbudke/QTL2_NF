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