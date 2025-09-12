process GENOTYPE_PROCESS {
    tag "Processing genotypes for ${prefix}"
    publishDir "${params.outdir}/02_genotype_processing", mode: 'copy'
    
    input:
    path(finalreport_file)
    val(prefix)
    
    output:
    path("${prefix}_geno*.csv"), emit: geno_files
    path("genotype_validation_report.txt"), emit: validation_report
    path("sample_id_mapping.txt"), emit: sample_mapping, optional: true
    path("genotype_summary.txt"), emit: summary
    path("GM_allelecodes.csv"), emit: allele_codes
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(qtl2)
        library(data.table)
    })
    
    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== Genotype Processing Validation Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, paste("FinalReport File:", "${finalreport_file}"))
    validation_log <- c(validation_log, "")
    
    # Download GigaMUGA reference files (build 39)
    cat("Downloading GigaMUGA reference files (build 39)...\\n")
    figshare_url <- "https://figshare.com/ndownloader/files/40233652"
    
    validation_log <- c(validation_log, "=== Reference Data Download ===")
    validation_log <- c(validation_log, paste("Downloading from:", figshare_url))
    
    tryCatch({
        download.file(figshare_url, 
                     destfile = "GM_processed_files_build39.zip", 
                     method = "auto",
                     mode = "wb")
        
        validation_log <- c(validation_log, "✓ Successfully downloaded GM_processed_files_build39.zip")
        
        # Extract the zip file
        unzip("GM_processed_files_build39.zip", exdir = ".")
        validation_log <- c(validation_log, "✓ Successfully extracted reference files")
        
        # Look for allele codes file in subdirectories
        allele_files <- list.files(pattern = ".*allelecodes\\\\.csv", recursive = TRUE, full.names = TRUE)
        
        if (length(allele_files) > 0) {
            # Copy the allele codes file to working directory
            file.copy(allele_files[1], "GM_allelecodes.csv")
            validation_log <- c(validation_log, paste("✓ Found allele codes file:", allele_files[1]))
        } else {
            # List all files to debug
            all_files <- list.files(recursive = TRUE)
            validation_log <- c(validation_log, "Files found in extraction:")
            validation_log <- c(validation_log, paste("  ", all_files))
            stop("No allele codes file found in extracted contents")
        }
        
        # Clean up zip file
        file.remove("GM_processed_files_build39.zip")
        
    }, error = function(e) {
        stop("Failed to download/extract reference files: ", e\$message)
    })
    
    # Validate input files exist
    if (!file.exists("${finalreport_file}")) {
        stop("FinalReport file not found: ${finalreport_file}")
    }
    
    if (!file.exists("GM_allelecodes.csv")) {
        stop("Allele codes file not found after download")
    }
    
    validation_log <- c(validation_log, "✓ Input files exist and are readable")
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Genotype Processing ===")
    
    # Read allele codes file (skip comment lines starting with #)
    allele_codes <- tryCatch({
        # Find the first non-comment line
        all_lines <- readLines("GM_allelecodes.csv")
        first_data_line <- which(!grepl("^#", all_lines))[1]
        
        read.csv("GM_allelecodes.csv", 
                stringsAsFactors = FALSE,
                skip = first_data_line - 1,
                header = TRUE,
                comment.char = "#")
    }, error = function(e) {
        stop("Error reading allele codes file: ", e\$message)
    })
    
    validation_log <- c(validation_log, paste("✓ Loaded allele codes:", nrow(allele_codes), "markers"))
    
    # Test mode: filter allele codes to chromosome 19 only BEFORE merging
    if ("${params.test_mode}" == "true") {
        allele_codes <- allele_codes[allele_codes\$chr == "19", ]
        validation_log <- c(validation_log, "⚠ TEST MODE: Filtered allele codes to chromosome 19 only")
        validation_log <- c(validation_log, paste("  Reduced to", nrow(allele_codes), "markers for processing"))
    }
    
    # Detect FinalReport format and read appropriately
    sample_lines <- readLines("${finalreport_file}", n = 20)
    data_start <- which(grepl("^\\\\[Data\\\\]", sample_lines))
    
    if (length(data_start) == 0) {
        data_start <- which(grepl("SNP.Name|SNP_Name", sample_lines))
        if (length(data_start) == 0) {
            data_start <- 1
            validation_log <- c(validation_log, "⚠ WARNING: Could not detect header structure, assuming data starts at line 1")
        } else {
            data_start <- data_start[1]
            validation_log <- c(validation_log, paste("✓ Detected header structure, data starts at line", data_start))
        }
    } else {
        data_start <- data_start[1] + 1
        validation_log <- c(validation_log, paste("✓ Found [Data] section at line", data_start-1))
    }
    
    # Read FinalReport data using data.table for memory efficiency
    cat("Reading FinalReport data...\\n")
    geno_data <- tryCatch({
        if (data_start > 1) {
            fread("${finalreport_file}", 
                  skip = data_start - 1,
                  header = TRUE,
                  sep = "\\t",
                  stringsAsFactors = FALSE,
                  data.table = FALSE)
        } else {
            fread("${finalreport_file}", 
                  header = TRUE,
                  sep = "\\t", 
                  stringsAsFactors = FALSE,
                  data.table = FALSE)
        }
    }, error = function(e) {
        stop("Error reading FinalReport data: ", e\$message)
    })
    
    validation_log <- c(validation_log, paste("✓ Successfully read FinalReport data"))
    validation_log <- c(validation_log, paste("  - Dimensions:", nrow(geno_data), "rows x", ncol(geno_data), "columns"))
    
    # Standardize column names
    colnames(geno_data) <- gsub("\\\\.", "_", colnames(geno_data))
    colnames(geno_data) <- gsub(" - ", "_", colnames(geno_data))
    colnames(geno_data) <- gsub(" ", "_", colnames(geno_data))
    
    # Identify key columns
    sample_col <- "Sample_ID"
    snp_col <- "SNP_Name"
    allele1_col <- "Allele1_Forward"
    allele2_col <- "Allele2_Forward"
    
    # Verify columns exist
    required_cols <- c(sample_col, snp_col, allele1_col, allele2_col)
    missing_cols <- setdiff(required_cols, colnames(geno_data))
    
    if (length(missing_cols) > 0) {
        validation_log <- c(validation_log, "Available columns:")
        validation_log <- c(validation_log, paste("  ", colnames(geno_data)))
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    validation_log <- c(validation_log, "✓ Successfully identified required columns")
    
    # Test mode: filter FinalReport data to only chromosome 19 markers BEFORE merging
    if ("${params.test_mode}" == "true") {
        # Get list of chromosome 19 markers from allele codes
        chr_19_markers <- allele_codes[allele_codes\$chr == "19", "marker"]
        
        # Filter FinalReport data to only these markers
        geno_data <- geno_data[geno_data[[snp_col]] %in% chr_19_markers, ]
        
        validation_log <- c(validation_log, paste("⚠ TEST MODE: Filtered FinalReport to", nrow(geno_data), "genotype calls for chromosome 19"))
    }
    
    # Get unique samples and SNPs
    unique_samples <- unique(geno_data[[sample_col]])
    unique_snps <- unique(geno_data[[snp_col]])
    
    validation_log <- c(validation_log, paste("✓ Found", length(unique_samples), "unique samples"))
    validation_log <- c(validation_log, paste("✓ Found", length(unique_snps), "unique SNPs"))
    
    # Apply sample ID transformations if needed
    if ("${params.auto_prefix_samples}" == "true") {
        geno_data[[sample_col]] <- paste0("${prefix}_", geno_data[[sample_col]])
        validation_log <- c(validation_log, "✓ Applied automatic sample ID prefixing")
    } else {
        validation_log <- c(validation_log, "✓ Kept original sample IDs")
    }
    
    # Merge with allele codes
    marker_col <- colnames(allele_codes)[1]
    geno_merged <- merge(geno_data, allele_codes, 
                        by.x = snp_col, by.y = marker_col, 
                        all.x = TRUE)
    
    validation_log <- c(validation_log, paste("✓ Merged with allele codes, retained", nrow(geno_merged), "genotype calls"))
    
    # Find chromosome column
    chr_col <- NULL
    chr_patterns <- c("chr", "chromosome", "CHR", "Chromosome")
    
    for (pattern in chr_patterns) {
        if (pattern %in% colnames(geno_merged)) {
            chr_col <- pattern
            break
        }
    }
    
    if (is.null(chr_col)) {
        validation_log <- c(validation_log, "Available columns in merged data:")
        validation_log <- c(validation_log, paste("  ", colnames(geno_merged)))
        stop("Could not find chromosome information in allele codes file")
    }
    
    validation_log <- c(validation_log, paste("✓ Found chromosome column:", chr_col))
    
    # Convert genotype calls to A/H/B format
    cat("Converting genotype calls to A/H/B format...\\n")
    
    geno_merged\$genotype_call <- paste0(geno_merged[[allele1_col]], geno_merged[[allele2_col]])
    geno_merged\$genotype_call[geno_merged[[allele1_col]] == "-" | geno_merged[[allele2_col]] == "-"] <- "-"
    geno_merged\$genotype_call[is.na(geno_merged[[allele1_col]]) | is.na(geno_merged[[allele2_col]])] <- "-"
    
    # Convert to A/H/B based on founder alleles
    if ("A" %in% colnames(geno_merged) && "B" %in% colnames(geno_merged)) {
        geno_merged\$qtl2_call <- "-"
        
        aa_pattern <- paste0(geno_merged\$A, geno_merged\$A)
        bb_pattern <- paste0(geno_merged\$B, geno_merged\$B)
        ab_pattern1 <- paste0(geno_merged\$A, geno_merged\$B)
        ab_pattern2 <- paste0(geno_merged\$B, geno_merged\$A)
        
        geno_merged\$qtl2_call[geno_merged\$genotype_call == aa_pattern] <- "A"
        geno_merged\$qtl2_call[geno_merged\$genotype_call == bb_pattern] <- "B"
        geno_merged\$qtl2_call[geno_merged\$genotype_call == ab_pattern1 | geno_merged\$genotype_call == ab_pattern2] <- "H"
        
        validation_log <- c(validation_log, "✓ Converted genotype calls to A/H/B format")
    } else {
        validation_log <- c(validation_log, "⚠ WARNING: Using simplified genotype conversion")
        geno_merged\$qtl2_call <- "H"
    }
    
    # Create chromosome-specific files using DO_pipe matrix approach
    chromosomes <- sort(unique(geno_merged[[chr_col]]))
    validation_log <- c(validation_log, paste("✓ Found chromosomes:", paste(chromosomes, collapse = ", ")))
    
    cat("Creating chromosome-specific files...\\n")
    
    files_created <- c()
    for (chr in chromosomes) {
        chr_data <- geno_merged[geno_merged[[chr_col]] == chr, ]
        
        if (nrow(chr_data) == 0) next
        
        # Get unique markers and samples for this chromosome
        chr_markers <- unique(chr_data[[snp_col]])
        chr_samples <- unique(chr_data[[sample_col]])
        
        # Create matrix like DO_pipe (handles duplicates by overwriting)
        geno_matrix <- matrix(nrow = length(chr_markers), ncol = length(chr_samples))
        dimnames(geno_matrix) <- list(chr_markers, chr_samples)
        
        # Fill matrix - this handles duplicates by overwriting like DO_pipe
        for (i in seq_along(chr_samples)) {
            sample_data <- chr_data[chr_data[[sample_col]] == chr_samples[i], ]
            geno_matrix[sample_data[[snp_col]], i] <- sample_data\$qtl2_call
        }
        
        # Replace NAs with "-"
        geno_matrix[is.na(geno_matrix)] <- "-"
        
        chr_filename <- paste0("${prefix}_geno", chr, ".csv")
        write.csv(geno_matrix, chr_filename, row.names = TRUE)
        files_created <- c(files_created, chr_filename)
        
        validation_log <- c(validation_log, paste("✓ Created", chr_filename, "with", nrow(geno_matrix), "markers"))
    }
    
    # Generate summary
    summary_info <- c()
    summary_info <- c(summary_info, "=== Genotype Processing Summary ===")
    summary_info <- c(summary_info, paste("Study Prefix:", "${prefix}"))
    summary_info <- c(summary_info, paste("Total Samples:", length(unique_samples)))
    summary_info <- c(summary_info, paste("Total SNPs:", length(unique_snps)))
    summary_info <- c(summary_info, paste("Files Created:", length(files_created)))
    
    writeLines(validation_log, "genotype_validation_report.txt")
    writeLines(summary_info, "genotype_summary.txt")
    
    cat("Genotype processing completed successfully\\n")
    """
}