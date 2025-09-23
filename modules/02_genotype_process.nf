process GENOTYPE_PROCESS {
    tag "Processing genotypes for ${prefix}"
    publishDir "${params.outdir}/02_genotype_processing", mode: 'copy'
    
    input:
    path(finalreport_files)
    path(valid_samples_file)
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
validation_log <- c(validation_log, paste("FinalReport Files:", "${finalreport_files}"))
validation_log <- c(validation_log, paste("Valid Samples File:", "${valid_samples_file}"))
validation_log <- c(validation_log, "")

# Read valid sample list from Module 1
if (!file.exists("${valid_samples_file}")) {
    stop("Valid samples file not found: ${valid_samples_file}")
}

valid_samples <- readLines("${valid_samples_file}")
validation_log <- c(validation_log, "=== Sample Intersection Logic ===")
validation_log <- c(validation_log, paste("✓ Loaded valid sample list:", length(valid_samples), "samples"))
validation_log <- c(validation_log, paste("✓ Sample filtering will be applied to FinalReport data"))
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
    cat("DEBUG: Searching for allele codes file with pattern: '.*allelecodes\\\\.csv'\\n")
    allele_files <- list.files(pattern = ".*allelecodes\\\\.csv", recursive = TRUE, full.names = TRUE)
    
    if (length(allele_files) > 0) {
        # Copy the allele codes file to working directory
        file.copy(allele_files[1], "GM_allelecodes.csv")
        validation_log <- c(validation_log, paste("✓ Found allele codes file:", allele_files[1]))
    } else {
        stop("No allele codes file found in extracted contents")
    }
    
    # Clean up zip file
    file.remove("GM_processed_files_build39.zip")
    
}, error = function(e) {
    stop("Failed to download/extract reference files: ", e\$message)
})

# Validate input files exist
finalreport_list <- strsplit("${finalreport_files}", " ")[[1]]
for (fr_file in finalreport_list) {
    if (!file.exists(fr_file)) {
        stop(paste("FinalReport file not found:", fr_file))
    }
}
validation_log <- c(validation_log, paste("Found", length(finalreport_list), "FinalReport file(s)"))

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

# Test mode: filter allele codes to chromosome 2 only BEFORE merging
if ("${params.test_mode}" == "true") {
    allele_codes <- allele_codes[allele_codes\$chr == "2", ]
    validation_log <- c(validation_log, "⚠ TEST MODE: Filtered allele codes to chromosome 2 only")
    validation_log <- c(validation_log, paste("  Reduced to", nrow(allele_codes), "markers for processing"))
}

# --- CRITICAL FIX START ---
# Replicating the DO_pipe approach for reading and filtering data
cat("Reading FinalReport data from multiple files...\\n")

# Initialize empty data frame for combined genotype data
geno_data <- NULL

# Process each FinalReport file
for (fr_file in finalreport_list) {
    cat(paste("Processing:", fr_file, "\\n"))

    # Detect FinalReport format and read appropriately
    sample_lines <- readLines(fr_file, n = 20)
    data_start <- which(grepl("SNP.Name|SNP_Name", sample_lines))

    if (length(data_start) == 0) {
        stop(paste("Could not detect header line in FinalReport file:", fr_file))
    }
    data_start <- data_start[1]

    cat("DEBUG: FinalReport data starts at line:", data_start, "\\n")
    cat("DEBUG: Reading FinalReport data with fread. Sep value is '\\\\t'\\n")
    current_data <- fread(fr_file,
                          skip = data_start - 1,
                          header = TRUE,
                          sep = "\\t",
                          stringsAsFactors = FALSE,
                          data.table = FALSE,
                          fill = TRUE)

    # Combine data from multiple files
    if (is.null(geno_data)) {
        geno_data <- current_data
    } else {
        # Bind rows, assuming same column structure
        geno_data <- rbind(geno_data, current_data)
    }

    cat(paste("  Added", nrow(current_data), "rows from", fr_file, "\\n"))
}

validation_log <- c(validation_log, paste("✓ Successfully read and combined FinalReport data from", length(finalreport_list), "file(s)"))
validation_log <- c(validation_log, paste("  - Dimensions:", nrow(geno_data), "rows x", ncol(geno_data), "columns"))

# Standardize column names
cat("DEBUG: Column names BEFORE gsub:", paste(colnames(geno_data), collapse=", "), "\\n")
cat("DEBUG: Standardizing column names (period to underscore)...\\n")
colnames(geno_data) <- gsub("\\\\.", "_", colnames(geno_data))
cat("DEBUG: Standardizing column names (' - ' to underscore)...\\n")
colnames(geno_data) <- gsub(" - ", "_", colnames(geno_data))
cat("DEBUG: Standardizing column names (space to underscore)...\\n")
colnames(geno_data) <- gsub(" ", "_", colnames(geno_data))
cat("DEBUG: Column names AFTER gsub:", paste(colnames(geno_data), collapse=", "), "\\n")

# Identify key columns
sample_col <- "Sample_ID"
snp_col <- "SNP_Name"
allele1_col <- "Allele1_Forward"
allele2_col <- "Allele2_Forward"

# Verify columns exist
required_cols <- c(sample_col, snp_col, allele1_col, allele2_col)
missing_cols <- setdiff(required_cols, colnames(geno_data))

if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

validation_log <- c(validation_log, "✓ Successfully identified required columns")

# Use DO_pipe strategy: get all unique samples from the combined FinalReport data
original_samples <- unique(geno_data[[sample_col]])
validation_log <- c(validation_log, paste("✓ Combined FinalReport data contains", length(original_samples), "unique samples"))

# Filter to get a final, consistent list of samples
# This list is the single source of truth for matrix columns
all_samples <- intersect(as.character(original_samples), as.character(valid_samples))

if (length(all_samples) == 0) {
    stop("No samples remain after intersection filtering! Check sample ID consistency between phenotype and genotype files.")
}

# Now, filter the data frame ONCE based on this list
geno_data <- geno_data[geno_data[[sample_col]] %in% all_samples, ]

filtered_samples <- unique(geno_data[[sample_col]])
samples_removed <- length(original_samples) - length(filtered_samples)

validation_log <- c(validation_log, paste("✓ Applied sample intersection filtering"))
validation_log <- c(validation_log, paste("  - Original samples in combined FinalReport data:", length(original_samples)))
validation_log <- c(validation_log, paste("  - Valid samples from phenotype file:", length(valid_samples)))
validation_log <- c(validation_log, paste("  - Samples after filtering:", length(filtered_samples)))
validation_log <- c(validation_log, paste("  - Samples removed:", samples_removed))

# Get all unique markers from the allele codes file (the master list of markers)
all_markers <- allele_codes\$marker

# Create an empty matrix for ALL markers and ALL filtered samples (DO_pipe style)
geno_matrix <- matrix(NA, nrow = length(all_markers), ncol = length(all_samples))
dimnames(geno_matrix) <- list(all_markers, all_samples)
validation_log <- c(validation_log, paste("✓ Created master genotype matrix:", nrow(geno_matrix), "x", ncol(geno_matrix)))

# --- CRITICAL FIX: Loop through samples to fill the master matrix ---
cat("Filling master matrix by sample...\\n")
for(i in seq_along(all_samples)) {
    sample_id <- all_samples[i]
    # Get all rows for this specific sample
    wh <- (geno_data[[sample_col]] == sample_id)

    if (any(wh)) {
        # Get SNP names and genotype calls for this sample
        snp_names <- geno_data[[snp_col]][wh]
        allele1 <- geno_data[[allele1_col]][wh]
        allele2 <- geno_data[[allele2_col]][wh]

        # Combine alleles into a single string (like DO_pipe)
        genotype_call <- paste0(allele1, allele2)
        
        # Only assign SNPs that exist in the master matrix rows
        # This handles cases where some markers are missing for a sample
        valid_snps <- snp_names %in% rownames(geno_matrix)

        if (any(valid_snps)) {
            # Assign the calls to the correct row (SNP name) and column (sample index i)
            geno_matrix[snp_names[valid_snps], i] <- genotype_call[valid_snps]
        }
    }
}
validation_log <- c(validation_log, "✓ Successfully populated master genotype matrix")

# --- CONVERT TO qtl2 FORMAT ---
cat("Converting to qtl2 A/H/B format...\\n")

# Ensure allele_codes and geno_matrix have the same marker order
# The dimnames() call above should ensure this, but a re-check is good.
# Let's reorder allele_codes to match geno_matrix row order
allele_codes_ordered <- allele_codes[match(rownames(geno_matrix), allele_codes\$marker), ]

# Use qtl2convert::encode_geno on the filled matrix
geno_matrix_qtl2 <- qtl2convert::encode_geno(geno_matrix, as.matrix(allele_codes_ordered[, c("A", "B")]))
validation_log <- c(validation_log, "✓ Converted genotypes to A/H/B calls")

# --- SPLIT AND WRITE FILES ---
cat("Splitting and writing chromosome-specific files...\\n")
files_created <- c()

# We need the chromosome information from the ordered allele_codes
chromosomes <- sort(unique(allele_codes_ordered\$chr))

for (chr in chromosomes) {
    markers_for_chr <- allele_codes_ordered\$marker[allele_codes_ordered\$chr == chr]
    
    if (length(markers_for_chr) == 0) next
    
    # Extract the subset of the master matrix for this chromosome
    g_chr <- geno_matrix_qtl2[markers_for_chr, , drop = FALSE]

    chr_filename <- paste0("${prefix}_geno", chr, ".csv")
    
    # Prepare data for write.csv
    g_chr_out <- as.data.frame(g_chr)
    g_chr_out\$marker <- rownames(g_chr_out)
    g_chr_out <- g_chr_out[, c("marker", all_samples)]
    
    write.csv(g_chr_out, chr_filename, row.names = FALSE)
    files_created <- c(files_created, chr_filename)
    
    validation_log <- c(validation_log, paste("✓ Created", chr_filename, "with", nrow(g_chr), "markers"))
}

# Generate summary
summary_info <- c()
summary_info <- c(summary_info, "=== Genotype Processing Summary ===")
summary_info <- c(summary_info, paste("Study Prefix:", "${prefix}"))
summary_info <- c(summary_info, paste("Total Samples:", length(all_samples)))
summary_info <- c(summary_info, paste("Total SNPs:", length(all_markers)))
summary_info <- c(summary_info, paste("Files Created:", length(files_created)))

writeLines(validation_log, "genotype_validation_report.txt")
writeLines(summary_info, "genotype_summary.txt")

cat("Genotype processing completed successfully\\n")
    """
}