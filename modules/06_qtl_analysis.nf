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