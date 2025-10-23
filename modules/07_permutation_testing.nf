#!/usr/bin/env nextflow

// Simplified Single-Phenotype Permutation Testing Module
// New Architecture: 1 batch = 1 phenotype × 200 permutations = 1 SLURM job
// Benefits: Real-time progress monitoring, faster failure detection, simpler code, optimized resource usage

process PERM_SETUP {
    tag "Setting up permutation testing"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '10m'

    input:
    path(cross2_file)
    path(filtered_phenotypes_file)
    val(study_prefix)

    output:
    path("${study_prefix}_phenotype_list.txt"), emit: phenotype_list
    path("${study_prefix}_filtered_cross2.rds"), emit: filtered_cross2
    path("${study_prefix}_setup_log.txt"), emit: setup_log

    script:
    """
    #!/usr/bin/env Rscript

    log_file <- "${study_prefix}_setup_log.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message("=== SIMPLIFIED PERMUTATION SETUP ===")
    log_message("New Architecture: 1 batch = 1 phenotype × 200 perms")
    log_message("Loading data...")

    # Load full cross2 object
    cross2_full <- readRDS("${cross2_file}")

    # Filter phenotypes based on Module 6 LOD threshold results
    filtered_phenos <- readLines("${filtered_phenotypes_file}")
    log_message(paste("Original phenotypes:", ncol(cross2_full\$pheno)))
    log_message(paste("Filtered phenotypes:", length(filtered_phenos)))

    # Subset cross2 object to only include filtered phenotypes
    pheno_indices <- which(colnames(cross2_full\$pheno) %in% filtered_phenos)
    if(length(pheno_indices) == 0) {
        log_message("ERROR: No phenotypes match between cross2 and filtered list")
        stop("No matching phenotypes found")
    }

    cross2 <- cross2_full
    cross2\$pheno <- cross2_full\$pheno[, pheno_indices, drop = FALSE]

    log_message(paste("Proceeding with", ncol(cross2\$pheno), "filtered phenotypes"))

    # Save filtered cross2 for downstream processes
    saveRDS(cross2, file = "${study_prefix}_filtered_cross2.rds")

    # Calculate batch structure
    num_phenotypes <- ncol(cross2\$pheno)
    perms_per_batch <- 200  # Fixed: 200 perms per batch (4x original 50)
    batches_per_pheno <- 5   # Fixed: 5 batches per phenotype (200 × 5 = 1000 total)
    total_batches <- num_phenotypes * batches_per_pheno

    log_message("")
    log_message("=== BATCH CONFIGURATION ===")
    log_message(paste("Phenotypes to process:", num_phenotypes))
    log_message(paste("Permutations per batch:", perms_per_batch))
    log_message(paste("Batches per phenotype:", batches_per_pheno))
    log_message(paste("Total batches:", total_batches))
    log_message(paste("Total permutations:", num_phenotypes * 1000))
    log_message("")
    log_message("=== RESOURCE ALLOCATION ===")
    log_message(paste("CPUs per batch: 48 (fork-based parallelism)"))
    log_message(paste("Memory per batch: 200 GB (optimized for 200 perms, ~164GB avg usage)"))
    log_message(paste("Max parallel batches: ~35-40 (observed on ceres cluster)"))
    log_message(paste("Estimated time per batch: 6-8 minutes (200 perms each, ~2x 100-perm time)"))
    log_message(paste("Estimated total time:", round(total_batches * 7 / 60 / 35, 1), "hours (with 35 parallel)"))
    log_message("")

    # Create phenotype list (write each phenotype name for iteration)
    writeLines(colnames(cross2\$pheno), "${study_prefix}_phenotype_list.txt")

    log_message("=== SETUP COMPLETE ===")
    log_message(paste("Ready to launch", total_batches, "batch jobs"))
    log_message("Progress will be visible as batches complete")

    close(log_conn)
    """
}

process PERM_BATCH_JOB {
    tag "Perm: \${phenotype_name} batch \${batch_num}/5"
    publishDir "${params.outdir}/07_permutation_testing/batches", mode: 'copy'

    // Single phenotype, single batch
    // 48 CPUs and 200GB per batch (optimized for 200 perms per batch, ~164GB avg usage)

    input:
    val(study_prefix)
    val(phenotype_name)
    val(batch_num)  // 1-5 (200 perms per batch)

    output:
    path("${study_prefix}_${phenotype_name}_batch_${batch_num}.rds"), emit: perm_result

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(qtl2)
    })

    batch_start_time <- Sys.time()

    cat("===========================================\\n")
    cat("BATCH:", "${batch_num}", "/ 5\\n")
    cat("PHENOTYPE:", "${phenotype_name}\\n")
    cat("PERMUTATIONS: 200\\n")
    cat("START TIME:", format(batch_start_time), "\\n")
    cat("===========================================\\n\\n")

    # Load data from published results directories
    cat("Loading data files...\\n")
    load_start <- Sys.time()

    project_dir <- "${projectDir}"
    cross2 <- readRDS(file.path(project_dir, "${params.outdir}/07_permutation_testing/${study_prefix}_filtered_cross2.rds"))
    alleleprob <- readRDS(file.path(project_dir, "${params.outdir}/05_genome_scan_preparation/${study_prefix}_alleleprob.rds"))
    kinship_loco <- readRDS(file.path(project_dir, "${params.outdir}/05_genome_scan_preparation/${study_prefix}_kinship_loco.rds"))

    load_time <- as.numeric(difftime(Sys.time(), load_start, units = "secs"))
    cat("Data loaded in", round(load_time, 1), "seconds\\n\\n")

    # Extract additive covariates from cross2 object
    if(!is.null(cross2\$covar)) {
        covar_data <- cross2\$covar
        if("coat_color" %in% colnames(covar_data)) {
            covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
        }
        covar_formula <- paste("~", paste(colnames(covar_data), collapse = " + "))
        addcovar <- model.matrix(as.formula(covar_formula), data = covar_data)[, -1, drop = FALSE]
    } else {
        addcovar <- NULL
    }

    # Get X chromosome covariates
    Xcovar <- NULL
    if("X" %in% names(cross2\$gmap)) {
        tryCatch({
            Xcovar <- get_x_covar(cross2)
        }, error = function(e) {
            cat("Warning: Could not generate X covariates:", e\$message, "\\n")
        })
    }

    # Prepare autosomes only
    is_x_chr <- attr(alleleprob, "is_x_chr")
    if (!is.null(is_x_chr)) {
        kinship_auto <- kinship_loco[names(is_x_chr)[!is_x_chr]]
        alleleprob_auto <- alleleprob[,!is_x_chr]
    } else {
        kinship_auto <- kinship_loco
        alleleprob_auto <- alleleprob
    }

    # Extract single phenotype
    pheno_data <- cross2\$pheno[, "${phenotype_name}", drop = FALSE]

    cat("Running permutation test...\\n")
    cat("Phenotype:", "${phenotype_name}\\n")
    cat("Permutations: 200\\n")
    cat("CPUs: 48 (fork-based parallelism with shared memory)\\n\\n")

    perm_start <- Sys.time()

    # Run permutations with fork-based parallelism (cores=48)
    # Fork shares memory on Linux - no data serialization like PSOCK
    # All 48 cores share the same alleleprob/kinship data in memory
    cat("Starting scan1perm() with 200 permutations using 48 cores (fork-based)...\\n")
    flush.console()

    perm <- scan1perm(alleleprob_auto,
                      pheno_data,
                      kinship = kinship_auto,
                      addcovar = addcovar,
                      Xcovar = Xcovar,
                      cores = 48,  # Fork-based parallelism - shares memory
                      n_perm = 200)

    cat("✓ scan1perm() completed successfully\\n")

    perm_time <- as.numeric(difftime(Sys.time(), perm_start, units = "secs"))

    # Save result
    out <- data.frame(perm, check.names = FALSE)
    saveRDS(out, file = "${study_prefix}_${phenotype_name}_batch_${batch_num}.rds")

    batch_end_time <- Sys.time()
    total_time <- as.numeric(difftime(batch_end_time, batch_start_time, units = "secs"))

    cat("\\n===========================================\\n")
    cat("BATCH COMPLETE\\n")
    cat("===========================================\\n")
    cat("Phenotype:", "${phenotype_name}\\n")
    cat("Batch:", "${batch_num}", "/ 5\\n")
    cat("Total time:", round(total_time / 60, 2), "minutes\\n")
    cat("Permutation time:", round(perm_time / 60, 2), "minutes\\n")
    cat("Data loading time:", round(load_time, 1), "seconds\\n")
    cat("End time:", format(batch_end_time), "\\n")
    cat("===========================================\\n")
    """
}

process PERM_VALIDATE_BATCHES {
    tag "Validating permutation batches"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 4
    memory '16 GB'
    time '30m'

    input:
    path(phenotype_list)
    val(study_prefix)

    output:
    path("${study_prefix}_batch_validation_report.txt"), emit: validation_report
    path("${study_prefix}_corrupt_batches.txt"), emit: corrupt_batches optional true
    path("${study_prefix}_missing_batches.txt"), emit: missing_batches optional true
    env VALIDATION_STATUS, emit: status

    script:
    """
    #!/usr/bin/env Rscript

    log_file <- "${study_prefix}_batch_validation_report.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message("=== PERMUTATION BATCH VALIDATION ===")

    # Read expected phenotypes
    phenotypes <- readLines("${phenotype_list}")
    num_phenotypes <- length(phenotypes)
    batches_per_pheno <- 5
    expected_total <- num_phenotypes * batches_per_pheno

    log_message(paste("Expected phenotypes:", num_phenotypes))
    log_message(paste("Batches per phenotype:", batches_per_pheno))
    log_message(paste("Expected total batches:", expected_total))
    log_message("")

    # Get batch directory from project directory
    project_dir <- "${projectDir}"
    batch_dir <- file.path(project_dir, "${params.outdir}/07_permutation_testing/batches")

    log_message(paste("Checking batch directory:", batch_dir))

    if (!dir.exists(batch_dir)) {
        log_message("ERROR: Batch directory does not exist!")
        writeLines("FAILED", "VALIDATION_STATUS")
        close(log_conn)
        quit(status = 1)
    }

    # Get all batch files
    batch_files <- list.files(batch_dir, pattern = "${study_prefix}_.*_batch_[0-9]+\\\\.rds\$", full.names = TRUE)
    log_message(paste("Found", length(batch_files), "batch files"))
    log_message("")

    # Track missing and corrupt batches
    missing_batches <- data.frame(
        phenotype = character(),
        batch = integer(),
        stringsAsFactors = FALSE
    )

    corrupt_batches <- data.frame(
        phenotype = character(),
        batch = integer(),
        file = character(),
        size = integer(),
        error = character(),
        stringsAsFactors = FALSE
    )

    # Check each phenotype has 5 batches
    log_message("=== CHECKING FOR MISSING BATCHES ===")
    for (pheno in phenotypes) {
        for (batch_num in 1:5) {
            expected_file <- file.path(batch_dir, paste0("${study_prefix}_", pheno, "_batch_", batch_num, ".rds"))
            if (!file.exists(expected_file)) {
                missing_batches <- rbind(missing_batches, data.frame(
                    phenotype = pheno,
                    batch = batch_num,
                    stringsAsFactors = FALSE
                ))
                log_message(paste("MISSING:", basename(expected_file)))
            }
        }
    }

    if (nrow(missing_batches) > 0) {
        log_message(paste("Total missing batches:", nrow(missing_batches)))
    } else {
        log_message("✓ All expected batch files exist")
    }
    log_message("")

    # Check all batch files for corruption
    log_message("=== CHECKING FOR CORRUPT BATCHES ===")
    log_message("Testing readability of all batch files...")

    corrupt_count <- 0
    pb <- txtProgressBar(min = 0, max = length(batch_files), style = 3)

    for (i in seq_along(batch_files)) {
        bf <- batch_files[i]
        file_size <- file.size(bf)
        file_name <- basename(bf)

        # Extract phenotype and batch from filename
        pheno <- sub("${study_prefix}_(.*)_batch_[0-9]+\\\\.rds", "\\\\1", file_name)
        batch_num <- as.integer(sub(".*_batch_([0-9]+)\\\\.rds", "\\\\1", file_name))

        # Check for suspiciously small files
        if (file_size < 1000) {
            corrupt_count <- corrupt_count + 1
            corrupt_batches <- rbind(corrupt_batches, data.frame(
                phenotype = pheno,
                batch = batch_num,
                file = file_name,
                size = file_size,
                error = paste("File too small:", file_size, "bytes"),
                stringsAsFactors = FALSE
            ))
            log_message(paste("CORRUPT (too small):", file_name, "-", file_size, "bytes"))
            next
        }

        # Try to read the file
        tryCatch({
            obj <- readRDS(bf)
            # Verify it's the correct structure
            if (!is.data.frame(obj) && !is.matrix(obj)) {
                corrupt_count <- corrupt_count + 1
                corrupt_batches <- rbind(corrupt_batches, data.frame(
                    phenotype = pheno,
                    batch = batch_num,
                    file = file_name,
                    size = file_size,
                    error = "Invalid object type",
                    stringsAsFactors = FALSE
                ))
                log_message(paste("CORRUPT (invalid type):", file_name))
            }
        }, error = function(e) {
            corrupt_count <- corrupt_count + 1
            corrupt_batches <- rbind(corrupt_batches, data.frame(
                phenotype = pheno,
                batch = batch_num,
                file = file_name,
                size = file_size,
                error = e\$message,
                stringsAsFactors = FALSE
            ))
            log_message(paste("CORRUPT (read error):", file_name, "-", e\$message))
        })

        setTxtProgressBar(pb, i)
    }
    close(pb)

    if (corrupt_count > 0) {
        log_message(paste("Total corrupt batches:", corrupt_count))
    } else {
        log_message("✓ All batch files are valid and readable")
    }
    log_message("")

    # Summary
    log_message("=== VALIDATION SUMMARY ===")
    log_message(paste("Expected batches:", expected_total))
    log_message(paste("Found batches:", length(batch_files)))
    log_message(paste("Missing batches:", nrow(missing_batches)))
    log_message(paste("Corrupt batches:", nrow(corrupt_batches)))

    total_issues <- nrow(missing_batches) + nrow(corrupt_batches)

    if (total_issues > 0) {
        log_message("")
        log_message(paste("VALIDATION FAILED:", total_issues, "batches need repair"))

        # Save lists for repair process
        if (nrow(missing_batches) > 0) {
            write.table(missing_batches, "${study_prefix}_missing_batches.txt",
                       sep = "\\t", quote = FALSE, row.names = FALSE)
            log_message("Missing batches list saved")
        }

        if (nrow(corrupt_batches) > 0) {
            write.table(corrupt_batches, "${study_prefix}_corrupt_batches.txt",
                       sep = "\\t", quote = FALSE, row.names = FALSE)
            log_message("Corrupt batches list saved")
        }

        writeLines("FAILED", "VALIDATION_STATUS")
    } else {
        log_message("")
        log_message("✓ VALIDATION PASSED: All batches present and valid")
        writeLines("PASSED", "VALIDATION_STATUS")
    }

    close(log_conn)
    """
}

process PERM_REPAIR_BATCH {
    tag "Repair: \${phenotype} batch \${batch_num}"
    publishDir "${params.outdir}/07_permutation_testing/batches", mode: 'copy', overwrite: true

    cpus 48
    memory '200 GB'
    time '2h'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val(study_prefix)
    val(phenotype)
    val(batch_num)

    output:
    path("${study_prefix}_${phenotype}_batch_${batch_num}.rds"), emit: repaired_batch

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(qtl2)
    })

    batch_start_time <- Sys.time()

    cat("===========================================\\n")
    cat("REPAIR BATCH JOB\\n")
    cat("===========================================\\n")
    cat("PHENOTYPE:", "${phenotype}", "\\n")
    cat("BATCH:", "${batch_num}", "/ 5\\n")
    cat("PERMUTATIONS: 200\\n")
    cat("START TIME:", format(batch_start_time), "\\n")
    cat("===========================================\\n\\n")

    # Load data from published results directories
    cat("Loading data files...\\n")
    load_start <- Sys.time()

    project_dir <- "${projectDir}"
    cross2 <- readRDS(file.path(project_dir, "${params.outdir}/07_permutation_testing/${study_prefix}_filtered_cross2.rds"))
    alleleprob <- readRDS(file.path(project_dir, "${params.outdir}/05_genome_scan_preparation/${study_prefix}_alleleprob.rds"))
    kinship_loco <- readRDS(file.path(project_dir, "${params.outdir}/05_genome_scan_preparation/${study_prefix}_kinship_loco.rds"))

    load_time <- as.numeric(difftime(Sys.time(), load_start, units = "secs"))
    cat("Data loaded in", round(load_time, 1), "seconds\\n\\n")

    # Extract additive covariates from cross2 object
    if(!is.null(cross2\$covar)) {
        covar_data <- cross2\$covar
        if("coat_color" %in% colnames(covar_data)) {
            covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
        }
        covar_formula <- paste("~", paste(colnames(covar_data), collapse = " + "))
        addcovar <- model.matrix(as.formula(covar_formula), data = covar_data)[, -1, drop = FALSE]
    } else {
        addcovar <- NULL
    }

    # Get X chromosome covariates
    Xcovar <- NULL
    if("X" %in% names(cross2\$gmap)) {
        tryCatch({
            Xcovar <- get_x_covar(cross2)
        }, error = function(e) {
            cat("Warning: Could not generate X covariates:", e\$message, "\\n")
        })
    }

    # Prepare autosomes only
    is_x_chr <- attr(alleleprob, "is_x_chr")
    if (!is.null(is_x_chr)) {
        kinship_auto <- kinship_loco[names(is_x_chr)[!is_x_chr]]
        alleleprob_auto <- alleleprob[,!is_x_chr]
    } else {
        kinship_auto <- kinship_loco
        alleleprob_auto <- alleleprob
    }

    # Extract single phenotype
    pheno_data <- cross2\$pheno[, "${phenotype}", drop = FALSE]

    cat("Running permutation test...\\n")
    cat("Phenotype:", "${phenotype}", "\\n")
    cat("Permutations: 200\\n")
    cat("CPUs: 48 (fork-based parallelism with shared memory)\\n\\n")

    perm_start <- Sys.time()

    # Run permutations with fork-based parallelism (cores=48)
    cat("Starting scan1perm() with 200 permutations using 48 cores (fork-based)...\\n")
    flush.console()

    perm <- scan1perm(alleleprob_auto,
                      pheno_data,
                      kinship = kinship_auto,
                      addcovar = addcovar,
                      Xcovar = Xcovar,
                      cores = 48,
                      n_perm = 200)

    cat("✓ scan1perm() completed successfully\\n")

    perm_time <- as.numeric(difftime(Sys.time(), perm_start, units = "secs"))

    # Save result
    out <- data.frame(perm, check.names = FALSE)
    saveRDS(out, file = "${study_prefix}_${phenotype}_batch_${batch_num}.rds")

    batch_end_time <- Sys.time()
    total_time <- as.numeric(difftime(batch_end_time, batch_start_time, units = "secs"))

    cat("\\n===========================================\\n")
    cat("REPAIR COMPLETE\\n")
    cat("===========================================\\n")
    cat("Phenotype:", "${phenotype}", "\\n")
    cat("Batch:", "${batch_num}", "/ 5\\n")
    cat("Total time:", round(total_time / 60, 2), "minutes\\n")
    cat("Permutation time:", round(perm_time / 60, 2), "minutes\\n")
    cat("End time:", format(batch_end_time), "\\n")
    cat("===========================================\\n")
    """
}

process PERM_AGGREGATE {
    tag "Aggregating permutation results"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '1h'

    input:
    path(perm_results)
    val(study_prefix)
    val(lod_threshold)

    output:
    path("${study_prefix}_fullPerm.rds"), emit: full_perm_matrix
    path("${study_prefix}_permThresh.rds"), emit: perm_thresholds
    path("${study_prefix}_filtered_phenotypes_lod${lod_threshold}.txt"), emit: significant_phenotypes
    path("${study_prefix}_aggregation_log.txt"), emit: aggregation_log

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(utils)
    })

    log_file <- "${study_prefix}_aggregation_log.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message("=== AGGREGATING PERMUTATION RESULTS ===")

    # Get all permutation result files
    perm_files <- list.files(".", pattern = "${study_prefix}_.*_batch_[0-9]+\\\\.rds\$")
    log_message(paste("Found", length(perm_files), "result files"))

    if(length(perm_files) == 0) {
        log_message("ERROR: No result files found")
        stop("No permutation results to aggregate")
    }

    log_message("Combining permutation matrices...")

    # Aggregate results by phenotype
    quiltR <- function(files){
        outList <- list()
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)

        for(i in seq_along(files)){
            tmp <- readRDS(files[i])

            for(j in seq_along(colnames(tmp))){
                col_name <- colnames(tmp)[j]
                if(!col_name %in% names(outList)){
                    outList[[col_name]] <- tmp[[j]]
                } else {
                    outList[[col_name]] <- c(outList[[col_name]], tmp[[j]])
                }
            }
            setTxtProgressBar(pb, i)
        }
        close(pb)

        out <- do.call("cbind", outList)
        return(out)
    }

    permMat <- quiltR(perm_files)

    log_message(paste("Final matrix dimensions:", nrow(permMat), "x", ncol(permMat)))
    log_message(paste("Each phenotype should have 1000 permutations"))

    # Verify we have 1000 perms per phenotype
    if(nrow(permMat) != 1000) {
        log_message(paste("WARNING: Expected 1000 permutations, got", nrow(permMat)))
    }

    log_message("Calculating permutation thresholds...")

    # Calculate thresholds at 4 significance levels
    # 63% = suggestive (1 false positive per genome scan expected)
    # 90% = alpha 0.10
    # 95% = alpha 0.05
    # 99% = alpha 0.01
    permThresh <- apply(permMat, 2, quantile, probs = c(0.63, 0.90, 0.95, 0.99))
    rownames(permThresh) <- c("63%", "90%", "95%", "99%")

    log_message("Threshold summary:")
    log_message(paste("Suggestive 63%: median =", round(median(permThresh[1,]), 2)))
    log_message(paste("Alpha 0.10 90%: median =", round(median(permThresh[2,]), 2)))
    log_message(paste("Alpha 0.05 95%: median =", round(median(permThresh[3,]), 2)))
    log_message(paste("Alpha 0.01 99%: median =", round(median(permThresh[4,]), 2)))

    # Filter phenotypes - use the most stringent threshold
    lod_thresh <- ${lod_threshold}
    # Use the 99% (alpha 0.01) threshold for filtering
    significant_phenos <- colnames(permMat)[permThresh[4,] >= lod_thresh]
    log_message(paste("Phenotypes passing LOD", lod_thresh, "at 99% threshold:", length(significant_phenos)))

    # Save outputs
    saveRDS(permMat, file = "${study_prefix}_fullPerm.rds")
    saveRDS(permThresh, file = "${study_prefix}_permThresh.rds")
    writeLines(significant_phenos, "${study_prefix}_filtered_phenotypes_lod${lod_threshold}.txt")

    log_message("=== AGGREGATION COMPLETE ===")
    log_message(paste("Saved full permutation matrix:", "${study_prefix}_fullPerm.rds"))
    log_message(paste("Saved thresholds:", "${study_prefix}_permThresh.rds"))
    log_message(paste("Saved filtered phenotypes:", "${study_prefix}_filtered_phenotypes_lod${lod_threshold}.txt"))

    close(log_conn)
    """
}

// Main simplified permutation workflow
workflow CHUNKED_PERMUTATION_TESTING {
    take:
        cross2_file
        genoprob_file
        kinship_file
        filtered_phenotypes_file
        study_prefix
        lod_threshold
        run_benchmark  // Ignored in new design
        perm_per_chunk  // Ignored - always 200
        chunks_per_batch  // Ignored - always 1

    main:
        // Check if batch results already exist in published directory
        def batch_results_dir = "${params.outdir}/07_permutation_testing/batches"
        def aggregated_perm_file = "${params.outdir}/07_permutation_testing/${params.study_prefix}_fullPerm.rds"

        // Step 1: Setup - prepare filtered cross2 and phenotype list (always run to get phenotype count)
        PERM_SETUP(
            cross2_file,
            filtered_phenotypes_file,
            study_prefix
        )

        if (file(aggregated_perm_file).exists()) {
            log.info "Aggregated permutation results already exist - skipping batch jobs and aggregation"
            ch_perm_matrix = Channel.fromPath(aggregated_perm_file)
            ch_perm_thresholds = Channel.fromPath("${params.outdir}/07_permutation_testing/${params.study_prefix}_permThresh.rds")
            ch_filtered_phenotypes = Channel.fromPath("${params.outdir}/07_permutation_testing/${params.study_prefix}_filtered_phenotypes_lod${params.lod_threshold}.txt")
            ch_aggregation_log = Channel.fromPath("${params.outdir}/07_permutation_testing/${params.study_prefix}_aggregation_log.txt")

        } else {
            // Check for existing batch results
            def existing_batches = []
            if (file(batch_results_dir).exists()) {
                existing_batches = file("${batch_results_dir}/${params.study_prefix}_*_batch_*.rds").collect { it.name }
            }

            def num_existing = existing_batches.size()

            log.info "=== PERMUTATION BATCH RESUME STATUS ==="
            log.info "Found ${num_existing} existing batch results"

            if (num_existing > 0) {
                log.info "Resuming - will skip ${num_existing} completed batches"
            } else {
                log.info "No existing batch results found - running full permutation workflow"
            }

            // Step 2: Create channel for all batches, then filter out existing ones
            // For each phenotype, create 5 batch jobs (200 perms each)
            def all_batches = PERM_SETUP.out.phenotype_list
                .splitText()
                .map { it.trim() }
                .combine(Channel.from(1..5))  // 5 batches per phenotype

            // Filter to only include batches that don't already exist
            phenotype_batches = all_batches
                .filter { phenotype_name, batch_num ->
                    def batch_file_name = "${params.study_prefix}_${phenotype_name}_batch_${batch_num}.rds"
                    def exists = existing_batches.contains(batch_file_name)
                    return !exists
                }

            // Step 3: Execute permutation batch jobs (only missing ones)
            PERM_BATCH_JOB(
                study_prefix,                     // study_prefix (value channel)
                phenotype_batches.map { it[0] },  // phenotype_name
                phenotype_batches.map { it[1] }   // batch_num
            )

            // Step 4: Validate all batch files after job completion
            log.info "Validating batch results..."

            // Wait for batch jobs to complete before validation
            def validation_trigger = PERM_BATCH_JOB.out.perm_result.collect()
                .ifEmpty(file("${batch_results_dir}"))  // If no new batches, use existing directory

            PERM_VALIDATE_BATCHES(
                PERM_SETUP.out.phenotype_list,
                study_prefix
            )

            // Step 5: Repair corrupt/missing batches if validation failed
            def needs_repair = PERM_VALIDATE_BATCHES.out.status
                .map { status -> status == "FAILED" }
                .first()

            // Create repair channel from missing and corrupt batch lists
            def ch_missing = PERM_VALIDATE_BATCHES.out.missing_batches
                .splitText()
                .filter { !it.startsWith("phenotype") }  // Skip header
                .map { line ->
                    def parts = line.trim().split("\\t")
                    return [parts[0], parts[1].toInteger()]
                }

            def ch_corrupt = PERM_VALIDATE_BATCHES.out.corrupt_batches
                .splitText()
                .filter { !it.startsWith("phenotype") }  // Skip header
                .map { line ->
                    def parts = line.trim().split("\\t")
                    return [parts[0], parts[1].toInteger()]
                }

            def ch_repair_batches = ch_missing.mix(ch_corrupt)

            // Execute repair jobs only if validation failed
            PERM_REPAIR_BATCH(
                study_prefix,
                ch_repair_batches.map { it[0] },  // phenotype
                ch_repair_batches.map { it[1] }   // batch_num
            )

            // Step 6: Combine all batches for aggregation
            // After repairs (if any), collect all batch files from published directory
            def ch_all_batches = PERM_REPAIR_BATCH.out.repaired_batch
                .collect()
                .ifEmpty([])  // Empty if no repairs needed
                .mix(Channel.fromPath("${batch_results_dir}/${params.study_prefix}_*_batch_*.rds").collect())
                .collect()

            PERM_AGGREGATE(
                ch_all_batches,
                study_prefix,
                lod_threshold
            )

            ch_perm_matrix = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenotypes = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log = PERM_AGGREGATE.out.aggregation_log
        }

    emit:
        permutation_matrix = ch_perm_matrix
        permutation_thresholds = ch_perm_thresholds
        filtered_phenotypes = ch_filtered_phenotypes
        filtered_cross2 = PERM_SETUP.out.filtered_cross2
        setup_log = PERM_SETUP.out.setup_log
        aggregation_log = ch_aggregation_log
}
