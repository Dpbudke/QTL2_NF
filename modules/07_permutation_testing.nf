process PERMUTATION_SETUP {
    tag "Chunking ${n_pheno} phenotypes for permutation testing"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy', pattern: "*.txt"

    cpus 4
    memory '32 GB'
    time '1h'

    input:
    path(filtered_phenotypes_file)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_perm_chunks.txt"), emit: chunk_file
    path("${prefix}_perm_batch_assignments.txt"), emit: batch_file
    path("${prefix}_perm_chunk_summary.txt"), emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # Load filtered phenotypes list
    filtered_phenotypes <- readLines("${filtered_phenotypes_file}")
    n_pheno <- length(filtered_phenotypes)

    chunk_size <- 200  # Optimized chunk size for permutation testing
    batch_size <- 6    # Process chunks in batches of 6 for 8 batches total (optimized for 13TB memory)
    n_chunks <- ceiling(n_pheno / chunk_size)
    n_batches <- ceiling(n_chunks / batch_size)

    # Create chunk assignments for filtered phenotypes
    chunk_assignments <- data.frame(
        chunk_id = rep(1:n_chunks, each = chunk_size)[1:n_pheno],
        phenotype = filtered_phenotypes
    )

    write.table(chunk_assignments, "${prefix}_perm_chunks.txt",
                quote = FALSE, row.names = FALSE, sep = "\\t")

    # Create batch assignments
    batch_assignments <- data.frame(
        batch_id = rep(1:n_batches, each = batch_size)[1:n_chunks],
        chunk_id = 1:n_chunks
    )

    write.table(batch_assignments, "${prefix}_perm_batch_assignments.txt",
                quote = FALSE, row.names = FALSE, sep = "\\t")

    summary_text <- c(
        "=== Permutation Testing Chunking Summary ===",
        paste("Total filtered phenotypes:", n_pheno),
        paste("Chunk size:", chunk_size, "phenotypes per chunk"),
        paste("Number of chunks:", n_chunks),
        paste("Batch size:", batch_size, "chunks per batch"),
        paste("Number of batches:", n_batches),
        paste("Resources per batch: 48 CPUs, 1.6TB memory"),
        paste("Target batches: 8 (optimized for 13TB total memory)"),
        paste("Permutations per phenotype: 1000"),
        paste("Batched submission for optimal HPC utilization"),
        paste("LOD threshold applied:", ${lod_threshold})
    )
    writeLines(summary_text, "${prefix}_perm_chunk_summary.txt")
    """
}

process PERMUTATION_BATCH {
    tag "Permutation Batch ${batch_id} for ${prefix}"
    publishDir "${params.outdir}/07_permutation_testing/batches", mode: 'copy'

    cpus 48
    memory '1.6 TB'
    time '12h'  // Time for permutation batch processing

    input:
    each batch_id
    path(cross2_file)
    path(genoprob_file)
    path(kinship_file)
    path(chunk_file)
    path(batch_file)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_perm_batch${batch_id}_results.rds"), emit: batch_perm_results
    path("${prefix}_perm_batch${batch_id}_thresholds.csv"), emit: batch_thresholds
    path("${prefix}_perm_batch${batch_id}_log.txt"), emit: batch_log

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(qtl2))

    # Initialize logging
    log_file <- "${prefix}_perm_batch${batch_id}_log.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n")
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message(paste("=== PERMUTATION BATCH ${batch_id} STARTED ==="))
    log_message("Study: ${prefix}")
    log_message("LOD threshold: ${lod_threshold}")
    log_message(paste("Node:", Sys.info()["nodename"]))
    log_message(paste("SLURM Job ID:", Sys.getenv("SLURM_JOB_ID", "not_set")))

    # Load input data
    log_message("Loading input files...")
    cross2 <- readRDS("${cross2_file}")
    genoprob <- readRDS("${genoprob_file}")
    kinship_loco <- readRDS("${kinship_file}")

    # Load assignments
    chunk_assignments <- read.table("${chunk_file}", header = TRUE, sep = "\\t")
    batch_assignments <- read.table("${batch_file}", header = TRUE, sep = "\\t")

    # Get chunks for this specific batch
    chunks_in_batch <- batch_assignments\$chunk_id[batch_assignments\$batch_id == ${batch_id}]
    log_message(paste("Chunks in this batch:", paste(chunks_in_batch, collapse = ", ")))
    log_message(paste("Total chunks in batch:", length(chunks_in_batch)))

    # Initialize result storage for this batch
    batch_perm_results <- list()
    batch_thresholds <- data.frame()

    batch_start_time <- Sys.time()

    # Process chunks in this batch sequentially
    for (chunk_id in chunks_in_batch) {
        log_message(paste("Processing permutation chunk", chunk_id))

        # Get phenotypes for this chunk
        phenos_in_chunk <- chunk_assignments\$phenotype[chunk_assignments\$chunk_id == chunk_id]
        pheno_subset <- cross2\$pheno[, phenos_in_chunk, drop = FALSE]

        chunk_start_time <- Sys.time()

        # Run permutation testing for this chunk
        tryCatch({
            # Check if we have X chromosome and get chromosome lengths
            has_x_chr <- "X" %in% names(cross2\$gmap)
            chr_lengths <- NULL
            if (has_x_chr) {
                # Extract chromosome lengths from genetic map
                chr_lengths <- sapply(cross2\$gmap, function(x) max(x, na.rm=TRUE))
            }

            # Run 1000 permutations using all available cores (48)
            perm_result <- scan1perm(genoprob, pheno_subset, kinship_loco,
                                   n_perm=1000, cores=48, perm_Xsp=has_x_chr,
                                   chr_lengths=chr_lengths, perm_strata=NULL,
                                   reindex=TRUE)

            chunk_end_time <- Sys.time()
            chunk_duration <- round(as.numeric(chunk_end_time - chunk_start_time, units = "mins"), 2)

            # Calculate significance thresholds for this chunk
            chunk_thresholds <- data.frame()
            for (i in 1:ncol(perm_result)) {
                phenotype_name <- colnames(perm_result)[i]
                perm_lods <- perm_result[, i]

                # Calculate thresholds at different significance levels
                threshold_95 <- quantile(perm_lods, 0.95, na.rm = TRUE)
                threshold_99 <- quantile(perm_lods, 0.99, na.rm = TRUE)
                threshold_999 <- quantile(perm_lods, 0.999, na.rm = TRUE)

                chunk_thresholds <- rbind(chunk_thresholds, data.frame(
                    phenotype = phenotype_name,
                    threshold_95 = threshold_95,
                    threshold_99 = threshold_99,
                    threshold_999 = threshold_999,
                    stringsAsFactors = FALSE
                ))
            }

            # Store results
            batch_perm_results[[as.character(chunk_id)]] <- perm_result
            batch_thresholds <- rbind(batch_thresholds, chunk_thresholds)

            log_message(paste("Chunk", chunk_id, "completed in", chunk_duration, "minutes,",
                            ncol(perm_result), "phenotypes processed"))

        }, error = function(e) {
            log_message(paste("ERROR in permutation chunk", chunk_id, ":", conditionMessage(e)))
            stop(paste("Permutation chunk", chunk_id, "failed"))
        })
    }

    # Combine batch permutation results
    batch_combined <- do.call(cbind, batch_perm_results)

    # Save batch results
    saveRDS(batch_combined, "${prefix}_perm_batch${batch_id}_results.rds")
    write.csv(batch_thresholds, "${prefix}_perm_batch${batch_id}_thresholds.csv", row.names = FALSE)

    batch_end_time <- Sys.time()
    batch_duration <- round(as.numeric(batch_end_time - batch_start_time, units = "mins"), 2)

    log_message(paste("=== PERMUTATION BATCH ${batch_id} COMPLETED ==="))
    log_message(paste("Batch duration:", batch_duration, "minutes"))
    log_message(paste("Chunks processed:", length(chunks_in_batch)))
    log_message(paste("Phenotypes processed:", ncol(batch_combined)))
    log_message(paste("Permutations completed:", nrow(batch_combined) * ncol(batch_combined)))

    close(log_conn)
    """
}

process COMBINE_PERMUTATION_RESULTS {
    tag "Combining ${prefix} permutation results"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 16
    memory '256 GB'
    time '2h'

    input:
    path(batch_perm_results)
    path(batch_thresholds)
    path(batch_logs)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_permutation_results.rds"), emit: perm_results
    path("${prefix}_significance_thresholds.csv"), emit: thresholds
    path("permutation_validation_report.txt"), emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(qtl2))

    cat("=== COMBINING PERMUTATION RESULTS ===\\n")
    cat("Study: ${prefix}\\n\\n")

    # Load and combine batch permutation results
    perm_files <- list.files(pattern = ".*_perm_batch[0-9]+_results\\\\.rds")
    perm_files <- perm_files[order(as.numeric(gsub(".*batch([0-9]+).*", "\\\\1", perm_files)))]

    cat("Loading", length(perm_files), "permutation batch results...\\n")
    perm_list <- lapply(perm_files, readRDS)
    combined_perms <- do.call(cbind, perm_list)

    saveRDS(combined_perms, "${prefix}_permutation_results.rds")
    cat("Combined permutation results saved\\n\\n")

    # Combine threshold files
    threshold_files <- list.files(pattern = ".*_perm_batch[0-9]+_thresholds\\\\.csv")
    all_thresholds <- do.call(rbind, lapply(threshold_files, read.csv))
    write.csv(all_thresholds, "${prefix}_significance_thresholds.csv", row.names = FALSE)

    # Create validation report
    validation_log <- c(
        "=== Permutation Testing Validation Report ===",
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Study: ${prefix}"),
        paste("LOD threshold applied: ${lod_threshold}"),
        "",
        "=== Parallel Batch Processing Strategy ===",
        paste("Batches processed in parallel:", length(perm_files)),
        paste("Chunks per batch: ~6"),
        paste("Resources per batch: 48 CPUs, 1.6TB memory"),
        paste("Permutations per phenotype: 1000"),
        paste("Processing approach: Parallel batches, sequential chunks within batch"),
        "",
        "=== Results Summary ===",
        paste("Total phenotypes tested:", ncol(combined_perms)),
        paste("Total permutations completed:", nrow(combined_perms) * ncol(combined_perms)),
        paste("Average 95% threshold:", round(mean(all_thresholds\$threshold_95, na.rm=TRUE), 2)),
        paste("Average 99% threshold:", round(mean(all_thresholds\$threshold_99, na.rm=TRUE), 2)),
        paste("Range of 95% thresholds:", paste(range(all_thresholds\$threshold_95, na.rm=TRUE), collapse=" - ")),
        "",
        "=== Threshold Distribution ===",
        paste("95% significance level - Min:", round(min(all_thresholds\$threshold_95, na.rm=TRUE), 2)),
        paste("95% significance level - Max:", round(max(all_thresholds\$threshold_95, na.rm=TRUE), 2)),
        paste("99% significance level - Min:", round(min(all_thresholds\$threshold_99, na.rm=TRUE), 2)),
        paste("99% significance level - Max:", round(max(all_thresholds\$threshold_99, na.rm=TRUE), 2))
    )

    writeLines(validation_log, "permutation_validation_report.txt")

    cat("=== PERMUTATION COMBINATION COMPLETED ===\\n")
    cat("Total phenotypes tested:", ncol(combined_perms), "\\n")
    cat("Total permutations:", nrow(combined_perms) * ncol(combined_perms), "\\n")
    cat("Average 95% threshold:", round(mean(all_thresholds\$threshold_95, na.rm=TRUE), 2), "\\n")
    """
}