process GENOME_SCAN_SETUP {
    tag "Chunking ${n_pheno} phenotypes for ${prefix}"
    publishDir "${params.outdir}/06_qtl_analysis", mode: 'copy', pattern: "*.txt"

    cpus 4
    memory '32 GB'
    time '1h'

    input:
    path(cross2_file)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_pheno_chunks.txt"), emit: chunk_file
    path("${prefix}_batch_assignments.txt"), emit: batch_file
    path("${prefix}_chunk_summary.txt"), emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(qtl2))

    cross2 <- readRDS("${cross2_file}")
    n_pheno <- ncol(cross2\$pheno)

    chunk_size <- 200  # Optimized chunk size for 1TB memory
    batch_size <- 10   # Process chunks in batches of 10 to avoid infrastructure issues
    n_chunks <- ceiling(n_pheno / chunk_size)
    n_batches <- ceiling(n_chunks / batch_size)

    # Create chunk assignments
    chunk_assignments <- data.frame(
        chunk_id = rep(1:n_chunks, each = chunk_size)[1:n_pheno],
        phenotype = colnames(cross2\$pheno)
    )

    write.table(chunk_assignments, "${prefix}_pheno_chunks.txt",
                quote = FALSE, row.names = FALSE, sep = "\\t")

    # Create batch assignments
    batch_assignments <- data.frame(
        batch_id = rep(1:n_batches, each = batch_size)[1:n_chunks],
        chunk_id = 1:n_chunks
    )

    write.table(batch_assignments, "${prefix}_batch_assignments.txt",
                quote = FALSE, row.names = FALSE, sep = "\\t")

    summary_text <- c(
        "=== Phenotype Chunking Summary ===",
        paste("Total phenotypes:", n_pheno),
        paste("Chunk size:", chunk_size, "phenotypes per chunk"),
        paste("Number of chunks:", n_chunks),
        paste("Batch size:", batch_size, "chunks per batch"),
        paste("Number of batches:", n_batches),
        paste("Resources per batch: 48 CPUs, 1.6TB memory"),
        paste("Batched submission to avoid infrastructure issues"),
        paste("LOD threshold for filtering:", ${lod_threshold})
    )
    writeLines(summary_text, "${prefix}_chunk_summary.txt")
    """
}

process GENOME_SCAN_BATCH {
    tag "Batch ${batch_id} for ${prefix}"
    publishDir "${params.outdir}/06_qtl_analysis/batches", mode: 'copy'

    cpus 48
    memory '1.6 TB'
    time '12h'  // Time for one batch (10 chunks)

    input:
    each batch_id
    path(cross2_file)
    path(genoprob_file)
    path(alleleprob_file)
    path(kinship_file)
    path(chunk_file)
    path(batch_file)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_batch${batch_id}_results.rds"), emit: batch_results
    path("${prefix}_batch${batch_id}_peaks.csv"), emit: batch_peaks
    path("${prefix}_batch${batch_id}_filtered_phenos.txt"), emit: batch_filtered_phenos
    path("${prefix}_batch${batch_id}_log.txt"), emit: batch_log

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(qtl2))

    # Initialize logging
    log_file <- "${prefix}_batch${batch_id}_log.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n")
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message(paste("=== BATCH ${batch_id} STARTED ==="))
    log_message("Study: ${prefix}")
    log_message("LOD threshold: ${lod_threshold}")
    log_message(paste("Node:", Sys.info()["nodename"]))
    log_message(paste("SLURM Job ID:", Sys.getenv("SLURM_JOB_ID", "not_set")))

    # Load input data
    log_message("Loading input files...")
    cross2 <- readRDS("${cross2_file}")
    genoprob <- readRDS("${genoprob_file}")
    alleleprob <- readRDS("${alleleprob_file}")
    kinship <- readRDS("${kinship_file}")

    # Load assignments
    chunk_assignments <- read.table("${chunk_file}", header = TRUE, sep = "\\t")
    batch_assignments <- read.table("${batch_file}", header = TRUE, sep = "\\t")

    # Function to find chromosome and position for a marker
    find_marker_position <- function(marker_id, cross2_obj) {
        for (chr in names(cross2_obj\$gmap)) {
            if (marker_id %in% names(cross2_obj\$gmap[[chr]])) {
                pos <- cross2_obj\$gmap[[chr]][marker_id]
                return(list(chr = chr, pos = as.numeric(pos)))
            }
        }
        return(list(chr = NA, pos = NA))
    }

    # Get chunks for this specific batch
    chunks_in_batch <- batch_assignments\$chunk_id[batch_assignments\$batch_id == ${batch_id}]
    log_message(paste("Chunks in this batch:", paste(chunks_in_batch, collapse = ", ")))
    log_message(paste("Total chunks in batch:", length(chunks_in_batch)))

    # Initialize result storage for this batch
    batch_results <- list()
    batch_peaks <- data.frame()
    batch_filtered_phenos <- character()

    batch_start_time <- Sys.time()

    # Process chunks in this batch sequentially
    for (chunk_id in chunks_in_batch) {
        log_message(paste("Processing chunk", chunk_id))

        # Get phenotypes for this chunk
        phenos_in_chunk <- chunk_assignments\$phenotype[chunk_assignments\$chunk_id == chunk_id]
        pheno_subset <- cross2\$pheno[, phenos_in_chunk, drop = FALSE]

        chunk_start_time <- Sys.time()

        # Run scan for this chunk
        tryCatch({
            scan_result <- scan1(
                genoprob = alleleprob,
                pheno = pheno_subset,
                kinship = kinship,
                cores = 0  # Use all available cores
            )

            chunk_end_time <- Sys.time()
            chunk_duration <- round(as.numeric(chunk_end_time - chunk_start_time, units = "mins"), 2)

            # Find peaks with proper position extraction
            chunk_peaks <- data.frame()
            for (i in 1:ncol(scan_result)) {
                max_lod <- max(scan_result[, i], na.rm = TRUE)
                if (max_lod >= ${lod_threshold}) {
                    max_idx <- which.max(scan_result[, i])
                    marker_id <- rownames(scan_result)[max_idx]
                    pos_info <- find_marker_position(marker_id, cross2)
                    chunk_peaks <- rbind(chunk_peaks, data.frame(
                        chr = pos_info\$chr,
                        phenotype = colnames(scan_result)[i],
                        pos = pos_info\$pos,
                        lod = max_lod,
                        stringsAsFactors = FALSE
                    ))
                }
            }

            # Store results
            batch_results[[as.character(chunk_id)]] <- scan_result
            batch_peaks <- rbind(batch_peaks, chunk_peaks)

            phenos_with_peaks <- unique(chunk_peaks\$phenotype)
            batch_filtered_phenos <- c(batch_filtered_phenos, phenos_with_peaks)

            log_message(paste("Chunk", chunk_id, "completed in", chunk_duration, "minutes,",
                            nrow(chunk_peaks), "peaks found"))

        }, error = function(e) {
            log_message(paste("ERROR in chunk", chunk_id, ":", conditionMessage(e)))
            stop(paste("Chunk", chunk_id, "failed"))
        })
    }

    # Combine batch results
    batch_combined <- do.call(cbind, batch_results)

    # Save batch results
    saveRDS(batch_combined, "${prefix}_batch${batch_id}_results.rds")
    write.csv(batch_peaks, "${prefix}_batch${batch_id}_peaks.csv", row.names = FALSE)

    batch_filtered_phenos <- unique(batch_filtered_phenos)
    writeLines(batch_filtered_phenos, "${prefix}_batch${batch_id}_filtered_phenos.txt")

    batch_end_time <- Sys.time()
    batch_duration <- round(as.numeric(batch_end_time - batch_start_time, units = "mins"), 2)

    log_message(paste("=== BATCH ${batch_id} COMPLETED ==="))
    log_message(paste("Batch duration:", batch_duration, "minutes"))
    log_message(paste("Chunks processed:", length(chunks_in_batch)))
    log_message(paste("Peaks found:", nrow(batch_peaks)))
    log_message(paste("Phenotypes with peaks:", length(batch_filtered_phenos)))

    close(log_conn)
    """
}

process COMBINE_BATCH_RESULTS {
    tag "Combining ${prefix} batch results"
    publishDir "${params.outdir}/06_qtl_analysis", mode: 'copy'

    cpus 16
    memory '256 GB'
    time '2h'

    input:
    path(batch_results)
    path(batch_peaks)
    path(batch_filtered_phenos)
    path(batch_logs)
    val(prefix)
    val(lod_threshold)

    output:
    path("${prefix}_scan_results.rds"), emit: scan_results
    path("${prefix}_all_peaks.csv"), emit: peaks
    path("${prefix}_filtered_phenotypes.txt"), emit: filtered_phenotypes
    path("${prefix}_batch_processing_summary.txt"), emit: batch_summary

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(qtl2))

    cat("=== COMBINING BATCH RESULTS ===\\n")
    cat("Study: ${prefix}\\n\\n")

    # Load and combine batch results
    batch_files <- list.files(pattern = ".*_batch[0-9]+_results\\\\.rds")
    batch_files <- batch_files[order(as.numeric(gsub(".*batch([0-9]+).*", "\\\\1", batch_files)))]

    cat("Loading", length(batch_files), "batch results...\\n")
    scan_list <- lapply(batch_files, readRDS)
    combined_scan <- do.call(cbind, scan_list)

    saveRDS(combined_scan, "${prefix}_scan_results.rds")
    cat("Combined scan results saved\\n\\n")

    # Combine peak files
    peak_files <- list.files(pattern = ".*_batch[0-9]+_peaks\\\\.csv")
    all_peaks <- do.call(rbind, lapply(peak_files, read.csv))
    write.csv(all_peaks, "${prefix}_all_peaks.csv", row.names = FALSE)

    # Combine filtered phenotype files
    pheno_files <- list.files(pattern = ".*_batch[0-9]+_filtered_phenos\\\\.txt")
    all_phenos <- unique(unlist(lapply(pheno_files, readLines)))
    writeLines(all_phenos, "${prefix}_filtered_phenotypes.txt")

    # Create summary report
    report <- c(
        "=== Parallel Batch Processing Summary ===",
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Study: ${prefix}"),
        paste("LOD threshold: ${lod_threshold}"),
        "",
        "=== Processing Strategy ===",
        paste("Batches processed in parallel:", length(batch_files)),
        paste("Chunks per batch: ~10"),
        paste("Resources per batch: 48 CPUs, 1TB memory"),
        paste("Processing approach: Parallel batches, sequential chunks within batch"),
        "",
        "=== Results ===",
        paste("Total phenotypes scanned:", ncol(combined_scan)),
        paste("Phenotypes with peaks >", ${lod_threshold}, "LOD:", length(all_phenos)),
        paste("Total peaks found:", nrow(all_peaks)),
        "",
        "=== Per-Chromosome Summary ===",
        paste(capture.output(table(all_peaks\$chr)), collapse = "\\n")
    )

    writeLines(report, "${prefix}_batch_processing_summary.txt")

    cat("=== BATCH COMBINATION COMPLETED ===\\n")
    cat("Total peaks:", nrow(all_peaks), "\\n")
    cat("Phenotypes for permutation:", length(all_phenos), "\\n")
    """
}

