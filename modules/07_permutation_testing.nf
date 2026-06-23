#!/usr/bin/env nextflow

// Permutation Testing Module — Two-Stage Architecture
//
// Design follows Huda et al. 2020 (paper this pipeline is built from):
//
//   Stage 1 (SCREEN): 50 perms × 1 batch × all filtered phenotypes
//     → conservative threshold per pheno = quantile(perms, 0.90) − bootstrap_SE(q90)
//     → screen-passers = phenos whose observed max LOD (from module 6) exceeds
//       the conservative threshold
//
//   Stage 2 (FULL):   200 perms × 5 batches × screen-passers (= 1000 perms total)
//     → final 4-quantile permThresh.rds (63%, 90%, 95%, 99%)
//
// Each stage runs as an independent 30-day SLURM coordinator job. Coordinators
// survive Nextflow driver death; on resume the workflow detects which stage's
// completion signal is on disk and skips ahead.
//
// Flow:
//   PERM_SETUP
//     → PERM_COORDINATOR_SCREEN
//     → PERM_AGGREGATE_SCREEN
//     → PERM_COORDINATOR_FULL
//     → PERM_AGGREGATE

include {
    PERM_COORDINATOR as PERM_COORDINATOR_SCREEN;
    PERM_COORDINATOR as PERM_COORDINATOR_FULL
} from './07b_perm_coordinator.nf'


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

    log_message("=== PERMUTATION SETUP ===")
    log_message("Loading data...")

    cross2_full <- readRDS("${cross2_file}")
    filtered_phenos <- readLines("${filtered_phenotypes_file}")
    log_message(paste("Original phenotypes:", ncol(cross2_full\$pheno)))
    log_message(paste("Filtered phenotypes:", length(filtered_phenos)))

    pheno_indices <- which(colnames(cross2_full\$pheno) %in% filtered_phenos)
    if (length(pheno_indices) == 0) {
        log_message("ERROR: No phenotypes match between cross2 and filtered list")
        stop("No matching phenotypes found")
    }

    cross2 <- cross2_full
    cross2\$pheno <- cross2_full\$pheno[, pheno_indices, drop = FALSE]
    log_message(paste("Proceeding with", ncol(cross2\$pheno), "filtered phenotypes"))

    saveRDS(cross2, file = "${study_prefix}_filtered_cross2.rds")

    num_phenotypes <- ncol(cross2\$pheno)

    log_message("")
    log_message("=== TWO-STAGE PERMUTATION CONFIG ===")
    log_message(paste("Filtered phenotypes:", num_phenotypes))
    log_message("Stage 1 (SCREEN): 50 perms × 1 batch × all phenos")
    log_message(paste("  Stage 1 batch count:", num_phenotypes))
    log_message("Stage 2 (FULL):   200 perms × 5 batches × screen-passers (= 1000 perms total)")
    log_message("  Stage 2 batch count: depends on stage 1 pass rate (typically 30-50%)")
    log_message("")
    log_message("=== RESOURCE ALLOCATION (per batch SLURM job) ===")
    log_message("CPUs: 48 (fork-based parallelism)")
    log_message("Memory: 200 GB")
    log_message("Wall time: 2 hours per batch")

    writeLines(colnames(cross2\$pheno), "${study_prefix}_phenotype_list.txt")

    log_message("=== SETUP COMPLETE ===")

    close(log_conn)
    """
}


process PERM_AGGREGATE_SCREEN {
    tag "Aggregating stage 1 (screen) permutations"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 4
    memory '128 GB'
    time '1h'

    input:
    path(completion_signal)        // COORDINATOR_screen_COMPLETE trigger
    val(study_prefix)

    output:
    path("${study_prefix}_screen_perms.rds"),                  emit: screen_perms
    path("${study_prefix}_screen_thresholds.rds"),             emit: screen_thresholds
    path("${study_prefix}_screen_passers.txt"),                emit: screen_passers
    path("${study_prefix}_screen_aggregation_log.txt"),        emit: aggregation_log

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(utils))

    log_file <- "${study_prefix}_screen_aggregation_log.txt"
    log_conn <- file(log_file, "w")
    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn); flush(log_conn)
    }

    log_message("=== STAGE 1 (SCREEN) AGGREGATION ===")
    log_message(paste("Trigger file:", "${completion_signal}"))

    # ── Read 50-perm batch files (1 file per pheno) ──────────────────────────
    batch_dir <- file.path("${projectDir}", "${params.outdir}",
                           "07_permutation_testing", "batches_screen")
    log_message(paste("Batch directory:", batch_dir))

    perm_files <- list.files(batch_dir,
                             pattern = paste0("^${study_prefix}_.*_batch_1\\\\.rds\$"),
                             full.names = TRUE)
    log_message(paste("Found", length(perm_files), "screen batch files"))
    if (length(perm_files) == 0) stop("No screen batches to aggregate")

    # Validate completeness against full phenotype list
    phenotype_list_path <- file.path("${projectDir}", "${params.outdir}",
                                     "07_permutation_testing",
                                     paste0("${study_prefix}", "_phenotype_list.txt"))
    expected_phenos <- readLines(phenotype_list_path)
    expected_phenos <- expected_phenos[nchar(trimws(expected_phenos)) > 0]
    if (length(perm_files) < length(expected_phenos)) {
        msg <- paste0("Incomplete screen data: ", length(perm_files), " of ",
                      length(expected_phenos), " phenotype batches present")
        log_message(paste("ERROR:", msg))
        close(log_conn); stop(msg)
    }
    log_message(paste("Completeness check passed:", length(perm_files), "/", length(expected_phenos)))

    # ── Combine into matrix (50 perms × n_phenos) ────────────────────────────
    log_message("Combining 50-perm matrices...")
    perm_list <- lapply(perm_files, readRDS)
    pheno_names <- sapply(perm_list, function(x) colnames(x)[1])
    perm_mat <- do.call(cbind, lapply(perm_list, function(x) x[, 1, drop = FALSE]))
    perm_mat <- as.matrix(perm_mat)
    colnames(perm_mat) <- pheno_names
    log_message(paste("Screen matrix:", nrow(perm_mat), "perms ×", ncol(perm_mat), "phenos"))

    # ── Per-pheno conservative threshold = q90 - bootstrap_SE(q90) ───────────
    # Cox & Hinkley closed-form quantile SE depends on a density estimate at
    # the tail of an n=50 sample (unstable). Bootstrap is more robust here.
    log_message("Computing conservative thresholds (q90 - bootstrap SE) per pheno...")
    set.seed(42)
    B <- 1000
    n <- nrow(perm_mat)
    conservative <- numeric(ncol(perm_mat))
    q90_vec      <- numeric(ncol(perm_mat))
    se_vec       <- numeric(ncol(perm_mat))
    for (j in seq_len(ncol(perm_mat))) {
        x   <- perm_mat[, j]
        q90 <- quantile(x, 0.90, names = FALSE)
        boot_q90 <- replicate(B, quantile(sample(x, n, replace = TRUE), 0.90, names = FALSE))
        se  <- sd(boot_q90)
        q90_vec[j]      <- q90
        se_vec[j]       <- se
        conservative[j] <- q90 - se
    }
    names(conservative) <- colnames(perm_mat)
    names(q90_vec)      <- colnames(perm_mat)
    names(se_vec)       <- colnames(perm_mat)
    log_message(paste("Conservative threshold: median =", round(median(conservative), 3),
                      "  range =", round(min(conservative), 3), "to", round(max(conservative), 3)))
    log_message(paste("q90: median =", round(median(q90_vec), 3),
                      "  bootstrap SE: median =", round(median(se_vec), 3)))

    # ── Compare against module 6 observed LODs ───────────────────────────────
    scan_results_path <- file.path("${projectDir}", "${params.outdir}",
                                    "06_qtl_analysis",
                                    paste0("${study_prefix}", "_scan_results.rds"))
    log_message(paste("Reading observed scan_results from:", scan_results_path))
    if (!file.exists(scan_results_path)) {
        msg <- paste("scan_results.rds not found at", scan_results_path)
        log_message(paste("ERROR:", msg)); close(log_conn); stop(msg)
    }
    scan_results <- readRDS(scan_results_path)
    obs_max_lod  <- apply(scan_results, 2, max, na.rm = TRUE)
    log_message(paste("Observed scan covers", length(obs_max_lod), "phenotypes"))

    # ── Identify screen-passers ──────────────────────────────────────────────
    common_phenos <- intersect(names(conservative), names(obs_max_lod))
    log_message(paste("Phenos with both screen perm + observed LOD:", length(common_phenos)))
    if (length(common_phenos) < length(conservative)) {
        log_message(paste("WARNING:", length(conservative) - length(common_phenos),
                          "screened phenos missing from scan_results — excluded"))
    }

    passers <- common_phenos[obs_max_lod[common_phenos] > conservative[common_phenos]]
    pass_pct <- if (length(common_phenos) > 0) round(100 * length(passers) / length(common_phenos), 1) else 0
    log_message(paste("Screen-passers:", length(passers), "of", length(common_phenos),
                      "(", pass_pct, "% pass rate )"))

    # ── Save artifacts ───────────────────────────────────────────────────────
    saveRDS(perm_mat, "${study_prefix}_screen_perms.rds")

    threshold_df <- data.frame(
        phenotype       = names(conservative),
        q90             = q90_vec[names(conservative)],
        bootstrap_se    = se_vec[names(conservative)],
        conservative_th = conservative,
        obs_max_lod     = obs_max_lod[names(conservative)],
        passed_screen   = names(conservative) %in% passers,
        stringsAsFactors = FALSE
    )
    saveRDS(threshold_df, "${study_prefix}_screen_thresholds.rds")
    writeLines(passers, "${study_prefix}_screen_passers.txt")

    log_message("=== SCREEN AGGREGATION COMPLETE ===")
    log_message(paste("Saved:", "${study_prefix}_screen_perms.rds"))
    log_message(paste("Saved:", "${study_prefix}_screen_thresholds.rds"))
    log_message(paste("Saved:", "${study_prefix}_screen_passers.txt"))

    close(log_conn)
    """
}


process PERM_AGGREGATE {
    tag "Aggregating stage 2 (full) permutations"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '1h'

    input:
    path(completion_signal)
    val(study_prefix)
    val(lod_threshold)

    output:
    path("${study_prefix}_fullPerm.rds"),                              emit: full_perm_matrix
    path("${study_prefix}_permThresh.rds"),                            emit: perm_thresholds
    path("${study_prefix}_filtered_phenotypes_lod${lod_threshold}.txt"), emit: significant_phenotypes
    path("${study_prefix}_aggregation_log.txt"),                       emit: aggregation_log

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(utils))

    log_file <- "${study_prefix}_aggregation_log.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message("=== STAGE 2 (FULL) AGGREGATION ===")
    log_message(paste("Trigger file:", "${completion_signal}"))

    # Read full-stage batch files (5 per pheno × 200 perms)
    batch_dir <- file.path("${projectDir}", "${params.outdir}",
                           "07_permutation_testing", "batches_full")
    log_message(paste("Batch directory:", batch_dir))

    perm_files <- list.files(batch_dir,
                             pattern = paste0("^${study_prefix}_.*_batch_[0-9]+\\\\.rds\$"),
                             full.names = TRUE)
    log_message(paste("Found", length(perm_files), "stage-2 batch files"))

    if (length(perm_files) == 0) {
        log_message("ERROR: No batch result files found in batch directory")
        stop("No permutation results to aggregate")
    }

    # ── Completeness gate against screen_passers.txt ──────────────────────────
    batches_per_pheno   <- 5L
    passers_path <- file.path("${projectDir}", "${params.outdir}",
                              "07_permutation_testing",
                              paste0("${study_prefix}", "_screen_passers.txt"))

    if (file.exists(passers_path)) {
        expected_phenos  <- readLines(passers_path)
        expected_phenos  <- expected_phenos[nchar(trimws(expected_phenos)) > 0]
        n_phenos         <- length(expected_phenos)
        expected_n_files <- n_phenos * batches_per_pheno

        log_message(paste("Screen-passers list:", passers_path))
        log_message(paste("Expected phenotypes (screen-passers):", n_phenos))
        log_message(paste("Expected batch files (", n_phenos, "x", batches_per_pheno, "):", expected_n_files))

        if (length(perm_files) < expected_n_files) {
            missing_n <- expected_n_files - length(perm_files)
            log_message(paste("ERROR:", missing_n, "of", expected_n_files, "batch files are missing."))
            log_message("Permutation jobs are likely still running or some SLURM jobs failed.")
            log_message("Check COORDINATOR_STATUS_full.txt and coordinator_job_logs_full/ for details.")
            close(log_conn)
            stop(paste0("Incomplete permutation data: ", length(perm_files), " of ", expected_n_files,
                        " batch files present. Refusing to aggregate partial results."))
        }

        # Per-phenotype completeness check (vectorized)
        batch_basenames  <- basename(perm_files)
        pheno_from_file  <- gsub("_batch_[0-9]+\\\\.rds\$", "",
                                 gsub(paste0("^${study_prefix}_"), "", batch_basenames))
        batch_counts     <- table(pheno_from_file)
        incomplete       <- names(batch_counts)[batch_counts < batches_per_pheno]
        missing_entirely <- setdiff(expected_phenos, names(batch_counts))
        all_incomplete   <- union(incomplete, missing_entirely)

        if (length(all_incomplete) > 0) {
            log_message(paste("ERROR:", length(all_incomplete), "phenotypes have incomplete batch files:"))
            for (p in head(all_incomplete, 20)) {
                n_found <- if (p %in% names(batch_counts)) batch_counts[[p]] else 0L
                log_message(paste0("  - ", p, " (", n_found, "/", batches_per_pheno, " batches)"))
            }
            if (length(all_incomplete) > 20)
                log_message(paste("  ... and", length(all_incomplete) - 20, "more"))
            close(log_conn)
            stop(paste0(length(all_incomplete), " phenotypes have incomplete permutation batches."))
        }

        log_message(paste("Completeness check passed:", expected_n_files,
                          "batch files present for all", n_phenos, "screen-passers"))
    } else {
        log_message(paste("WARNING: screen_passers list not found at:", passers_path))
        log_message("Cannot validate batch completeness — proceeding anyway")
    }

    log_message("Combining permutation matrices (1000 perms per screen-passer)...")

    quiltR <- function(files) {
        outList <- list()
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
            tmp <- readRDS(files[i])
            for (j in seq_along(colnames(tmp))) {
                col_name <- colnames(tmp)[j]
                if (!col_name %in% names(outList)) {
                    outList[[col_name]] <- tmp[[j]]
                } else {
                    outList[[col_name]] <- c(outList[[col_name]], tmp[[j]])
                }
            }
            setTxtProgressBar(pb, i)
        }
        close(pb)
        do.call("cbind", outList)
    }

    permMat <- quiltR(perm_files)

    log_message(paste("Final matrix dimensions:", nrow(permMat), "x", ncol(permMat)))
    if (nrow(permMat) != 1000) {
        log_message(paste("WARNING: Expected 1000 permutations per phenotype, got", nrow(permMat)))
    }

    log_message("Calculating permutation thresholds...")
    permThresh <- apply(permMat, 2, quantile, probs = c(0.63, 0.90, 0.95, 0.99))
    rownames(permThresh) <- c("63%", "90%", "95%", "99%")

    log_message("Threshold summary:")
    log_message(paste("  Suggestive 63%: median =", round(median(permThresh[1,]), 2)))
    log_message(paste("  Alpha 0.10 90%: median =", round(median(permThresh[2,]), 2)))
    log_message(paste("  Alpha 0.05 95%: median =", round(median(permThresh[3,]), 2)))
    log_message(paste("  Alpha 0.01 99%: median =", round(median(permThresh[4,]), 2)))

    lod_thresh         <- ${lod_threshold}
    significant_phenos <- colnames(permMat)[permThresh[4,] >= lod_thresh]
    log_message(paste("Phenotypes passing LOD", lod_thresh, "at 99% threshold:",
                      length(significant_phenos)))

    saveRDS(permMat,            file = "${study_prefix}_fullPerm.rds")
    saveRDS(permThresh,         file = "${study_prefix}_permThresh.rds")
    writeLines(significant_phenos,
               "${study_prefix}_filtered_phenotypes_lod${lod_threshold}.txt")

    log_message("=== AGGREGATION COMPLETE ===")
    log_message(paste("Saved:", "${study_prefix}_fullPerm.rds"))
    log_message(paste("Saved:", "${study_prefix}_permThresh.rds"))
    log_message(paste("Saved:", "${study_prefix}_filtered_phenotypes_lod${lod_threshold}.txt"))

    close(log_conn)
    """
}


// ── Main permutation workflow ────────────────────────────────────────────────
//
// Resume-state detection: each stage publishes its completion signal to disk
// (COORDINATOR_screen_COMPLETE, screen_passers.txt, COORDINATOR_full_COMPLETE,
// fullPerm.rds). The workflow checks these in priority order and skips ahead.
//
//   Case A — fullPerm.rds exists:               skip everything
//   Case B — COORDINATOR_full_COMPLETE exists:  run final PERM_AGGREGATE only
//   Case C — screen_passers.txt exists:         run PERM_COORDINATOR_FULL → PERM_AGGREGATE
//   Case D — COORDINATOR_screen_COMPLETE exists: run PERM_AGGREGATE_SCREEN onward
//   Case E — fresh:                             run everything from PERM_COORDINATOR_SCREEN
//
workflow CHUNKED_PERMUTATION_TESTING {
    take:
        cross2_file
        genoprob_file            // Not used in coordinator design (kept for API compat)
        kinship_file             // Not used in coordinator design (kept for API compat)
        filtered_phenotypes_file
        study_prefix
        lod_threshold
        run_benchmark            // Not used (kept for API compat)
        perm_per_chunk           // Not used (always 200)
        chunks_per_batch         // Not used (always 5 stage 2 / 1 stage 1)
        interactive_covar        // "null" or covariate name; matches scan1 in module 6

    main:
        def perm_dir          = "${params.outdir}/07_permutation_testing"
        def aggregated_perm   = "${perm_dir}/${params.study_prefix}_fullPerm.rds"
        def screen_done_path  = "${perm_dir}/COORDINATOR_screen_COMPLETE"
        def screen_passers_p  = "${perm_dir}/${params.study_prefix}_screen_passers.txt"
        def full_done_path    = "${perm_dir}/COORDINATOR_full_COMPLETE"

        // Step 0: Setup — always runs (idempotent, fast)
        PERM_SETUP(
            cross2_file,
            filtered_phenotypes_file,
            study_prefix
        )

        if (file(aggregated_perm).exists()) {
            // ── Case A: fully done ───────────────────────────────────────
            log.info "fullPerm.rds exists — skipping all permutation steps"
            ch_perm_matrix     = Channel.fromPath(aggregated_perm)
            ch_perm_thresholds = Channel.fromPath("${perm_dir}/${params.study_prefix}_permThresh.rds")
            ch_filtered_phenos = Channel.fromPath("${perm_dir}/${params.study_prefix}_filtered_phenotypes_lod${params.lod_threshold}.txt")
            ch_aggregation_log = Channel.fromPath("${perm_dir}/${params.study_prefix}_aggregation_log.txt")

        } else if (file(full_done_path).exists()) {
            // ── Case B: stage 2 done, only need final aggregation ────────
            log.info "COORDINATOR_full_COMPLETE found — running PERM_AGGREGATE only"
            PERM_AGGREGATE(
                Channel.fromPath(full_done_path),
                study_prefix,
                lod_threshold
            )
            ch_perm_matrix     = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenos = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log = PERM_AGGREGATE.out.aggregation_log

        } else if (file(screen_passers_p).exists()) {
            // ── Case C: screen passers exist → resume stage 2 ────────────
            log.info "screen_passers.txt found — submitting/resuming PERM_COORDINATOR_FULL"
            PERM_COORDINATOR_FULL(
                Channel.fromPath(screen_passers_p),
                study_prefix,
                interactive_covar,
                Channel.value("full"),
                Channel.value(200),
                Channel.value(5)
            )
            PERM_AGGREGATE(
                PERM_COORDINATOR_FULL.out.completion_signal,
                study_prefix,
                lod_threshold
            )
            ch_perm_matrix     = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenos = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log = PERM_AGGREGATE.out.aggregation_log

        } else if (file(screen_done_path).exists()) {
            // ── Case D: screen done, need aggregation + stage 2 ──────────
            log.info "COORDINATOR_screen_COMPLETE found — running PERM_AGGREGATE_SCREEN onward"
            PERM_AGGREGATE_SCREEN(
                Channel.fromPath(screen_done_path),
                study_prefix
            )
            PERM_COORDINATOR_FULL(
                PERM_AGGREGATE_SCREEN.out.screen_passers,
                study_prefix,
                interactive_covar,
                Channel.value("full"),
                Channel.value(200),
                Channel.value(5)
            )
            PERM_AGGREGATE(
                PERM_COORDINATOR_FULL.out.completion_signal,
                study_prefix,
                lod_threshold
            )
            ch_perm_matrix     = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenos = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log = PERM_AGGREGATE.out.aggregation_log

        } else {
            // ── Case E: fresh run — full two-stage flow ──────────────────
            log.info "Fresh run — submitting PERM_COORDINATOR_SCREEN (stage 1)"
            PERM_COORDINATOR_SCREEN(
                PERM_SETUP.out.phenotype_list,
                study_prefix,
                interactive_covar,
                Channel.value("screen"),
                Channel.value(50),
                Channel.value(1)
            )
            PERM_AGGREGATE_SCREEN(
                PERM_COORDINATOR_SCREEN.out.completion_signal,
                study_prefix
            )
            PERM_COORDINATOR_FULL(
                PERM_AGGREGATE_SCREEN.out.screen_passers,
                study_prefix,
                interactive_covar,
                Channel.value("full"),
                Channel.value(200),
                Channel.value(5)
            )
            PERM_AGGREGATE(
                PERM_COORDINATOR_FULL.out.completion_signal,
                study_prefix,
                lod_threshold
            )
            ch_perm_matrix     = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenos = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log = PERM_AGGREGATE.out.aggregation_log
        }

    emit:
        permutation_matrix     = ch_perm_matrix
        permutation_thresholds = ch_perm_thresholds
        filtered_phenotypes    = ch_filtered_phenos
        filtered_cross2        = PERM_SETUP.out.filtered_cross2
        setup_log              = PERM_SETUP.out.setup_log
        aggregation_log        = ch_aggregation_log
}
