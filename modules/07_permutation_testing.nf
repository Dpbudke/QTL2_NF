#!/usr/bin/env nextflow

// Permutation Testing Module - Coordinator Architecture
// Design: Nextflow submits ONE long-lived SLURM coordinator job (7-day wall time).
// The coordinator manages all batch jobs internally via sbatch/sacct, handling
// submission, monitoring, and retries — completely independent of the Nextflow
// driver session. This prevents driver-death (e.g., VS_Code wall-time) from
// interrupting a multi-day permutation run.
//
// Flow: PERM_SETUP → PERM_COORDINATOR → PERM_AGGREGATE

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

    num_phenotypes    <- ncol(cross2\$pheno)
    perms_per_batch   <- 200
    batches_per_pheno <- 5
    total_batches     <- num_phenotypes * batches_per_pheno

    log_message("")
    log_message("=== BATCH CONFIGURATION ===")
    log_message(paste("Phenotypes to process:", num_phenotypes))
    log_message(paste("Permutations per batch:", perms_per_batch))
    log_message(paste("Batches per phenotype:", batches_per_pheno))
    log_message(paste("Total batches:", total_batches))
    log_message(paste("Total permutations:", num_phenotypes * 1000))
    log_message("")
    log_message("=== RESOURCE ALLOCATION (per batch SLURM job) ===")
    log_message("CPUs: 48 (fork-based parallelism)")
    log_message("Memory: 200 GB")
    log_message("Wall time: 2 hours per batch")
    log_message(paste("Estimated total time: ~",
                      round(total_batches * 7 / 60 / 35, 1),
                      "hours (with 35 parallel batches)"))
    log_message("")
    log_message("=== COORDINATOR ARCHITECTURE ===")
    log_message("One 7-day SLURM coordinator job manages all batch submissions.")
    log_message("Independent of Nextflow driver session - survives any driver interruption.")

    writeLines(colnames(cross2\$pheno), "${study_prefix}_phenotype_list.txt")

    log_message("=== SETUP COMPLETE ===")
    log_message(paste("Ready to launch", total_batches, "batch jobs via coordinator"))

    close(log_conn)
    """
}


process PERM_COORDINATOR {
    tag "Coordinating all permutation batches (7-day SLURM job)"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    // No container — this process needs host SLURM tools (sbatch, sacct).
    // Batch worker jobs submitted from here use the container explicitly.
    container null

    cpus   1
    memory '4 GB'
    time   '7d'
    errorStrategy 'retry'
    maxRetries 0   // Internal retry logic handles batch failures

    input:
    path(phenotype_list)
    val(study_prefix)

    output:
    path("${study_prefix}_coordinator_log.txt"), emit: coordinator_log
    path("COORDINATOR_COMPLETE"),               emit: completion_signal

    script:
    // Nextflow resolves: ${study_prefix}, ${projectDir}, ${params.outdir}
    // Bash variables use \${...} to survive Nextflow template substitution
    """
    #!/bin/bash
    set -euo pipefail

    # ── Paths ────────────────────────────────────────────────────────────────
    STUDY_PREFIX="${study_prefix}"
    PROJECT_DIR="${projectDir}"
    OUTDIR="${params.outdir}"
    SIF="\${PROJECT_DIR}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"
    PERM_DIR="\${PROJECT_DIR}/\${OUTDIR}/07_permutation_testing"
    BATCH_DIR="\${PERM_DIR}/batches"
    JOBS_LOG_DIR="\${PERM_DIR}/coordinator_job_logs"
    WORKER_SCRIPT="\${PERM_DIR}/perm_batch_worker.R"
    LOG="${study_prefix}_coordinator_log.txt"

    mkdir -p "\${BATCH_DIR}" "\${JOBS_LOG_DIR}"

    log() { echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$*" | tee -a "\${LOG}"; }

    log "========================================================"
    log "PERMUTATION COORDINATOR STARTED"
    log "Study:       \${STUDY_PREFIX}"
    log "Project:     \${PROJECT_DIR}"
    log "Output dir:  \${OUTDIR}"
    log "Batch dir:   \${BATCH_DIR}"
    log "Container:   \${SIF}"
    log "========================================================"

    # ── Sanity checks ────────────────────────────────────────────────────────
    if ! command -v sbatch &>/dev/null; then
        log "ERROR: sbatch not found. Nested SLURM submission requires SLURM"
        log "client tools to be available on the compute node."
        exit 1
    fi
    if [[ ! -f "\${SIF}" ]]; then
        log "ERROR: Singularity image not found: \${SIF}"
        exit 1
    fi

    # ── Write shared R worker script ─────────────────────────────────────────
    # Written to PERM_DIR (shared filesystem) so all batch compute nodes can
    # access it. Uses command-line args so no per-phenotype file needed.
    cat > "\${WORKER_SCRIPT}" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: perm_batch_worker.R <pheno> <batch_num> <study_prefix> <project_dir> <outdir>")
}
pheno      <- args[1]
batch_num  <- as.integer(args[2])
study_pfx  <- args[3]
proj_dir   <- args[4]
outdir     <- args[5]

suppressPackageStartupMessages(library(qtl2))

batch_start <- Sys.time()
cat("===========================================\\n")
cat("BATCH:", batch_num, "/ 5\\n")
cat("PHENOTYPE:", pheno, "\\n")
cat("PERMUTATIONS: 200\\n")
cat("START TIME:", format(batch_start), "\\n")
cat("===========================================\\n\\n")

perm_dir  <- file.path(proj_dir, outdir, "07_permutation_testing")
batch_dir <- file.path(perm_dir,  "batches")
scan_dir  <- file.path(proj_dir, outdir, "05_genome_scan_preparation")
out_file  <- file.path(batch_dir,
                        paste0(study_pfx, "_", pheno, "_batch_", batch_num, ".rds"))

# Skip if already done (handles duplicate submissions on restart)
if (file.exists(out_file) && file.size(out_file) > 1000) {
    cat("Output already exists, skipping:", out_file, "\\n")
    quit(status = 0)
}

cat("Loading data files...\\n")
load_start   <- Sys.time()
cross2       <- readRDS(file.path(perm_dir, paste0(study_pfx, "_filtered_cross2.rds")))
alleleprob   <- readRDS(file.path(scan_dir,  paste0(study_pfx, "_alleleprob.rds")))
kinship_loco <- readRDS(file.path(scan_dir,  paste0(study_pfx, "_kinship_loco.rds")))
load_time    <- as.numeric(difftime(Sys.time(), load_start, units = "secs"))
cat("Data loaded in", round(load_time, 1), "seconds\\n\\n")

# Covariates
if (!is.null(cross2\$covar)) {
    covar_data <- cross2\$covar
    if ("coat_color" %in% colnames(covar_data)) {
        covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
    }
    covar_formula <- paste("~", paste(colnames(covar_data), collapse = " + "))
    addcovar <- model.matrix(as.formula(covar_formula), data = covar_data)[, -1, drop = FALSE]
} else {
    addcovar <- NULL
}

Xcovar <- NULL
if ("X" %in% names(cross2\$gmap)) {
    tryCatch({
        Xcovar <- get_x_covar(cross2)
    }, error = function(e) {
        cat("Warning: Could not get X covariates:", e\$message, "\\n")
    })
}

is_x_chr <- attr(alleleprob, "is_x_chr")
if (!is.null(is_x_chr)) {
    kinship_auto    <- kinship_loco[names(is_x_chr)[!is_x_chr]]
    alleleprob_auto <- alleleprob[, !is_x_chr]
} else {
    kinship_auto    <- kinship_loco
    alleleprob_auto <- alleleprob
}

pheno_data <- cross2\$pheno[, pheno, drop = FALSE]

cat("Running scan1perm(): 200 permutations, 48 cores (fork)...\\n")
perm_start <- Sys.time()

perm <- scan1perm(alleleprob_auto, pheno_data,
                  kinship  = kinship_auto,
                  addcovar = addcovar,
                  Xcovar   = Xcovar,
                  cores    = 48,
                  n_perm   = 200)

perm_time <- as.numeric(difftime(Sys.time(), perm_start, units = "mins"))
cat("scan1perm() complete in", round(perm_time, 2), "minutes\\n")

out <- data.frame(perm, check.names = FALSE)
saveRDS(out, out_file)

total_time <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))
cat("\\n===========================================\\n")
cat("BATCH COMPLETE\\n")
cat("Phenotype:", pheno, "  Batch:", batch_num, "/ 5\\n")
cat("Total time:", round(total_time, 2), "minutes\\n")
cat("End time:", format(Sys.time()), "\\n")
cat("===========================================\\n")
RSCRIPT

    log "Worker script written to: \${WORKER_SCRIPT}"

    # ── Read phenotype list ──────────────────────────────────────────────────
    PHENOTYPES=()
    while IFS= read -r line; do
        [[ -n "\${line}" ]] && PHENOTYPES+=("\${line}")
    done < "${phenotype_list}"

    NUM_PHENOS=\${#PHENOTYPES[@]}
    TOTAL_BATCHES=\$(( NUM_PHENOS * 5 ))
    log "Phenotypes: \${NUM_PHENOS}  |  Total batches: \${TOTAL_BATCHES}"

    # ── Constants ────────────────────────────────────────────────────────────
    MAX_CONCURRENT=200    # Max concurrent array tasks at once (%N in --array)
    MAX_ARRAY_SIZE=10000  # Ceres MaxArraySize limit (scontrol show config)
    MAX_RETRIES=2         # Max retry attempts per failed batch

    # ── Paths ────────────────────────────────────────────────────────────────
    STATUS_FILE="\${PERM_DIR}/COORDINATOR_STATUS.txt"
    ARRAY_WORKER="\${PERM_DIR}/perm_array_worker.sh"

    # ── Helper: check if a batch output file is complete ─────────────────────
    batch_done() {
        local F="\$1"
        [[ -f "\${F}" ]] && [[ \$(stat -c%s "\${F}" 2>/dev/null || echo 0) -gt 1000 ]]
    }

    # ── Helper: count completed .rds files in batch dir (fast, one ls call) ──
    count_done() {
        ls "\${BATCH_DIR}"/*.rds 2>/dev/null | wc -l
    }

    # ── Helper: update COORDINATOR_STATUS.txt ────────────────────────────────
    update_status_file() {
        local DONE="\$1" RUNNING="\$2" QUEUED="\$3" UNSUBMITTED="\$4" FAILED="\$5"
        local PCT=0
        [[ \${TOTAL_BATCHES} -gt 0 ]] && PCT=\$(( DONE * 100 / TOTAL_BATCHES ))
        cat > "\${STATUS_FILE}" <<STATUSEOF
=== PERMUTATION COORDINATOR STATUS ===
Updated:      \$(date)
Progress:     \${DONE} / \${TOTAL_BATCHES} batches complete (\${PCT}%)
------------------------------------------
Running:      \${RUNNING} (executing on compute nodes)
Queued:       \${QUEUED} (waiting in SLURM queue)
Unsubmitted:  \${UNSUBMITTED} (pending future array waves)
Perm failed:  \${FAILED} (exhausted retries)
STATUSEOF
    }

    # ── Helper: get running/queued task counts for an array job ──────────────
    get_array_counts() {
        local JID="\$1"
        local RUNNING=0 QUEUED=0
        while IFS= read -r ST; do
            case "\${ST}" in
                R|CG) RUNNING=\$((RUNNING+1)) ;;
                PD)   QUEUED=\$((QUEUED+1))   ;;
            esac
        done < <(squeue -j "\${JID}" -h --format="%t" 2>/dev/null)
        echo "\${RUNNING} \${QUEUED}"
    }

    # ── Write the array worker script ─────────────────────────────────────────
    # Each SLURM array task reads its PHENO/BATCH from the wave manifest using
    # \$SLURM_ARRAY_TASK_ID, then runs perm_batch_worker.R via apptainer.
    cat > "\${ARRAY_WORKER}" << 'ARRAYEOF'
#!/bin/bash
MANIFEST="\$1"
SIF="\$2"
WORKER="\$3"
STUDY="\$4"
PROJDIR="\$5"
OUTDIR="\$6"

PHENO=\$(awk -v id="\$SLURM_ARRAY_TASK_ID" '\$1==id{print \$2}' "\$MANIFEST")
BATCH=\$(awk -v id="\$SLURM_ARRAY_TASK_ID" '\$1==id{print \$3}' "\$MANIFEST")

if [[ -z "\$PHENO" || -z "\$BATCH" ]]; then
    echo "ERROR: task \$SLURM_ARRAY_TASK_ID not found in \$MANIFEST" >&2
    exit 1
fi

echo "Array task \$SLURM_ARRAY_TASK_ID -> \$PHENO batch \$BATCH"
apptainer exec --bind /project:/project "\$SIF" Rscript "\$WORKER" "\$PHENO" "\$BATCH" "\$STUDY" "\$PROJDIR" "\$OUTDIR"
ARRAYEOF
    chmod +x "\${ARRAY_WORKER}"
    log "Array worker written to: \${ARRAY_WORKER}"

    # ── Helper: submit one array wave ─────────────────────────────────────────
    # Args: wave_label manifest_file n_tasks
    # Prints the SLURM array job ID on success, empty string on failure
    submit_wave() {
        local LABEL="\$1" MANIFEST="\$2" N="\$3"
        local JOB_NAME="perm_\${STUDY_PREFIX}_\${LABEL}"
        JOB_NAME="\${JOB_NAME:0:63}"

        local JOB_ID
        JOB_ID=\$(sbatch \
            --account=do2_projects \
            --partition=ceres \
            --cpus-per-task=48 \
            --mem=200G \
            --time=2:00:00 \
            --job-name="\${JOB_NAME}" \
            --output="\${JOBS_LOG_DIR}/\${STUDY_PREFIX}_\${LABEL}_%a.log" \
            --array="1-\${N}%\${MAX_CONCURRENT}" \
            --wrap="bash '\${ARRAY_WORKER}' '\${MANIFEST}' '\${SIF}' '\${WORKER_SCRIPT}' '\${STUDY_PREFIX}' '\${PROJECT_DIR}' '\${OUTDIR}'" \
            2>/dev/null | awk '{print \$NF}')

        [[ "\${JOB_ID}" =~ ^[0-9]+\$ ]] && echo "\${JOB_ID}" || echo ""
    }

    # ── Build initial pending list ─────────────────────────────────────────────
    log "--- Scanning for completed/pending batches ---"
    PENDING=()
    declare -A RETRY_COUNT
    ALREADY_DONE=0

    for PHENO in "\${PHENOTYPES[@]}"; do
        for BATCH in 1 2 3 4 5; do
            if batch_done "\${BATCH_DIR}/\${STUDY_PREFIX}_\${PHENO}_batch_\${BATCH}.rds"; then
                ALREADY_DONE=\$((ALREADY_DONE+1))
            else
                PENDING+=("\${PHENO} \${BATCH}")
            fi
        done
    done

    log "Scan: \${ALREADY_DONE} already done, \${#PENDING[@]} need submitting"

    if [[ \${#PENDING[@]} -eq 0 ]]; then
        log "All batches complete — exiting"
        touch COORDINATOR_COMPLETE
        exit 0
    fi

    # ── Wave loop: submit PENDING in chunks of MAX_ARRAY_SIZE ─────────────────
    # Each wave is one sbatch --array submission. After it drains, failures are
    # retried in the next wave. Continues until PENDING is empty.
    WAVE=0
    PERM_FAILED=0

    while [[ \${#PENDING[@]} -gt 0 ]]; do
        WAVE=\$((WAVE+1))
        N_PENDING=\${#PENDING[@]}
        SLICE=\$(( N_PENDING < MAX_ARRAY_SIZE ? N_PENDING : MAX_ARRAY_SIZE ))

        # Write manifest for this wave: "task_id PHENO BATCH" per line
        MANIFEST="\${PERM_DIR}/wave_\${WAVE}_manifest.txt"
        for (( i=0; i<SLICE; i++ )); do
            echo "\$((i+1)) \${PENDING[\$i]}"
        done > "\${MANIFEST}"

        log "--- Wave \${WAVE}: \${SLICE} array tasks (of \${N_PENDING} remaining), \${MAX_CONCURRENT} concurrent ---"
        ARRAY_JID=\$(submit_wave "w\${WAVE}" "\${MANIFEST}" "\${SLICE}")

        if [[ -z "\${ARRAY_JID}" ]]; then
            log "ERROR: sbatch failed for wave \${WAVE} — aborting"
            exit 1
        fi
        log "Wave \${WAVE} submitted: array job \${ARRAY_JID}"

        # ── Monitor this wave until all tasks drain ────────────────────────────
        while true; do
            sleep 300
            read -r WV_RUNNING WV_QUEUED <<< "\$(get_array_counts "\${ARRAY_JID}")"
            DONE=\$(count_done)
            UNSUBMITTED=\$(( N_PENDING - SLICE ))

            PCT=0
            [[ \${TOTAL_BATCHES} -gt 0 ]] && PCT=\$(( DONE * 100 / TOTAL_BATCHES ))
            log "Wave \${WAVE} (job \${ARRAY_JID}) | \${DONE}/\${TOTAL_BATCHES} (\${PCT}%) | running: \${WV_RUNNING} | queued: \${WV_QUEUED} | unsubmitted: \${UNSUBMITTED} | failed: \${PERM_FAILED}"
            update_status_file "\${DONE}" "\${WV_RUNNING}" "\${WV_QUEUED}" "\${UNSUBMITTED}" "\${PERM_FAILED}"

            [[ \${WV_RUNNING} -eq 0 && \${WV_QUEUED} -eq 0 ]] && break
        done
        log "Wave \${WAVE} drained (job \${ARRAY_JID})"

        # ── Rebuild PENDING: failures from this wave + remainder for next wave ─
        NEXT_PENDING=()
        for (( i=0; i<SLICE; i++ )); do
            PHENO="\${PENDING[\$i]%% *}"
            BATCH="\${PENDING[\$i]##* }"
            if ! batch_done "\${BATCH_DIR}/\${STUDY_PREFIX}_\${PHENO}_batch_\${BATCH}.rds"; then
                KEY="\${PHENO}__b\${BATCH}"
                RETRIES=\${RETRY_COUNT["\${KEY}"]:-0}
                if [[ \${RETRIES} -lt \${MAX_RETRIES} ]]; then
                    RETRY_COUNT["\${KEY}"]=\$((RETRIES+1))
                    log "  Queue retry \$((RETRIES+1))/\${MAX_RETRIES}: \${PHENO} b\${BATCH}"
                    NEXT_PENDING+=("\${PHENO} \${BATCH}")
                else
                    log "  PERMANENTLY FAILED: \${PHENO} b\${BATCH}"
                    PERM_FAILED=\$((PERM_FAILED+1))
                fi
            fi
        done
        # Items beyond this wave's slice carry forward unchanged
        for (( i=SLICE; i<N_PENDING; i++ )); do
            NEXT_PENDING+=("\${PENDING[\$i]}")
        done

        if [[ \${#NEXT_PENDING[@]} -gt 0 ]]; then
            PENDING=("\${NEXT_PENDING[@]}")
        else
            PENDING=()
        fi

        log "After wave \${WAVE}: \${#PENDING[@]} batches remain (\${PERM_FAILED} permanent failures so far)"
    done

    if [[ \${PERM_FAILED} -gt 0 ]]; then
        log "ERROR: \${PERM_FAILED} batches permanently failed. Check logs in \${JOBS_LOG_DIR}"
        exit 1
    fi

    log "========================================================"
    log "ALL \${TOTAL_BATCHES} BATCHES COMPLETE"
    log "========================================================"
    echo "COORDINATOR_COMPLETE" > COORDINATOR_COMPLETE
    """
}


process PERM_AGGREGATE {
    tag "Aggregating permutation results"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '1h'

    input:
    path(completion_signal)   // COORDINATOR_COMPLETE trigger (or any existing file for perm_aggregate resume)
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

    log_message("=== AGGREGATING PERMUTATION RESULTS ===")
    log_message(paste("Trigger file:", "${completion_signal}"))

    # Read batch files directly from the published batch directory
    batch_dir <- file.path("${projectDir}", "${params.outdir}",
                           "07_permutation_testing", "batches")
    log_message(paste("Batch directory:", batch_dir))

    perm_files <- list.files(batch_dir,
                             pattern = paste0("^${study_prefix}_.*_batch_[0-9]+\\\\.rds\$"),
                             full.names = TRUE)
    log_message(paste("Found", length(perm_files), "batch result files"))

    if (length(perm_files) == 0) {
        log_message("ERROR: No batch result files found in batch directory")
        stop("No permutation results to aggregate")
    }

    log_message("Combining permutation matrices...")

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
// Session-independence strategy (no external submit scripts needed):
//
//   1. PERM_COORDINATOR is a 7-day SLURM job — completely independent of the
//      Nextflow driver session. It keeps running even if the driver dies.
//
//   2. On coordinator completion, it publishes COORDINATOR_COMPLETE to the
//      analysis output directory (persistent on shared filesystem).
//
//   3. If the Nextflow driver dies mid-run, simply re-run the same command
//      with -resume. The workflow checks for COORDINATOR_COMPLETE first:
//        - Found  → skip coordinator, run PERM_AGGREGATE immediately
//        - Missing → submit/resume PERM_COORDINATOR (its internal logic skips
//                    already-finished batches, so duplicate submission is safe)
//
//   Usage (any analysis type, no wrapper script needed):
//     nextflow run main_resume.nf --resume_from permutation \
//         --study_prefix <prefix> --analysis_type <type> -profile standard -resume
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
        chunks_per_batch         // Not used (always 1)

    main:
        def perm_dir           = "${params.outdir}/07_permutation_testing"
        def aggregated_perm    = "${perm_dir}/${params.study_prefix}_fullPerm.rds"
        def coordinator_done   = "${perm_dir}/COORDINATOR_COMPLETE"

        // Step 1: Setup — always runs to provide phenotype list and filtered cross2
        PERM_SETUP(
            cross2_file,
            filtered_phenotypes_file,
            study_prefix
        )

        if (file(aggregated_perm).exists()) {
            // ── Case A: fully done ─────────────────────────────────────────
            log.info "Aggregated permutation results already exist — skipping all permutation steps"
            ch_perm_matrix      = Channel.fromPath(aggregated_perm)
            ch_perm_thresholds  = Channel.fromPath("${perm_dir}/${params.study_prefix}_permThresh.rds")
            ch_filtered_phenos  = Channel.fromPath("${perm_dir}/${params.study_prefix}_filtered_phenotypes_lod${params.lod_threshold}.txt")
            ch_aggregation_log  = Channel.fromPath("${perm_dir}/${params.study_prefix}_aggregation_log.txt")

        } else if (file(coordinator_done).exists()) {
            // ── Case B: coordinator finished but driver died before aggregation ──
            // COORDINATOR_COMPLETE is on disk — all batches are done.
            // Skip re-submitting the coordinator; go straight to aggregation.
            log.info "COORDINATOR_COMPLETE found — all batches done. Running PERM_AGGREGATE."
            ch_agg_trigger = Channel.fromPath(coordinator_done)

            PERM_AGGREGATE(
                ch_agg_trigger,
                study_prefix,
                lod_threshold
            )

            ch_perm_matrix      = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds  = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenos  = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log  = PERM_AGGREGATE.out.aggregation_log

        } else {
            // ── Case C: normal run or driver-died-mid-coordinator ─────────
            // Submit (or re-submit) the coordinator. Its internal resume logic
            // skips batches whose output files already exist, so re-submission
            // after a driver death is safe and efficient.
            PERM_COORDINATOR(
                PERM_SETUP.out.phenotype_list,
                study_prefix
            )

            PERM_AGGREGATE(
                PERM_COORDINATOR.out.completion_signal,
                study_prefix,
                lod_threshold
            )

            ch_perm_matrix      = PERM_AGGREGATE.out.full_perm_matrix
            ch_perm_thresholds  = PERM_AGGREGATE.out.perm_thresholds
            ch_filtered_phenos  = PERM_AGGREGATE.out.significant_phenotypes
            ch_aggregation_log  = PERM_AGGREGATE.out.aggregation_log
        }

    emit:
        permutation_matrix     = ch_perm_matrix
        permutation_thresholds = ch_perm_thresholds
        filtered_phenotypes    = ch_filtered_phenos
        filtered_cross2        = PERM_SETUP.out.filtered_cross2
        setup_log              = PERM_SETUP.out.setup_log
        aggregation_log        = ch_aggregation_log
}
