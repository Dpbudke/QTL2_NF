#!/usr/bin/env nextflow

// Parameterized permutation coordinator
// Runs as a single long-lived SLURM job (30-day wall time) that submits
// and monitors batch worker jobs for one stage of the two-stage permutation flow.
//
// Two callers (in modules/07_permutation_testing.nf):
//   PERM_COORDINATOR_SCREEN: stage_label='screen', n_perms=50,  n_batches=1
//   PERM_COORDINATOR_FULL:   stage_label='full',   n_perms=200, n_batches=5
//
// Each call gets its own batch dir, worker scripts, status/log files,
// and completion signal (COORDINATOR_${stage_label}_COMPLETE).

process PERM_COORDINATOR {
    tag "Coordinating ${stage_label} permutations (30-day SLURM job)"
    publishDir "${params.outdir}/07_permutation_testing", mode: 'copy'

    // No container — this process needs host SLURM tools (sbatch, sacct).
    container null

    cpus   1
    memory '4 GB'
    time   '30d'
    errorStrategy 'retry'
    maxRetries 0

    input:
    path(phenotype_list)
    val(study_prefix)
    val(interactive_covar)
    val(stage_label)         // 'screen' or 'full'
    val(n_perms_per_batch)   // 50 (screen) or 200 (full)
    val(n_batches_per_pheno) // 1  (screen) or 5   (full)

    output:
    path("${study_prefix}_${stage_label}_coordinator_log.txt"), emit: coordinator_log
    path("COORDINATOR_${stage_label}_COMPLETE"),                emit: completion_signal

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # ── Paths ────────────────────────────────────────────────────────────────
    STUDY_PREFIX="${study_prefix}"
    PROJECT_DIR="${projectDir}"
    OUTDIR="${params.outdir}"
    INTERACTIVE_COVAR="${interactive_covar}"
    STAGE_LABEL="${stage_label}"
    N_PERMS=${n_perms_per_batch}
    N_BATCHES=${n_batches_per_pheno}
    SIF="\${PROJECT_DIR}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"
    PERM_DIR="\${PROJECT_DIR}/\${OUTDIR}/07_permutation_testing"
    BATCH_DIR="\${PERM_DIR}/batches_\${STAGE_LABEL}"
    JOBS_LOG_DIR="\${PERM_DIR}/coordinator_job_logs_\${STAGE_LABEL}"
    WORKER_SCRIPT="\${PERM_DIR}/perm_batch_worker_\${STAGE_LABEL}.R"
    LOG="${study_prefix}_${stage_label}_coordinator_log.txt"

    mkdir -p "\${BATCH_DIR}" "\${JOBS_LOG_DIR}"

    log() { echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$*" | tee -a "\${LOG}"; }

    log "========================================================"
    log "PERMUTATION COORDINATOR STARTED — stage: \${STAGE_LABEL}"
    log "Study:       \${STUDY_PREFIX}"
    log "Project:     \${PROJECT_DIR}"
    log "Output dir:  \${OUTDIR}"
    log "Batch dir:   \${BATCH_DIR}"
    log "Container:   \${SIF}"
    log "Int covar:   \${INTERACTIVE_COVAR}"
    log "Perms/batch: \${N_PERMS}"
    log "Batches/ph:  \${N_BATCHES}"
    log "========================================================"

    # ── Sanity checks ────────────────────────────────────────────────────────
    if ! command -v sbatch &>/dev/null; then
        log "ERROR: sbatch not found"
        exit 1
    fi
    if [[ ! -f "\${SIF}" ]]; then
        log "ERROR: Singularity image not found: \${SIF}"
        exit 1
    fi

    # ── Write R worker script ────────────────────────────────────────────────
    # Stage-specific filename prevents the screen worker overwriting the full
    # worker (or vice versa) when both coordinators publish to PERM_DIR.
    cat > "\${WORKER_SCRIPT}" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
    stop("Usage: perm_batch_worker.R <pheno> <batch_num> <study_prefix> <project_dir> <outdir> <interactive_covar> <stage_label> <n_perms>")
}
pheno                  <- args[1]
batch_num              <- as.integer(args[2])
study_pfx              <- args[3]
proj_dir               <- args[4]
outdir                 <- args[5]
interactive_covar_name <- args[6]
stage_label            <- args[7]
n_perms                <- as.integer(args[8])

suppressPackageStartupMessages(library(qtl2))

batch_start <- Sys.time()
cat("===========================================\\n")
cat("STAGE:", stage_label, "  BATCH:", batch_num, "\\n")
cat("PHENOTYPE:", pheno, "\\n")
cat("PERMUTATIONS:", n_perms, "\\n")
cat("START TIME:", format(batch_start), "\\n")
cat("===========================================\\n\\n")

perm_dir  <- file.path(proj_dir, outdir, "07_permutation_testing")
batch_dir <- file.path(perm_dir,  paste0("batches_", stage_label))
scan_dir  <- file.path(proj_dir, outdir, "05_genome_scan_preparation")
out_file  <- file.path(batch_dir,
                        paste0(study_pfx, "_", pheno, "_batch_", batch_num, ".rds"))

# Skip if already done (handles duplicate submissions on restart)
if (file.exists(out_file) && file.size(out_file) > 100) {
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

# Covariates — split into additive and interactive matrices to match scan1 in module 6
addcovar <- NULL
intcovar <- NULL

if (!is.null(cross2\$covar)) {
    covar_data <- cross2\$covar
    if ("coat_color" %in% colnames(covar_data)) {
        covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
    }
    constant_cols <- sapply(covar_data, function(x) length(unique(na.omit(x))) < 2)
    if (any(constant_cols)) {
        cat("NOTE: Removing constant covariate(s) with no variation:", paste(names(which(constant_cols)), collapse=", "), "\\n")
        covar_data <- covar_data[, !constant_cols, drop = FALSE]
    }

    use_intcovar <- interactive_covar_name != "null" && interactive_covar_name != "" &&
                    interactive_covar_name %in% colnames(covar_data)

    if (use_intcovar) {
        intcovar_data <- covar_data[, interactive_covar_name, drop = FALSE]
        addcovar_data <- covar_data[, !colnames(covar_data) %in% interactive_covar_name, drop = FALSE]

        intcovar <- model.matrix(as.formula(paste("~", interactive_covar_name)),
                                 data = intcovar_data)[, -1, drop = FALSE]
        if (ncol(addcovar_data) > 0) {
            addcovar <- model.matrix(as.formula(paste("~", paste(colnames(addcovar_data), collapse = " + "))),
                                     data = addcovar_data)[, -1, drop = FALSE]
            if (ncol(addcovar) == 0) addcovar <- NULL
        }
        cat("Interactive covariate:", interactive_covar_name, "(", ncol(intcovar), "cols);",
            "additive covariates:", if (is.null(addcovar)) 0 else ncol(addcovar), "cols\\n")
    } else {
        if (ncol(covar_data) > 0) {
            addcovar <- model.matrix(as.formula(paste("~", paste(colnames(covar_data), collapse = " + "))),
                                     data = covar_data)[, -1, drop = FALSE]
            if (ncol(addcovar) == 0) addcovar <- NULL
        }
    }
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

cat("Running scan1perm():", n_perms, "permutations, 48 cores (fork)",
    if (!is.null(intcovar)) "(with intcovar)" else "(additive only)", "...\\n")
perm_start <- Sys.time()

perm <- scan1perm(alleleprob_auto, pheno_data,
                  kinship  = kinship_auto,
                  addcovar = addcovar,
                  intcovar = intcovar,
                  Xcovar   = Xcovar,
                  cores    = 48,
                  n_perm   = n_perms)

perm_time <- as.numeric(difftime(Sys.time(), perm_start, units = "mins"))
cat("scan1perm() complete in", round(perm_time, 2), "minutes\\n")

out <- data.frame(perm, check.names = FALSE)
saveRDS(out, out_file)

total_time <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))
cat("\\n===========================================\\n")
cat("BATCH COMPLETE — stage:", stage_label, "\\n")
cat("Phenotype:", pheno, "  Batch:", batch_num, "\\n")
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
    TOTAL_BATCHES=\$(( NUM_PHENOS * N_BATCHES ))
    log "Phenotypes: \${NUM_PHENOS}  |  Batches/pheno: \${N_BATCHES}  |  Total: \${TOTAL_BATCHES}"

    # ── Constants ────────────────────────────────────────────────────────────
    MAX_CONCURRENT=200
    MAX_ARRAY_SIZE=10000
    MAX_RETRIES=2

    # ── Stage-specific paths ─────────────────────────────────────────────────
    STATUS_FILE="\${PERM_DIR}/COORDINATOR_STATUS_\${STAGE_LABEL}.txt"
    ARRAY_WORKER="\${PERM_DIR}/perm_array_worker_\${STAGE_LABEL}.sh"

    # ── Helpers ──────────────────────────────────────────────────────────────
    batch_done() {
        local F="\$1"
        [[ -f "\${F}" ]] && [[ \$(stat -c%s "\${F}" 2>/dev/null || echo 0) -gt 100 ]]
    }

    count_done() {
        find "\${BATCH_DIR}" -maxdepth 1 -name "*.rds" 2>/dev/null | wc -l
    }

    update_status_file() {
        local DONE="\$1" RUNNING="\$2" QUEUED="\$3" UNSUBMITTED="\$4" FAILED="\$5"
        local PCT=0
        [[ \${TOTAL_BATCHES} -gt 0 ]] && PCT=\$(( DONE * 100 / TOTAL_BATCHES ))
        cat > "\${STATUS_FILE}" <<STATUSEOF
=== PERMUTATION COORDINATOR STATUS (\${STAGE_LABEL}) ===
Updated:      \$(date)
Progress:     \${DONE} / \${TOTAL_BATCHES} batches complete (\${PCT}%)
------------------------------------------
Running:      \${RUNNING} (executing on compute nodes)
Queued:       \${QUEUED} (waiting in SLURM queue)
Unsubmitted:  \${UNSUBMITTED} (pending future array waves)
Perm failed:  \${FAILED} (exhausted retries)
STATUSEOF
    }

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

    release_held_jobs() {
        local JID="\$1"
        local HELD
        HELD=\$(squeue -j "\${JID}" -h --format="%r" 2>/dev/null | grep -c "held" || true)
        if [[ \${HELD} -gt 0 ]]; then
            log "WARNING: \${HELD} task(s) in held state for job \${JID} — running scontrol release"
            scontrol release "\${JID}" 2>/dev/null || true
        fi
    }

    # ── Write the array worker script ─────────────────────────────────────────
    cat > "\${ARRAY_WORKER}" << 'ARRAYEOF'
#!/bin/bash
MANIFEST="\$1"
SIF="\$2"
WORKER="\$3"
STUDY="\$4"
PROJDIR="\$5"
OUTDIR="\$6"
INTCOVAR="\$7"
STAGE="\$8"
NPERMS="\$9"

PHENO=\$(awk -v id="\$SLURM_ARRAY_TASK_ID" '\$1==id{print \$2}' "\$MANIFEST")
BATCH=\$(awk -v id="\$SLURM_ARRAY_TASK_ID" '\$1==id{print \$3}' "\$MANIFEST")

if [[ -z "\$PHENO" || -z "\$BATCH" ]]; then
    echo "ERROR: task \$SLURM_ARRAY_TASK_ID not found in \$MANIFEST" >&2
    exit 1
fi

echo "Array task \$SLURM_ARRAY_TASK_ID -> \$PHENO batch \$BATCH (stage=\$STAGE intcovar=\$INTCOVAR n_perms=\$NPERMS)"
apptainer exec --bind /project:/project "\$SIF" Rscript "\$WORKER" "\$PHENO" "\$BATCH" "\$STUDY" "\$PROJDIR" "\$OUTDIR" "\$INTCOVAR" "\$STAGE" "\$NPERMS"
ARRAYEOF
    chmod +x "\${ARRAY_WORKER}"
    log "Array worker written to: \${ARRAY_WORKER}"

    # ── Helper: submit one array wave ─────────────────────────────────────────
    submit_wave() {
        local LABEL="\$1" MANIFEST="\$2" N="\$3"
        local JOB_NAME="perm_\${STUDY_PREFIX}_\${STAGE_LABEL}_\${LABEL}"
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
            --wrap="bash '\${ARRAY_WORKER}' '\${MANIFEST}' '\${SIF}' '\${WORKER_SCRIPT}' '\${STUDY_PREFIX}' '\${PROJECT_DIR}' '\${OUTDIR}' '\${INTERACTIVE_COVAR}' '\${STAGE_LABEL}' '\${N_PERMS}'" \
            2>/dev/null | awk '{print \$NF}')

        [[ "\${JOB_ID}" =~ ^[0-9]+\$ ]] && echo "\${JOB_ID}" || echo ""
    }

    # ── Build initial pending list ─────────────────────────────────────────────
    log "--- Scanning for completed/pending batches ---"
    PENDING=()
    declare -A RETRY_COUNT
    ALREADY_DONE=0

    for PHENO in "\${PHENOTYPES[@]}"; do
        for (( BATCH=1; BATCH<=N_BATCHES; BATCH++ )); do
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
        touch COORDINATOR_\${STAGE_LABEL}_COMPLETE
        exit 0
    fi

    # ── Wave loop ─────────────────────────────────────────────────────────────
    WAVE=0
    PERM_FAILED=0

    while [[ \${#PENDING[@]} -gt 0 ]]; do
        WAVE=\$((WAVE+1))
        N_PENDING=\${#PENDING[@]}
        SLICE=\$(( N_PENDING < MAX_ARRAY_SIZE ? N_PENDING : MAX_ARRAY_SIZE ))

        MANIFEST="\${PERM_DIR}/wave_\${STAGE_LABEL}_\${WAVE}_manifest.txt"
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

        while true; do
            sleep 300
            release_held_jobs "\${ARRAY_JID}"
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
    log "ALL \${TOTAL_BATCHES} \${STAGE_LABEL} BATCHES COMPLETE"
    log "========================================================"
    echo "COORDINATOR_\${STAGE_LABEL}_COMPLETE" > COORDINATOR_\${STAGE_LABEL}_COMPLETE

    # Eager-publish: copy completion signal + coordinator log directly to the
    # published results dir. Nextflow's publishDir requires the driver to be
    # alive at process exit; if the driver died during this 30-day run, these
    # files would otherwise be stranded in the work dir and the next resume
    # would not detect the stage as complete. Workers/manifests/batches/status
    # are already written directly into PERM_DIR, so this closes the only gap.
    cp -f "COORDINATOR_\${STAGE_LABEL}_COMPLETE" "\${PERM_DIR}/" 2>/dev/null && \
        log "Eager-published COORDINATOR_\${STAGE_LABEL}_COMPLETE to \${PERM_DIR}/" || \
        log "WARNING: failed to eager-publish COORDINATOR_\${STAGE_LABEL}_COMPLETE"
    cp -f "\${LOG}" "\${PERM_DIR}/" 2>/dev/null && \
        log "Eager-published \${LOG} to \${PERM_DIR}/" || \
        log "WARNING: failed to eager-publish \${LOG}"
    """
}
