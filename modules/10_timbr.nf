// Module 10: TIMBR Allelic Series Analysis
// Runs TIMBR (Crouse et al. 2020, Genetics 216:957-983) on significant QTLs
// to infer the number of functional alleles among the 8 DO founder strains.

process TIMBR_SETUP {
    tag "TIMBR setup for ${prefix} (sig_level: ${sig_level})"
    publishDir "${params.outdir}/10_timbr", mode: 'copy',
        saveAs: { fn -> fn.endsWith("_timbr_qtls.csv") || fn.endsWith("_timbr_setup_log.txt") ? fn : null }

    input:
    path(cross2_file)
    path(significant_qtls_file)
    val(prefix)
    val(sig_level)
    val(qtls_per_batch)

    output:
    path("${prefix}_timbr_qtls.csv"),      emit: filtered_qtls
    path("${prefix}_timbr_batch_map.tsv"), emit: batch_file
    path("${prefix}_addcovar.rds"),        emit: addcovar_rds
    path("${prefix}_timbr_setup_log.txt"), emit: setup_log
    stdout,                                emit: batch_ids

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })

    setup_log <- c()
    setup_log <- c(setup_log, "=== TIMBR Setup Report ===")
    setup_log <- c(setup_log, paste("Timestamp:", Sys.time()))
    setup_log <- c(setup_log, paste("Study Prefix:", "${prefix}"))
    setup_log <- c(setup_log, paste("Significance Filter:", "${sig_level}"))
    setup_log <- c(setup_log, paste("QTLs per Batch:", ${qtls_per_batch}))
    setup_log <- c(setup_log, "")

    cat("Loading data...\\n")
    cross2 <- readRDS("${cross2_file}")
    sig_qtls <- read.csv("${significant_qtls_file}", stringsAsFactors = FALSE)

    setup_log <- c(setup_log, paste("✓ Total QTL rows loaded:", nrow(sig_qtls)))

    # ── Filter to requested significance level ────────────────────────────────
    # significance_level order: 63% < 90% < 95% < 99%
    sig_order <- c("63%", "90%", "95%", "99%")
    requested_idx <- which(sig_order == "${sig_level}")
    if (length(requested_idx) == 0) stop("Invalid timbr_sig_level: ${sig_level}. Valid: 63%, 90%, 95%, 99%")
    keep_levels <- sig_order[requested_idx:length(sig_order)]

    filtered <- sig_qtls %>%
        filter(significance_level %in% keep_levels)

    setup_log <- c(setup_log, paste("✓ QTLs at or above ${sig_level}:", nrow(filtered)))

    if (nrow(filtered) == 0) {
        stop("No QTLs pass the significance filter '${sig_level}'. Check timbr_sig_level parameter.")
    }

    # ── Deduplicate: keep highest significance row per (lodcolumn, chr, pos) ──
    filtered <- filtered %>%
        mutate(.sig_rank = match(significance_level, rev(sig_order))) %>%
        group_by(lodcolumn, chr, pos) %>%
        slice_max(.sig_rank, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(-.sig_rank) %>%
        arrange(factor(chr, levels = c(as.character(1:19), "X")), pos)

    setup_log <- c(setup_log, paste("✓ Unique QTLs after deduplication:", nrow(filtered)))

    # ── Assign batch IDs ──────────────────────────────────────────────────────
    filtered <- filtered %>%
        mutate(batch_id = ceiling(row_number() / ${qtls_per_batch}))

    n_batches <- max(filtered\$batch_id)
    setup_log <- c(setup_log, paste("✓ Batches created:", n_batches, "(", ${qtls_per_batch}, "QTLs/batch)"))
    setup_log <- c(setup_log, "")

    # ── Pre-compute addcovar matrix ───────────────────────────────────────────
    if (!is.null(cross2\$covar) && ncol(cross2\$covar) > 0) {
        covar_data <- cross2\$covar
        # Strip coat_color (categorical, not a proper covariate for model)
        if ("coat_color" %in% colnames(covar_data)) {
            covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
        }
        covar_formula <- paste("~", paste(colnames(covar_data), collapse = " + "))
        addcovar <- model.matrix(as.formula(covar_formula), data = covar_data)[, -1, drop = FALSE]
        setup_log <- c(setup_log, paste("✓ Addcovar matrix built:", nrow(addcovar), "samples x", ncol(addcovar), "columns"))
        setup_log <- c(setup_log, paste("  Covariates:", paste(colnames(covar_data), collapse = ", ")))
    } else {
        addcovar <- NULL
        setup_log <- c(setup_log, "✓ No covariates found — addcovar = NULL")
    }
    saveRDS(addcovar, "${prefix}_addcovar.rds")

    # ── Write outputs ─────────────────────────────────────────────────────────
    write.csv(filtered, "${prefix}_timbr_qtls.csv", row.names = FALSE)
    write.table(filtered[, c("batch_id", "lodcolumn", "chr", "pos", "lod", "significance_level")],
                "${prefix}_timbr_batch_map.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    setup_log <- c(setup_log, "=== Setup Complete ===")
    setup_log <- c(setup_log, paste("✓ Filtered QTL list:", "${prefix}_timbr_qtls.csv"))
    setup_log <- c(setup_log, paste("✓ Batch map:", "${prefix}_timbr_batch_map.tsv"))
    setup_log <- c(setup_log, paste("✓ Addcovar RDS:", "${prefix}_addcovar.rds"))
    writeLines(setup_log, "${prefix}_timbr_setup_log.txt")

    # Emit batch IDs to stdout (one per line) for channel splitting
    cat(paste(1:n_batches, collapse = "\\n"), "\\n")
    """
}


process TIMBR_BATCH {
    tag "TIMBR batch ${batch_id} for ${prefix}"
    publishDir "${params.outdir}/10_timbr/batches", mode: 'copy'

    input:
    val(batch_id)
    path(genoprob_file)
    path(genetic_map_file)
    path(filtered_qtls_file)
    path(batch_map_file)
    path(cross2_file)
    path(addcovar_file)
    val(prefix)
    val(timbr_samples)

    output:
    path("${prefix}_timbr_batch${batch_id}_summary.csv"), emit: batch_summary
    path("${prefix}_timbr_batch${batch_id}_log.txt"),     emit: batch_log
    path("timbr_plots_batch${batch_id}/*.png"),           emit: timbr_plots, optional: true
    path("timbr_rds_batch${batch_id}/*.rds"),             emit: timbr_rds,   optional: true

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(qtl2)
        library(TIMBR)
    })

    set.seed(${batch_id} * 12345)

    batch_log <- c()
    batch_log <- c(batch_log, paste("=== TIMBR Batch ${batch_id} Log ==="))
    batch_log <- c(batch_log, paste("Timestamp:", Sys.time()))
    batch_log <- c(batch_log, paste("MCMC samples:", ${timbr_samples}))
    batch_log <- c(batch_log, "")

    cat("Loading data for batch ${batch_id}...\\n")
    genoprob    <- readRDS("${genoprob_file}")
    genetic_map <- readRDS("${genetic_map_file}")
    cross2      <- readRDS("${cross2_file}")
    addcovar    <- readRDS("${addcovar_file}")   # NULL if no covariates

    # Load batch-specific QTLs
    all_qtls  <- read.csv("${filtered_qtls_file}", stringsAsFactors = FALSE)
    batch_map <- read.table("${batch_map_file}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    this_batch_ids <- batch_map\$lodcolumn[batch_map\$batch_id == ${batch_id}]
    this_batch     <- all_qtls[all_qtls\$lodcolumn %in% this_batch_ids & all_qtls\$batch_id == ${batch_id}, ]

    batch_log <- c(batch_log, paste("✓ QTLs in this batch:", nrow(this_batch)))

    # Create output directories
    dir.create("timbr_plots_batch${batch_id}", showWarnings = FALSE)
    dir.create("timbr_rds_batch${batch_id}",   showWarnings = FALSE)

    # ── CRP prior (computed once, reused for all QTLs in batch) ───────────────
    hp      <- TIMBR::calc.concentration.prior(8, 0.05, 0.01)
    prior.M <- list(model.type        = "crp",
                    prior.alpha.type  = "gamma",
                    prior.alpha.shape = hp[1],
                    prior.alpha.rate  = hp[2])

    # ── Identify common samples across genoprob, phenotype, and addcovar ──────
    geno_samples  <- rownames(genoprob[[1]])
    pheno_samples <- rownames(cross2\$pheno)
    if (!is.null(addcovar)) {
        common_samples <- Reduce(intersect, list(geno_samples, pheno_samples, rownames(addcovar)))
    } else {
        common_samples <- intersect(geno_samples, pheno_samples)
    }
    batch_log <- c(batch_log, paste("✓ Common samples:", length(common_samples)))

    # ── Results container ─────────────────────────────────────────────────────
    results_list <- vector("list", nrow(this_batch))

    for (i in seq_len(nrow(this_batch))) {
        row       <- this_batch[i, ]
        pheno_id  <- row\$lodcolumn
        chr       <- as.character(row\$chr)
        pos_cM    <- row\$pos
        safe_id   <- gsub("[^A-Za-z0-9_-]", "_", pheno_id)

        cat(sprintf("  [%d/%d] %s  chr%s  %.2f cM\\n", i, nrow(this_batch), pheno_id, chr, pos_cM))

        result_row <- data.frame(
            prefix              = "${prefix}",
            lodcolumn           = pheno_id,
            chr                 = chr,
            pos                 = pos_cM,
            lod                 = row\$lod,
            significance_level  = row\$significance_level,
            nearest_marker      = NA_character_,
            ln_BF               = NA_real_,
            n_functional_alleles = NA_integer_,
            top_series_pattern  = NA_character_,
            top_series_prob     = NA_real_,
            n_samples_used      = NA_integer_,
            status              = "pending",
            stringsAsFactors    = FALSE
        )

        tryCatch({
            # ── Find nearest pseudomarker ──────────────────────────────────
            if (!chr %in% names(genetic_map)) stop(paste("Chr", chr, "not in genetic_map"))
            gmap_chr    <- genetic_map[[chr]]
            nearest_mk  <- names(which.min(abs(gmap_chr - pos_cM)))

            # ── Extract P matrix [individuals x 36 diplotype states] ───────
            if (!chr %in% names(genoprob)) stop(paste("Chr", chr, "not in genoprob"))
            P_full <- genoprob[[chr]][, , nearest_mk]

            # Align to common samples
            samples_here <- intersect(common_samples, rownames(P_full))
            if (length(samples_here) < 10) stop(paste("Too few common samples:", length(samples_here)))

            P_aligned <- P_full[samples_here, ]

            # ── Build prior.D ──────────────────────────────────────────────
            prior.D <- list(
                P           = P_aligned,
                A           = TIMBR::additive.design(8, "qtl2"),
                fixed.diplo = FALSE
            )

            # ── Align phenotype and addcovar ───────────────────────────────
            if (!pheno_id %in% colnames(cross2\$pheno)) stop(paste("Phenotype", pheno_id, "not in cross2"))
            y <- cross2\$pheno[samples_here, pheno_id]
            Z <- if (!is.null(addcovar)) addcovar[samples_here, , drop = FALSE] else NULL

            # Remove samples with NA phenotype
            valid <- !is.na(y)
            if (sum(valid) < 10) stop(paste("Too few non-NA samples:", sum(valid)))
            y <- y[valid]
            prior.D\$P <- prior.D\$P[valid, ]
            if (!is.null(Z)) Z <- Z[valid, , drop = FALSE]

            # ── Run TIMBR ─────────────────────────────────────────────────
            timbr_result <- TIMBR::TIMBR(y, prior.D, prior.M,
                                          Z       = Z,
                                          samples = ${timbr_samples},
                                          calc.lnBF = TRUE)

            # ── Extract summary statistics ─────────────────────────────────
            # MAP allelic series: the pattern with highest posterior probability
            map_pattern      <- names(which.max(timbr_result\$p.M.given.y))
            map_prob         <- max(timbr_result\$p.M.given.y)
            n_func_alleles   <- length(unique(as.integer(strsplit(map_pattern, ",")[[1]])))

            # ── Save RDS and plot ──────────────────────────────────────────
            rds_fname  <- file.path("timbr_rds_batch${batch_id}",
                                     paste0(safe_id, "_chr", chr, "_", round(pos_cM, 1), "cM.rds"))
            plot_fname <- file.path("timbr_plots_batch${batch_id}",
                                     paste0(safe_id, "_chr", chr, "_", round(pos_cM, 1), "cM.png"))

            saveRDS(timbr_result, rds_fname)

            png(plot_fname, width = 800, height = 600)
            TIMBR::TIMBR.plot.haplotypes(timbr_result,
                                          main = paste0(pheno_id, "  chr", chr, "  ", round(pos_cM, 1), " cM"))
            dev.off()

            # ── Populate result row ────────────────────────────────────────
            result_row\$nearest_marker       <- nearest_mk
            result_row\$ln_BF                <- timbr_result\$ln.BF
            result_row\$n_functional_alleles <- n_func_alleles
            result_row\$top_series_pattern   <- map_pattern
            result_row\$top_series_prob      <- map_prob
            result_row\$n_samples_used       <- sum(valid)
            result_row\$status               <- "success"

        }, error = function(e) {
            msg <- paste("ERROR:", conditionMessage(e))
            cat("  WARNING:", msg, "\\n")
            result_row\$status <<- msg
        })

        results_list[[i]] <- result_row
    }

    # ── Write batch summary ───────────────────────────────────────────────────
    batch_summary <- do.call(rbind, results_list)
    write.csv(batch_summary, "${prefix}_timbr_batch${batch_id}_summary.csv", row.names = FALSE)

    n_success <- sum(batch_summary\$status == "success")
    n_fail    <- sum(batch_summary\$status != "success")
    batch_log <- c(batch_log, paste("✓ Completed:", n_success, "successes,", n_fail, "failures"))
    batch_log <- c(batch_log, paste("✓ End time:", Sys.time()))
    writeLines(batch_log, "${prefix}_timbr_batch${batch_id}_log.txt")

    cat(sprintf("Batch ${batch_id} complete: %d/%d succeeded\\n", n_success, nrow(this_batch)))
    """
}


process TIMBR_AGGREGATE {
    tag "TIMBR aggregate for ${prefix}"
    publishDir "${params.outdir}/10_timbr", mode: 'copy'

    input:
    path(batch_summaries)
    path(batch_logs)
    val(prefix)

    output:
    path("${prefix}_timbr_master_summary.csv"), emit: master_summary
    path("${prefix}_timbr_run_report.txt"),     emit: run_report

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(dplyr))

    cat("Aggregating TIMBR batch results...\\n")

    # ── Read and combine all batch summary CSVs ───────────────────────────────
    csv_files <- sort(list.files(".", pattern = "_timbr_batch[0-9]+_summary\\.csv\$", full.names = TRUE))
    cat("Found", length(csv_files), "batch summary files\\n")

    all_results <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
    master <- do.call(rbind, all_results)

    # Sort by chromosome then position
    chr_order <- c(as.character(1:19), "X")
    master <- master %>%
        mutate(chr = factor(chr, levels = chr_order)) %>%
        arrange(chr, pos) %>%
        mutate(chr = as.character(chr))

    write.csv(master, "${prefix}_timbr_master_summary.csv", row.names = FALSE)
    cat("Master summary written:", nrow(master), "QTLs\\n")

    # ── Run report ────────────────────────────────────────────────────────────
    n_total   <- nrow(master)
    n_success <- sum(master\$status == "success", na.rm = TRUE)
    n_fail    <- n_total - n_success

    report <- c()
    report <- c(report, "=== TIMBR Run Report ===")
    report <- c(report, paste("Timestamp:", Sys.time()))
    report <- c(report, paste("Study prefix:", "${prefix}"))
    report <- c(report, "")
    report <- c(report, "=== Summary ===")
    report <- c(report, paste("Total QTLs attempted:", n_total))
    report <- c(report, paste("Successful:", n_success))
    report <- c(report, paste("Failed:", n_fail))
    report <- c(report, "")

    if (n_success > 0) {
        success <- master %>% filter(status == "success")
        report <- c(report, "=== Functional Allele Distribution ===")
        tbl <- table(success\$n_functional_alleles)
        for (k in names(tbl)) {
            report <- c(report, paste(" ", k, "allele(s):", tbl[[k]], "QTLs"))
        }
        report <- c(report, "")
        report <- c(report, "=== Results by Significance Level ===")
        sig_tbl <- success %>% count(significance_level) %>% arrange(significance_level)
        for (i in seq_len(nrow(sig_tbl))) {
            report <- c(report, paste(" ", sig_tbl\$significance_level[i], ":", sig_tbl\$n[i], "QTLs"))
        }
        report <- c(report, "")
        report <- c(report, paste("Median ln(BF):", round(median(success\$ln_BF, na.rm = TRUE), 3)))
        report <- c(report, paste("Median top series prob:", round(median(success\$top_series_prob, na.rm = TRUE), 3)))
    }

    writeLines(report, "${prefix}_timbr_run_report.txt")
    cat("Run report written.\\n")
    """
}


workflow TIMBR_ANALYSIS {
    take:
    cross2_ch
    genoprob_ch
    genetic_map_ch
    significant_qtls_ch
    prefix_ch
    sig_level_ch
    qtls_per_batch_ch
    timbr_samples_ch

    main:
    TIMBR_SETUP(
        cross2_ch,
        significant_qtls_ch,
        prefix_ch,
        sig_level_ch,
        qtls_per_batch_ch
    )

    TIMBR_SETUP.out.batch_ids
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
        .set { ch_batch_ids }

    TIMBR_BATCH(
        ch_batch_ids,
        genoprob_ch.first(),
        genetic_map_ch.first(),
        TIMBR_SETUP.out.filtered_qtls.first(),
        TIMBR_SETUP.out.batch_file.first(),
        cross2_ch.first(),
        TIMBR_SETUP.out.addcovar_rds.first(),
        prefix_ch,
        timbr_samples_ch
    )

    TIMBR_AGGREGATE(
        TIMBR_BATCH.out.batch_summary.collect(),
        TIMBR_BATCH.out.batch_log.collect(),
        prefix_ch
    )

    emit:
    master_summary = TIMBR_AGGREGATE.out.master_summary
    run_report     = TIMBR_AGGREGATE.out.run_report
    timbr_plots    = TIMBR_BATCH.out.timbr_plots
    setup_log      = TIMBR_SETUP.out.setup_log
}
