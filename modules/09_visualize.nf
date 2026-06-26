// Module 9: QTL Visualization (batched + parallelized)
// Splits significant QTLs into SLURM batches; each batch uses mclapply internally.

process VISUALIZE_SETUP {
    tag "Visualize setup for ${prefix}"
    publishDir "${params.outdir}/09_visualize", mode: 'copy'

    cpus   2
    memory '8 GB'
    time   '30m'

    input:
    path(significant_qtls_file)
    val(prefix)
    val(qtls_per_batch)

    output:
    path("viz_batch_qtls.csv"),  emit: batch_qtls
    path("viz_batch_ids.txt"),   emit: batch_ids
    path("viz_setup_log.txt"),   emit: setup_log

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dplyr))

    sig_qtls <- read.csv("${significant_qtls_file}", stringsAsFactors = FALSE)
    qtls_99  <- sig_qtls %>% filter(significance_level == "99%")

    if (nrow(qtls_99) == 0) stop("No QTLs at 99% significance found.")

    qtls_99 <- qtls_99 %>%
        mutate(viz_batch_id = ceiling(row_number() / ${qtls_per_batch}))

    n_batches <- max(qtls_99\$viz_batch_id)

    write.csv(qtls_99, "viz_batch_qtls.csv", row.names = FALSE)
    writeLines(as.character(1:n_batches), "viz_batch_ids.txt")

    log_lines <- c(
        "=== Visualize Setup Report ===",
        paste("Study prefix:", "${prefix}"),
        paste("Timestamp:", Sys.time()),
        paste("Total QTLs at 99%:", nrow(qtls_99)),
        paste("QTLs per batch:", ${qtls_per_batch}),
        paste("Total batches:", n_batches)
    )
    writeLines(log_lines, "viz_setup_log.txt")
    cat(paste(log_lines, collapse = "\\n"), "\\n")
    """
}


process VISUALIZE_BATCH {
    tag "Visualize batch ${batch_id} for ${prefix}"
    publishDir "${params.outdir}/09_visualize", mode: 'copy'

    input:
    val(batch_id)
    path(alleleprob_file)
    path(scan_results_file)
    path(cross2_file)
    path(kinship_file)
    path(batch_qtls_file)
    val(prefix)
    val(interactive_covar)

    output:
    path("chr*/*/*/*.png"), emit: qtl_plots, optional: true
    path("batch_${batch_id}_report.txt"), emit: batch_report

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
        library(parallel)
    })

    n_cores <- ${task.cpus}
    cat("Batch ${batch_id}: using", n_cores, "cores\\n")

    # ── Load shared data once ─────────────────────────────────────────────────
    cat("Loading shared data...\\n")
    alleleprobs  <- readRDS("${alleleprob_file}")
    scan_results <- readRDS("${scan_results_file}")
    cross2       <- readRDS("${cross2_file}")
    kinship      <- readRDS("${kinship_file}")

    pmap        <- cross2\$pmap
    gmap        <- cross2\$gmap
    HALF_WIN_MB <- 4

    # Order chromosomes naturally (1..19, X) for genome-wide plots. The stored
    # map is in lexical order (10,11,...,19,1,2,...) which otherwise places chr10
    # first on the genome-wide LOD plot. Per-chr/window lookups are name-based and
    # unaffected; only plot_scan1's genome-wide x-axis depends on map order.
    chr_order <- c(as.character(1:99), "X", "Y", "M")
    chr_order <- c(chr_order[chr_order %in% names(pmap)], setdiff(names(pmap), chr_order))
    pmap      <- pmap[chr_order]

    # ── Build covariate matrices ──────────────────────────────────────────────
    interactive_covar_name <- "${interactive_covar}"
    covar_data <- cross2\$covar

    addcovar <- NULL
    if (!is.null(covar_data) && ncol(covar_data) > 0) {
        variable_cols <- sapply(covar_data, function(x) length(unique(na.omit(x))) >= 2)
        covar_for_model <- covar_data[, variable_cols, drop = FALSE]
        if (ncol(covar_for_model) > 0) {
            covar_formula <- paste("~", paste(colnames(covar_for_model), collapse = " + "))
            addcovar <- model.matrix(as.formula(covar_formula), data = covar_for_model)[, -1, drop = FALSE]
        }
    }

    do_stratify <- interactive_covar_name != "null" && interactive_covar_name != "" &&
                   !is.null(covar_data) && interactive_covar_name %in% colnames(covar_data)

    addcovar_sex <- NULL
    intcovar_mat <- NULL
    diet_levels  <- NULL

    if (do_stratify) {
        diet_levels  <- sort(unique(covar_data[[interactive_covar_name]]))
        non_int_cols <- setdiff(colnames(covar_data), interactive_covar_name)
        non_int_cols <- non_int_cols[sapply(covar_data[, non_int_cols, drop = FALSE],
                                            function(x) length(unique(na.omit(x))) >= 2)]
        if (length(non_int_cols) > 0) {
            sex_formula  <- paste("~", paste(non_int_cols, collapse = " + "))
            addcovar_sex <- model.matrix(as.formula(sex_formula), data = covar_data)[, -1, drop = FALSE]
        }
        int_formula  <- paste("~", interactive_covar_name)
        intcovar_mat <- model.matrix(as.formula(int_formula), data = covar_data)[, -1, drop = FALSE]
    }

    # ── Load batch-specific QTLs ──────────────────────────────────────────────
    all_batch_qtls <- read.csv("${batch_qtls_file}", stringsAsFactors = FALSE)
    qtls_99 <- all_batch_qtls %>% filter(viz_batch_id == ${batch_id})
    cat("Batch ${batch_id}:", nrow(qtls_99), "QTLs to process\\n")

    # Pre-subset scan_results to only batch-relevant columns before forking.
    # scan_results for eQTL runs can be tens of GB; without this every mclapply
    # worker receives a full copy, causing OOM kills.
    needed_cols <- intersect(qtls_99\$lodcolumn, colnames(scan_results))
    scan_results <- scan_results[, needed_cols, drop = FALSE]
    cat("scan_results subsetted to", ncol(scan_results), "columns for this batch\\n")

    # Create chromosome output directories up front
    for (chr in unique(qtls_99\$chr)) {
        dir.create(paste0("chr", chr), showWarnings = FALSE, recursive = TRUE)
    }

    # ── Per-QTL plot function (called via mclapply) ───────────────────────────
    process_qtl <- function(i) {
        qtl     <- qtls_99[i, ]
        chr     <- as.character(qtl\$chr)
        gene_id <- qtl\$lodcolumn
        pos_cM  <- qtl\$pos
        lod     <- qtl\$lod
        safe_id <- gsub("[^A-Za-z0-9_-]", "_", gene_id)

        plots_ok   <- 0L
        plots_fail <- 0L

        tryCatch({
            if (!gene_id %in% colnames(cross2\$pheno)) stop(paste("Phenotype", gene_id, "not found"))
            if (!gene_id %in% colnames(scan_results))  stop(paste("Phenotype", gene_id, "not in scan_results"))

            # Peak marker and ±HALF_WIN_MB window
            gmap_chr    <- gmap[[chr]]
            peak_marker <- names(gmap_chr)[which.min(abs(gmap_chr - pos_cM))]
            peak_Mb     <- pmap[[chr]][peak_marker]
            win_lo      <- peak_Mb - HALF_WIN_MB
            win_hi      <- peak_Mb + HALF_WIN_MB

            win_markers <- names(pmap[[chr]])[pmap[[chr]] >= win_lo & pmap[[chr]] <= win_hi]
            if (length(win_markers) < 2) stop("Too few markers in window")

            ap_chr        <- alleleprobs[, chr]
            ap_win        <- ap_chr
            ap_win[[chr]] <- ap_chr[[chr]][, , win_markers, drop = FALSE]
            pmap_win      <- list(); pmap_win[[chr]] <- pmap[[chr]][win_markers]

            # On-disk layout: chr{N}/{phenotype}/{lod,coef,blup}/{cis,chr,genome}.png
            chr_dir   <- paste0("chr", chr)
            qtl_dir   <- file.path(chr_dir, safe_id)
            lod_dir   <- file.path(qtl_dir, "lod")
            coef_dir  <- file.path(qtl_dir, "coef")
            dir.create(lod_dir,  showWarnings = FALSE, recursive = TRUE)
            dir.create(coef_dir, showWarnings = FALSE, recursive = TRUE)

            pmap_full_chr <- pmap[chr]   # named list holding the single full chromosome

            # ── LOD plots: cis (±window) / chr / genome ───────────────────────
            png(file.path(lod_dir, "cis.png"), width = 900, height = 400)
            par(mar = c(4.5, 4.5, 2.5, 0.5))
            plot_scan1(scan_results, pmap_win, lodcolumn = gene_id, bgcolor = "gray95",
                       xlab = paste("Chr", chr, "position (Mb)"), ylab = "LOD score")
            abline(v = peak_Mb, col = "red", lty = 3, lwd = 1.5)
            title(main = sprintf("LOD [cis +/-%g Mb] — %s  Chr %s  (peak %.2f Mb, LOD %.2f)",
                                 HALF_WIN_MB, gene_id, chr, peak_Mb, lod))
            dev.off(); plots_ok <- plots_ok + 1L

            png(file.path(lod_dir, "chr.png"), width = 900, height = 400)
            par(mar = c(4.5, 4.5, 2.5, 0.5))
            plot_scan1(scan_results, pmap, lodcolumn = gene_id, chr = chr, bgcolor = "gray95",
                       xlab = paste("Chr", chr, "position (Mb)"), ylab = "LOD score")
            abline(v = peak_Mb, col = "red", lty = 3, lwd = 1.5)
            title(main = sprintf("LOD [chr] — %s  Chr %s  (peak %.2f Mb, LOD %.2f)",
                                 gene_id, chr, peak_Mb, lod))
            dev.off(); plots_ok <- plots_ok + 1L

            # genome: no chr= → all chromosomes; no peak line (Mb does not map
            # cleanly across the concatenated multi-chromosome x-axis).
            png(file.path(lod_dir, "genome.png"), width = 1400, height = 400)
            par(mar = c(4.5, 4.5, 2.5, 0.5))
            plot_scan1(scan_results, pmap, lodcolumn = gene_id, bgcolor = "gray95",
                       ylab = "LOD score")
            title(main = sprintf("LOD [genome] — %s  (peak Chr %s %.2f Mb, LOD %.2f)",
                                 gene_id, chr, peak_Mb, lod))
            dev.off(); plots_ok <- plots_ok + 1L

            # ── Effects: restrict to genotyped samples ────────────────────────
            # cross2\$covar may contain samples that were phenotyped but lack
            # genotype data (absent from the FinalReport chip files). Kinship and
            # alleleprob only cover genotyped samples, so subsetting with the full
            # covar rownames causes "subscript out of bounds" on kinship[[chr]].
            genotyped_ids <- rownames(alleleprobs[[names(alleleprobs)[1]]])
            ids   <- intersect(rownames(covar_data), genotyped_ids)
            ac_g  <- if (!is.null(addcovar) && nrow(addcovar) > 0)
                         addcovar[intersect(ids, rownames(addcovar)), , drop = FALSE] else NULL
            ic_g  <- if (do_stratify && !is.null(intcovar_mat) && nrow(intcovar_mat) > 0)
                         intcovar_mat[intersect(ids, rownames(intcovar_mat)), , drop = FALSE] else NULL
            kin_g       <- kinship[[chr]][ids, ids]
            pheno_mat   <- cross2\$pheno[ids, gene_id, drop = FALSE]
            model_label <- if (!is.null(ic_g)) "GxE intcovar" else "addcovar only"

            # Helper: render one founder-effect (coef/BLUP) plot
            effect_plot <- function(eff_obj, map_obj, fname, kind, level_lab) {
                png(fname, width = 900, height = 500)
                par(mar = c(4.5, 4.5, 3.5, 0.5))
                plot_coefCC(eff_obj, map_obj, bgcolor = "gray95", legend = "bottomleft",
                            xlab = paste("Chr", chr, "position (Mb)"))
                abline(v = peak_Mb, col = "red", lty = 3, lwd = 1.5)
                title(main = sprintf("%s [%s] — %s  (N=%d)", level_lab, kind, gene_id, length(ids)))
                mtext(sprintf("peak %s (%.2f Mb) | %s", peak_marker, peak_Mb, model_label),
                      side = 3, line = 0.1, cex = 0.65, col = "gray40", font = 3)
                dev.off()
            }

            # ── coef: cis (±window) / chr ─────────────────────────────────────
            # Founder effects are estimated per-chromosome; no genome-wide level.
            coef_cis <- scan1coef(ap_win[ids, chr], pheno_mat, kinship = kin_g, addcovar = ac_g, intcovar = ic_g)
            effect_plot(coef_cis, pmap_win, file.path(coef_dir, "cis.png"), "coef",
                        sprintf("cis +/-%g Mb", HALF_WIN_MB)); plots_ok <- plots_ok + 1L
            coef_chr <- scan1coef(ap_chr[ids, chr], pheno_mat, kinship = kin_g, addcovar = ac_g, intcovar = ic_g)
            effect_plot(coef_chr, pmap_full_chr, file.path(coef_dir, "chr.png"), "coef", "chr")
            plots_ok <- plots_ok + 1L

            # ── blup: cis (±window) / chr — additive runs only ────────────────
            # BLUP cannot use the interaction term; for interactive runs these are
            # produced by the diet-specific analyses (AIN76_only, HC_only, Diet_additive).
            if (!do_stratify) {
                blup_dir <- file.path(qtl_dir, "blup")
                dir.create(blup_dir, showWarnings = FALSE, recursive = TRUE)
                blup_cis <- scan1blup(ap_win[ids, chr], pheno_mat, kinship = kin_g, addcovar = ac_g)
                effect_plot(blup_cis, pmap_win, file.path(blup_dir, "cis.png"), "BLUP",
                            sprintf("cis +/-%g Mb", HALF_WIN_MB)); plots_ok <- plots_ok + 1L
                blup_chr <- scan1blup(ap_chr[ids, chr], pheno_mat, kinship = kin_g, addcovar = ac_g)
                effect_plot(blup_chr, pmap_full_chr, file.path(blup_dir, "chr.png"), "BLUP", "chr")
                plots_ok <- plots_ok + 1L
            }

        }, error = function(e) {
            plots_fail <<- plots_fail + 1L
            cat(sprintf("WARNING [batch ${batch_id}]: failed for %s chr%s: %s\\n", gene_id, chr, e\$message))
        })

        list(ok = plots_ok, fail = plots_fail)
    }

    # ── Run in parallel ───────────────────────────────────────────────────────
    cat("Running mclapply across", nrow(qtls_99), "QTLs on", n_cores, "cores...\\n")
    results <- parallel::mclapply(seq_len(nrow(qtls_99)), process_qtl, mc.cores = n_cores)

    total_ok   <- sum(sapply(results, `[[`, "ok"))
    total_fail <- sum(sapply(results, `[[`, "fail"))

    report <- c(
        paste("=== Visualize Batch ${batch_id} Report ==="),
        paste("Timestamp:", Sys.time()),
        paste("QTLs processed:", nrow(qtls_99)),
        paste("Plots generated:", total_ok),
        paste("Plots failed:", total_fail)
    )
    writeLines(report, "batch_${batch_id}_report.txt")
    cat(sprintf("Batch ${batch_id} complete: %d plots generated, %d failed\\n", total_ok, total_fail))
    """
}


process VISUALIZE_AGGREGATE {
    tag "Visualize aggregate for ${prefix}"
    publishDir "${params.outdir}/09_visualize", mode: 'copy'

    cpus   2
    memory '8 GB'
    time   '30m'

    input:
    path(batch_reports)
    val(prefix)

    output:
    path("visualization_report.txt"), emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript

    report_files <- sort(list.files(".", pattern = "^batch_[0-9]+_report[.]txt\$", full.names = TRUE))
    cat("Aggregating", length(report_files), "batch reports\\n")

    all_lines  <- unlist(lapply(report_files, readLines))
    total_ok   <- sum(as.integer(gsub(".*Plots generated: ([0-9]+).*", "\\\\1",
                                      grep("Plots generated", all_lines, value = TRUE))))
    total_fail <- sum(as.integer(gsub(".*Plots failed: ([0-9]+).*", "\\\\1",
                                      grep("Plots failed",    all_lines, value = TRUE))))

    summary <- c(
        "=== QTL Visualization Report ===",
        paste("Study Prefix:", "${prefix}"),
        paste("Timestamp:", Sys.time()),
        paste("Total batches:", length(report_files)),
        paste("Total plots generated:", total_ok),
        if (total_fail > 0) paste("Total plots failed:", total_fail) else NULL,
        "",
        "=== QTL Visualization Complete ==="
    )
    writeLines(summary, "visualization_report.txt")
    cat(paste(summary, collapse = "\\n"), "\\n")
    """
}


workflow VISUALIZE_QTLS {
    take:
    alleleprob_ch
    scan_results_ch
    cross2_ch
    significant_qtls_ch
    prefix_ch
    kinship_ch
    interactive_covar_ch

    main:
    // Clear previous visualization results so old and new plots do not mix.
    // Runs on the head node before any publishDir copy, so freshly published
    // outputs from this run are unaffected.
    def vizDir = file("${params.outdir}/09_visualize")
    if (vizDir.exists()) {
        log.info "Clearing previous visualization results at ${vizDir}"
        vizDir.deleteDir()
    }

    VISUALIZE_SETUP(
        significant_qtls_ch,
        prefix_ch,
        Channel.value(params.viz_qtls_per_batch)
    )

    VISUALIZE_SETUP.out.batch_ids
        .splitText()
        .map  { it.trim() }
        .filter { it != "" }
        .set  { ch_batch_ids }

    VISUALIZE_BATCH(
        ch_batch_ids,
        alleleprob_ch.first(),
        scan_results_ch.first(),
        cross2_ch.first(),
        kinship_ch.first(),
        VISUALIZE_SETUP.out.batch_qtls.first(),
        prefix_ch,
        interactive_covar_ch
    )

    VISUALIZE_AGGREGATE(
        VISUALIZE_BATCH.out.batch_report.collect(),
        prefix_ch
    )

    emit:
    qtl_plots         = VISUALIZE_BATCH.out.qtl_plots
    validation_report = VISUALIZE_AGGREGATE.out.validation_report
}
