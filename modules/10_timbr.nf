// Module 10: TIMBR Allelic Series Analysis
// Runs TIMBR (Crouse et al. 2020, Genetics 216:957-983) on significant QTLs
// to infer the number of functional alleles among the 8 DO founder strains.

process TIMBR_SETUP {
    tag "TIMBR setup for ${prefix} (sig_level: ${sig_level})"
    publishDir "${params.outdir}/10_timbr", mode: 'copy'

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
    path("batch_ids.txt"),                 emit: batch_ids

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
        # Strip coat_color (genetic phenotype, not a covariate for the model)
        if ("coat_color" %in% colnames(covar_data)) {
            covar_data <- covar_data[, !colnames(covar_data) %in% "coat_color", drop = FALSE]
        }
        # Remove constant covariates — model.matrix crashes with "contrasts can be
        # applied only to factors with 2 or more levels" for single-level factors
        # (e.g., Sex in an all-female or all-male study).
        constant_cols <- sapply(covar_data, function(x) length(unique(na.omit(x))) < 2)
        if (any(constant_cols)) {
            setup_log <- c(setup_log, paste("  Removing constant covariate(s):", paste(names(which(constant_cols)), collapse = ", ")))
            covar_data <- covar_data[, !constant_cols, drop = FALSE]
        }
        if (ncol(covar_data) > 0) {
            covar_formula <- paste("~", paste(colnames(covar_data), collapse = " + "))
            addcovar <- model.matrix(as.formula(covar_formula), data = covar_data)[, -1, drop = FALSE]
            if (ncol(addcovar) == 0) addcovar <- NULL
        } else {
            addcovar <- NULL
        }
        if (!is.null(addcovar)) {
            setup_log <- c(setup_log, paste("✓ Addcovar matrix built:", nrow(addcovar), "samples x", ncol(addcovar), "columns"))
            setup_log <- c(setup_log, paste("  Covariates:", paste(colnames(covar_data), collapse = ", ")))
        } else {
            setup_log <- c(setup_log, "✓ No variable covariates after filtering — addcovar = NULL")
        }
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

    # Write batch IDs to file (one per line) for channel splitting
    writeLines(as.character(1:n_batches), "batch_ids.txt")
    """
}


process TIMBR_BATCH {
    tag "TIMBR batch ${batch_id} for ${prefix}"
    publishDir "${params.outdir}/10_timbr", mode: 'copy'

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
    path("chr*/*/plots/*.png"),                           emit: timbr_plots, optional: true
    path("chr*/*/*.rds"),                                 emit: timbr_rds,   optional: true

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(qtl2)
        library(TIMBR)
        library(ape)
        library(parallel)
    })

    n_cores <- ${task.cpus}
    cat("Batch ${batch_id}: using", n_cores, "cores via mclapply\\n")

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
    pmap        <- cross2\$pmap                  # physical map (Mb), for cM→Mb lookup

    # Load batch-specific QTLs
    all_qtls  <- read.csv("${filtered_qtls_file}", stringsAsFactors = FALSE)
    batch_map <- read.table("${batch_map_file}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    this_batch_ids <- batch_map\$lodcolumn[batch_map\$batch_id == ${batch_id}]
    this_batch     <- all_qtls[all_qtls\$lodcolumn %in% this_batch_ids & all_qtls\$batch_id == ${batch_id}, ]

    batch_log <- c(batch_log, paste("✓ QTLs in this batch:", nrow(this_batch)))

    # Output dirs are created per-QTL inside the loop (chr-specific structure)

    # ── Tree-informed prior (DO/CC founder phylogeny, tip labels = qtl2 order) ─
    # A=A/J  B=C57BL/6J  C=129S1  D=NOD  E=NZO  F=CAST  G=PWK  H=WSB
    # Tree from Crouse et al. 2020 / TIMBR mcv.data example
    do_tree <- ape::read.tree(text = "((B:0.05537190085,(E:0.01310357241,D:0.01310357241):0.04226832844):1.227652821,((F:0.4288873004,((C:0.03940317437,A:0.03940317437):0.03135662222,H:0.07075979659):0.3581275038):0.200533186,G:0.6294204864):0.6536042357);")
    hp      <- TIMBR::calc.concentration.prior(8, 0.05, 0.01)
    prior.M <- list(model.type        = "tree",
                    tree              = do_tree,
                    prior.alpha.type  = "gamma",
                    prior.alpha.shape = hp[1],
                    prior.alpha.rate  = hp[2],
                    hash.names        = TRUE)

    # ── Identify common samples across genoprob, phenotype, and addcovar ──────
    geno_samples  <- rownames(genoprob[[1]])
    pheno_samples <- rownames(cross2\$pheno)
    if (!is.null(addcovar)) {
        common_samples <- Reduce(intersect, list(geno_samples, pheno_samples, rownames(addcovar)))
    } else {
        common_samples <- intersect(geno_samples, pheno_samples)
    }
    batch_log <- c(batch_log, paste("✓ Common samples:", length(common_samples)))

    # ── Per-QTL function (called via mclapply) ────────────────────────────────
    process_qtl <- function(i) {
        row       <- this_batch[i, ]
        pheno_id  <- row\$lodcolumn
        chr       <- as.character(row\$chr)
        pos_cM    <- row\$pos
        safe_id   <- gsub("[^A-Za-z0-9_-]", "_", pheno_id)

        cat(sprintf("  [%d/%d] %s  chr%s  %.2f cM\\n", i, nrow(this_batch), pheno_id, chr, pos_cM))

        result_row <- data.frame(
            prefix               = "${prefix}",
            lodcolumn            = pheno_id,
            chr                  = chr,
            pos                  = pos_cM,
            peak_Mb              = NA_real_,
            lod                  = row\$lod,
            significance_level   = row\$significance_level,
            nearest_marker       = NA_character_,
            ln_BF                = NA_real_,
            n_functional_alleles = NA_integer_,
            top_series_pattern   = NA_character_,
            top_series_prob      = NA_real_,
            n_samples_used       = NA_integer_,
            status               = "pending",
            stringsAsFactors     = FALSE
        )

        tryCatch({
            # ── Find nearest pseudomarker ──────────────────────────────────
            if (!chr %in% names(genetic_map)) stop(paste("Chr", chr, "not in genetic_map"))
            gmap_chr    <- genetic_map[[chr]]
            nearest_mk  <- names(which.min(abs(gmap_chr - pos_cM)))

            # ── Convert genetic position (cM) to physical position (Mb) ───
            if (!chr %in% names(pmap)) stop(paste("Chr", chr, "not in pmap"))
            peak_Mb   <- round(pmap[[chr]][nearest_mk], 2)
            qtl_label <- paste0(safe_id, "_chr", chr, "_", peak_Mb, "Mb")
            chr_dir   <- paste0("chr", chr)
            qtl_dir   <- file.path(chr_dir, qtl_label)
            plots_dir <- file.path(qtl_dir, "plots")
            dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

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

            # ── Residualize phenotype against covariates ───────────────────
            # Passing Z directly to TIMBR() can cause positive-definite errors;
            # residualizing first is the standard workaround.
            if (!is.null(Z)) {
                fit <- lm(y ~ Z)
                y   <- residuals(fit)
            }

            # ── Run TIMBR ─────────────────────────────────────────────────
            timbr_result <- TIMBR::TIMBR(y, prior.D, prior.M,
                                          Z         = NULL,
                                          samples   = ${timbr_samples},
                                          calc.lnBF = TRUE)

            # ── Extract summary statistics ─────────────────────────────────
            # MAP allelic series: the pattern with highest posterior probability
            map_pattern    <- names(which.max(timbr_result\$p.M.given.y))
            map_prob       <- max(timbr_result\$p.M.given.y)
            n_func_alleles <- length(unique(as.integer(strsplit(map_pattern, ",")[[1]])))

            # ── Save RDS ───────────────────────────────────────────────────
            rds_fname <- file.path(qtl_dir, paste0(qtl_label, ".rds"))
            saveRDS(timbr_result, rds_fname)

            # ── Plot title base ────────────────────────────────────────────
            plot_title <- paste0(pheno_id, "  chr", chr, "  ", peak_Mb, " Mb")

            # ── Plot 1: Posterior haplotype effects ────────────────────────
            tryCatch({
                png(file.path(plots_dir, "01_haplotype_effects.png"), width = 960, height = 480)
                TIMBR::TIMBR.plot.haplotypes(timbr_result)
                title(main = plot_title)
                dev.off()
            }, error = function(e) { try(dev.off(), silent=TRUE); cat("  Plot 1 failed:", conditionMessage(e), "\\n") })

            # ── Plot 2: Posterior distribution of number of functional alleles
            tryCatch({
                n_alleles    <- sapply(names(timbr_result\$p.M.given.y),
                                       function(x) max(as.numeric(unlist(strsplit(x, ",")))) + 1)
                allele_probs <- tapply(timbr_result\$p.M.given.y, n_alleles, sum)
                png(file.path(plots_dir, "02_n_alleles_distribution.png"), width = 600, height = 500)
                barplot(allele_probs,
                        xlab = "Number of functional alleles",
                        ylab = "Posterior probability",
                        main = plot_title,
                        col  = "steelblue")
                dev.off()
            }, error = function(e) { try(dev.off(), silent=TRUE); cat("  Plot 2 failed:", conditionMessage(e), "\\n") })

            # ── Plot 3: Top 10 posterior allelic series ────────────────────
            tryCatch({
                top_n   <- min(10, length(timbr_result\$p.M.given.y))
                top_ser <- head(timbr_result\$p.M.given.y, top_n)
                png(file.path(plots_dir, "03_top_allelic_series.png"), width = 800, height = 500)
                barplot(rev(top_ser), horiz = TRUE, las = 1, cex.names = 0.8,
                        xlab = "Posterior probability",
                        ylab = "Allelic series",
                        main = plot_title,
                        col  = "steelblue")
                dev.off()
            }, error = function(e) { try(dev.off(), silent=TRUE); cat("  Plot 3 failed:", conditionMessage(e), "\\n") })

            # ── Plot 4: TIMBR allelic series vs qtl2 founder effects ───────
            # Compares TIMBR posterior haplotype effects (colored by MAP allelic
            # series) against qtl2-style regression-on-probabilities (ROP) effects.
            # Both use the same residualized phenotype and haplotype dosage matrix.
            tryCatch({
                founder_names <- c("A/J", "C57BL/6J", "129S1", "NOD", "NZO", "CAST", "PWK", "WSB")

                # ── TIMBR posterior mean and 95% CI ──────────────────────────
                timbr_means <- colMeans(timbr_result\$post.hap.effects)
                timbr_lo    <- apply(timbr_result\$post.hap.effects, 2, quantile, 0.025)
                timbr_hi    <- apply(timbr_result\$post.hap.effects, 2, quantile, 0.975)
                # Centre around mean effect
                timbr_means <- timbr_means - mean(timbr_means)
                timbr_lo    <- timbr_lo    - mean(colMeans(timbr_result\$post.hap.effects))
                timbr_hi    <- timbr_hi    - mean(colMeans(timbr_result\$post.hap.effects))

                # ── qtl2 ROP founder effects via regression on dosages ────────
                # H = P %*% A: [n_valid x 8] expected haplotype dosages.
                # Rows of H sum to 1, so we omit the intercept (y is already
                # residualized / roughly zero-mean) to avoid perfect collinearity.
                H <- prior.D\$P %*% prior.D\$A
                colnames(H) <- founder_names
                fit_rop  <- lm(y ~ H - 1)
                rop_eff  <- coef(fit_rop)
                rop_se   <- summary(fit_rop)\$coefficients[, "Std. Error"]
                rop_eff  <- rop_eff - mean(rop_eff, na.rm = TRUE)   # centre

                # ── MAP allelic series colours (one colour per allele group) ─
                series_idx  <- as.integer(strsplit(map_pattern, ",")[[1]]) + 1L
                n_groups    <- max(series_idx)
                group_cols  <- rainbow(n_groups, s = 0.8, v = 0.85)
                founder_cols <- group_cols[series_idx]

                # ── Plot ──────────────────────────────────────────────────────
                png(file.path(plots_dir, "04_allelic_series_vs_founder_effects.png"),
                    width = 1000, height = 520)
                op <- par(mfrow = c(1, 2), mar = c(7, 4, 4, 1), oma = c(0, 0, 3, 0))

                ylim_all <- range(c(timbr_lo, timbr_hi,
                                    rop_eff - 2*rop_se, rop_eff + 2*rop_se), na.rm = TRUE)

                # Panel A: TIMBR posterior effects
                bp <- barplot(timbr_means, names.arg = founder_names,
                              col = founder_cols, las = 2,
                              ylab = "Effect (centred)", ylim = ylim_all,
                              main = "TIMBR posterior effects")
                arrows(bp, timbr_lo, bp, timbr_hi,
                       angle = 90, code = 3, length = 0.05, lwd = 1.5)
                abline(h = 0, lty = 2, col = "grey50")

                # Panel B: qtl2 ROP founder effects
                bp2 <- barplot(rop_eff, names.arg = founder_names,
                               col = founder_cols, las = 2,
                               ylab = "Effect (centred)", ylim = ylim_all,
                               main = "qtl2 founder effects (ROP)")
                arrows(bp2, rop_eff - 2*rop_se, bp2, rop_eff + 2*rop_se,
                       angle = 90, code = 3, length = 0.05, lwd = 1.5)
                abline(h = 0, lty = 2, col = "grey50")

                # Shared title and allele-group legend
                mtext(plot_title, outer = TRUE, cex = 1.1, font = 2)
                par(op)
                dev.off()
            }, error = function(e) { try(dev.off(), silent=TRUE); cat("  Plot 4 failed:", conditionMessage(e), "\\n") })

            # ── Populate result row ────────────────────────────────────────
            result_row\$peak_Mb              <- peak_Mb
            result_row\$nearest_marker       <- nearest_mk
            result_row\$ln_BF               <- timbr_result\$ln.BF
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

        result_row
    }

    # ── Run in parallel ───────────────────────────────────────────────────────
    results_list <- parallel::mclapply(seq_len(nrow(this_batch)), process_qtl, mc.cores = n_cores)

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
    csv_files <- sort(list.files(".", pattern = "_timbr_batch[0-9]+_summary[.]csv\$", full.names = TRUE))
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
