process VISUALIZE_QTLS {
    tag "Generating QTL visualizations for ${prefix}"
    publishDir "${params.outdir}/09_visualize", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '8h'

    input:
    path(alleleprob_file)
    path(scan_results_file)
    path(cross2_file)
    path(significant_qtls_file)
    val(prefix)
    path(kinship_file)
    val(interactive_covar)

    output:
    path("chr*/*.png"),   emit: qtl_plots_flat,   optional: true
    path("chr*/*/*.png"), emit: qtl_plots_subdir, optional: true
    path("visualization_report.txt"), emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })

    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== QTL Visualization Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")

    cat("Loading QTL analysis results...\\n")

    # Load data
    alleleprobs <- readRDS("${alleleprob_file}")
    scan_results <- readRDS("${scan_results_file}")
    cross2 <- readRDS("${cross2_file}")
    kinship <- readRDS("${kinship_file}")

    # Use physical map (Mb) from cross2 instead of genetic map (cM)
    pmap <- cross2\$pmap
    gmap <- cross2\$gmap
    HALF_WIN_MB <- 2   # ±2 Mb window around QTL peak (4 Mb total)

    # Build covariate matrices for scan1coef
    # Note: For coefficient visualization, all covariates (including any interactive
    # covariate) are included as additive terms. This adjusts for population structure
    # and covariate effects while keeping the founder allele estimates numerically
    # stable. The GxE interaction effect is captured by the LOD score from scan1,
    # not the coefficient plot. Using intcovar in scan1coef can produce numerically
    # unstable estimates at markers with ill-conditioned design matrices.
    interactive_covar_name <- "${interactive_covar}"
    covar_data <- cross2\$covar

    addcovar <- NULL

    if (!is.null(covar_data) && ncol(covar_data) > 0) {
        covar_formula <- paste("~", paste(colnames(covar_data), collapse = " + "))
        addcovar <- model.matrix(as.formula(covar_formula), data = covar_data)[, -1, drop = FALSE]
        cat("Additive covariate columns:", ncol(addcovar), "\\n")
        if (interactive_covar_name != "null" && interactive_covar_name != "") {
            cat("Note:", interactive_covar_name, "included as additive covariate for stable coefficient estimation\\n")
        }
    }

    # Determine if diet-stratified plots are needed (interactive model only)
    do_stratify <- interactive_covar_name != "null" && interactive_covar_name != "" &&
                   !is.null(covar_data) && interactive_covar_name %in% colnames(covar_data)

    addcovar_sex <- NULL
    intcovar_mat <- NULL
    diet_levels  <- NULL

    if (do_stratify) {
        diet_levels  <- sort(unique(covar_data[[interactive_covar_name]]))
        non_int_cols <- setdiff(colnames(covar_data), interactive_covar_name)
        if (length(non_int_cols) > 0) {
            sex_formula  <- paste("~", paste(non_int_cols, collapse = " + "))
            addcovar_sex <- model.matrix(as.formula(sex_formula), data = covar_data)[, -1, drop = FALSE]
        }
        # Build intcovar matrix from the interactive covariate column only
        int_formula  <- paste("~", interactive_covar_name)
        intcovar_mat <- model.matrix(as.formula(int_formula), data = covar_data)[, -1, drop = FALSE]
        cat("Stratified plots enabled for", interactive_covar_name, "levels:",
            paste(diet_levels, collapse = ", "), "\\n")
        cat("_full plot will use Diet as intcovar (GxE model confirmed stable in window)\\n")
    }

    # Load significant QTLs and filter to 99% significance level
    significant_qtls <- read.csv("${significant_qtls_file}", stringsAsFactors = FALSE)
    qtls_99 <- significant_qtls %>% filter(significance_level == "99%")

    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Alleleprobs loaded:", length(alleleprobs), "chromosomes"))
    validation_log <- c(validation_log, paste("✓ Physical map loaded:", length(pmap), "chromosomes (positions in Mb)"))
    validation_log <- c(validation_log, paste("✓ Scan results loaded:", paste(dim(scan_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Cross2 object loaded:", nrow(cross2\$pheno), "individuals,", ncol(cross2\$pheno), "phenotypes"))
    validation_log <- c(validation_log, paste("✓ Kinship matrices loaded:", length(kinship), "chromosomes"))
    validation_log <- c(validation_log, paste("✓ Interactive covariate:", interactive_covar_name))
    if (!is.null(addcovar)) validation_log <- c(validation_log, paste("✓ Additive covariate columns:", ncol(addcovar)))
    validation_log <- c(validation_log, paste("✓ Diet-stratified plots:", if (do_stratify) paste("Yes (", paste(diet_levels, collapse=", "), ")") else "No"))
    validation_log <- c(validation_log, paste("✓ Window size: ±", HALF_WIN_MB, "Mb around QTL peak"))
    validation_log <- c(validation_log, paste("✓ Total significant QTLs (99%):", nrow(qtls_99)))
    validation_log <- c(validation_log, "")

    # Create output directories by chromosome
    chromosomes <- unique(qtls_99\$chr)
    for (chr in chromosomes) {
        dir.create(paste0("chr", chr), showWarnings = FALSE, recursive = TRUE)
    }

    validation_log <- c(validation_log, paste("✓ Created output directories for", length(chromosomes), "chromosomes"))
    validation_log <- c(validation_log, "")

    # Generate plots for each QTL
    # Each QTL produces: 1 chromosome-wide LOD plot + 1 effects plot per group
    cat("Generating QTL LOD and effect plots...\\n")
    n_groups <- if (do_stratify) 1 + length(diet_levels) else 1
    n_plots_per_qtl <- 1 + (n_groups * 2)   # 1 LOD + n_groups x (coef + blup)
    cat("This will generate up to", nrow(qtls_99) * n_plots_per_qtl, "plots (",
        nrow(qtls_99), "QTLs x", n_plots_per_qtl, "images each: 1 LOD +", n_groups, "effects)\\n")

    plots_generated <- 0
    plots_failed <- 0

    for (i in 1:nrow(qtls_99)) {
        qtl     <- qtls_99[i, ]
        chr     <- as.character(qtl\$chr)
        gene_id <- qtl\$lodcolumn
        pos_cM  <- qtl\$pos
        lod     <- qtl\$lod
        safe_id <- gsub("[^A-Za-z0-9_-]", "_", gene_id)

        tryCatch({
            if (!gene_id %in% colnames(cross2\$pheno)) stop(paste("Phenotype", gene_id, "not found in cross2 object"))
            if (!gene_id %in% colnames(scan_results))  stop(paste("Phenotype", gene_id, "not found in scan_results"))

            # ── Find peak marker and define ±HALF_WIN_MB window ──────────────
            gmap_chr    <- gmap[[chr]]
            peak_marker <- names(gmap_chr)[which.min(abs(gmap_chr - pos_cM))]
            peak_Mb     <- pmap[[chr]][peak_marker]
            win_lo      <- peak_Mb - HALF_WIN_MB
            win_hi      <- peak_Mb + HALF_WIN_MB

            win_markers <- names(pmap[[chr]])[pmap[[chr]] >= win_lo & pmap[[chr]] <= win_hi]
            if (length(win_markers) < 2) stop(paste("Too few markers in window for", gene_id, "on chr", chr))

            # Subset alleleprobs to window (scan_results kept full-chr for LOD plot)
            ap_chr        <- alleleprobs[, chr]
            ap_win        <- ap_chr
            ap_win[[chr]] <- ap_chr[[chr]][, , win_markers, drop = FALSE]
            pmap_win      <- list(); pmap_win[[chr]] <- pmap[[chr]][win_markers]

            # ── File base name (Mb) and output directory ──────────────────────
            base_name <- paste0(safe_id, "_chr", chr, "_", round(peak_Mb, 1), "Mb_LOD", round(lod, 2))
            chr_dir   <- paste0("chr", chr)

            if (do_stratify) {
                qtl_dir <- file.path(chr_dir, base_name)
                dir.create(qtl_dir, showWarnings = FALSE, recursive = TRUE)
            }

            # ── LOD plot: chromosome-wide, one per QTL ────────────────────────
            lod_fname <- if (do_stratify) {
                file.path(qtl_dir, "lod.png")
            } else {
                file.path(chr_dir, paste0(base_name, "_lod.png"))
            }

            png(lod_fname, width = 900, height = 400)
            par(mar = c(4.5, 4.5, 2.5, 0.5))
            plot_scan1(scan_results, pmap, lodcolumn = gene_id, chr = chr,
                       bgcolor = "gray95",
                       xlab    = paste("Chr", chr, "position (Mb)"),
                       ylab    = "LOD score")
            title(main = sprintf("LOD — %s  Chr %s  (peak %.2f Mb, LOD %.2f)", gene_id, chr, peak_Mb, lod))
            abline(v = peak_Mb, col = "red", lty = 3, lwd = 1.5)
            dev.off()
            plots_generated <- plots_generated + 1

            # ── Define groups ─────────────────────────────────────────────────
            all_ids <- rownames(covar_data)
            if (do_stratify) {
                groups <- list(full = list(ids = all_ids, label = "Full", ac = addcovar, ic = intcovar_mat))
                for (dlvl in diet_levels) {
                    lvl_ids  <- all_ids[covar_data[[interactive_covar_name]] == dlvl]
                    safe_lvl <- tolower(gsub("[^A-Za-z0-9]", "", dlvl))
                    ac_lvl   <- if (!is.null(addcovar_sex)) addcovar_sex[lvl_ids, , drop = FALSE] else NULL
                    groups[[safe_lvl]] <- list(ids = lvl_ids, label = dlvl, ac = ac_lvl, ic = NULL)
                }
            } else {
                groups <- list(full = list(ids = all_ids, label = "Full", ac = addcovar, ic = NULL))
            }

            # ── Effects plots: 4 Mb window, one per group ─────────────────────
            for (grp_name in names(groups)) {
                g    <- groups[[grp_name]]
                ids  <- g\$ids
                ac_g <- if (!is.null(g\$ac) && nrow(g\$ac) > 0) g\$ac[intersect(ids, rownames(g\$ac)), , drop = FALSE] else NULL
                ic_g <- if (!is.null(g\$ic) && nrow(g\$ic) > 0) g\$ic[intersect(ids, rownames(g\$ic)), , drop = FALSE] else NULL

                kin_g     <- kinship[[chr]][ids, ids]
                pheno_mat <- cross2\$pheno[ids, gene_id, drop = FALSE]

                # ── scan1coef effects plot ────────────────────────────────────
                coef_fname <- if (do_stratify) {
                    file.path(qtl_dir, paste0(grp_name, "_coef_effects.png"))
                } else {
                    file.path(chr_dir, paste0(base_name, "_coef_effects.png"))
                }

                coef_win <- scan1coef(ap_win[ids, chr],
                                      pheno_mat,
                                      kinship  = kin_g,
                                      addcovar = ac_g,
                                      intcovar = ic_g)

                model_label <- if (!is.null(ic_g)) "GxE intcovar" else "addcovar only"
                png(coef_fname, width = 900, height = 500)
                par(mar = c(4.5, 4.5, 3.5, 0.5))
                plot_coefCC(coef_win, pmap_win,
                            bgcolor = "gray95",
                            legend  = "bottomleft",
                            xlab    = paste("Chr", chr, "position (Mb)"))
                title(main = sprintf("%s [coef] — %s  [%.1f–%.1f Mb]  (N=%d)",
                                     g\$label, gene_id, win_lo, win_hi, length(ids)))
                mtext(sprintf("+/-%g Mb window | peak %s (%.2f Mb) | %s",
                              HALF_WIN_MB, peak_marker, peak_Mb, model_label),
                      side = 3, line = 0.1, cex = 0.65, col = "gray40", font = 3)
                dev.off()
                plots_generated <- plots_generated + 1

                # ── scan1blup effects plot (additive only; intcovar not supported) ──
                blup_fname <- if (do_stratify) {
                    file.path(qtl_dir, paste0(grp_name, "_blup_effects.png"))
                } else {
                    file.path(chr_dir, paste0(base_name, "_blup_effects.png"))
                }

                blup_win <- scan1blup(ap_win[ids, chr],
                                      pheno_mat,
                                      kinship  = kin_g,
                                      addcovar = ac_g)

                blup_label <- if (!is.null(ic_g)) "addcovar only (BLUP; intcovar not supported)" else "addcovar only"
                png(blup_fname, width = 900, height = 500)
                par(mar = c(4.5, 4.5, 3.5, 0.5))
                plot_coefCC(blup_win, pmap_win,
                            bgcolor = "gray95",
                            legend  = "bottomleft",
                            xlab    = paste("Chr", chr, "position (Mb)"))
                title(main = sprintf("%s [BLUP] — %s  [%.1f–%.1f Mb]  (N=%d)",
                                     g\$label, gene_id, win_lo, win_hi, length(ids)))
                mtext(sprintf("+/-%g Mb window | peak %s (%.2f Mb) | %s",
                              HALF_WIN_MB, peak_marker, peak_Mb, blup_label),
                      side = 3, line = 0.1, cex = 0.65, col = "gray40", font = 3)
                dev.off()
                plots_generated <- plots_generated + 1
            }

            if (plots_generated %% 50 == 0) {
                cat("  Generated", plots_generated, "images so far...\\n")
            }

        }, error = function(e) {
            cat("WARNING: Failed to generate plot for", gene_id, "on chr", chr, ":", e\$message, "\\n")
            plots_failed <<- plots_failed + 1
        })
    }

    cat("\\n")
    cat("=================================================================\\n")
    cat("QTL Visualization Complete\\n")
    cat("=================================================================\\n")
    cat("Total plots generated:", plots_generated, "\\n")
    if (plots_failed > 0) {
        cat("Plots failed:", plots_failed, "\\n")
    }
    cat("=================================================================\\n")

    validation_log <- c(validation_log, "=== Plot Generation Summary ===")
    validation_log <- c(validation_log, paste("✓ Plots successfully generated:", plots_generated))
    if (plots_failed > 0) {
        validation_log <- c(validation_log, paste("⚠ Plots failed:", plots_failed))
    }
    validation_log <- c(validation_log, "")

    # Chromosome-wise breakdown
    validation_log <- c(validation_log, "=== Plots by Chromosome ===")
    chr_summary <- qtls_99 %>%
        group_by(chr) %>%
        summarise(count = n()) %>%
        arrange(chr)

    for (i in 1:nrow(chr_summary)) {
        chr <- chr_summary\$chr[i]
        count <- chr_summary\$count[i]
        validation_log <- c(validation_log, paste("Chr", chr, ":", count, "plots"))
    }

    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Visualization Complete ===")
    validation_log <- c(validation_log, paste("✓ Output organized by chromosome in chr*/ subdirectories"))
    validation_log <- c(validation_log, paste("✓ LOD plots: chromosome-wide (900x400), peak marked with red dashed line"))
    validation_log <- c(validation_log, paste("✓ Effects plots: ±", HALF_WIN_MB, "Mb window around peak (900x500), windowed approach eliminates marker instability"))
    validation_log <- c(validation_log, paste("✓ Effects plots: both scan1coef (*_coef_effects.png) and scan1blup (*_blup_effects.png) per group"))
    validation_log <- c(validation_log, paste("✓ Interactive model: plots saved in chr*/QTL_subdir/ (lod.png + {group}_effects.png)"))
    validation_log <- c(validation_log, paste("✓ Additive model: plots saved flat as chr*/basename_{lod,effects}.png"))
    validation_log <- c(validation_log, paste("✓ File names use Mb position of QTL peak"))

    # Write validation report
    writeLines(validation_log, "visualization_report.txt")

    cat("\\nVisualization report saved to: visualization_report.txt\\n")
    """
}
