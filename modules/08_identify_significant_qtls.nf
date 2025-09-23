process IDENTIFY_SIGNIFICANT_QTLS {
    tag "Identifying significant QTLs for ${prefix}"
    publishDir "${params.outdir}/08_significant_qtls", mode: 'copy'

    input:
    path(cross2_file)
    path(scan_results_file)
    path(perm_results_file)
    path(thresholds_file)
    val(prefix)

    output:
    path("${prefix}_significant_qtls.csv"), emit: significant_qtls
    path("${prefix}_qtl_summary.txt"), emit: qtl_summary
    path("qtl_identification_report.txt"), emit: validation_report

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
    validation_log <- c(validation_log, paste("=== QTL Identification Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")

    cat("Loading scan results and significance thresholds...\\n")

    # Load required data
    cross2 <- readRDS("${cross2_file}")
    scan_results <- readRDS("${scan_results_file}")
    perm_results <- readRDS("${perm_results_file}")
    thresholds <- read.csv("${thresholds_file}", stringsAsFactors = FALSE)

    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Scan results loaded:", paste(dim(scan_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Permutation results loaded:", paste(dim(perm_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Significance thresholds loaded:", nrow(thresholds), "entries"))
    validation_log <- c(validation_log, "")

    # Initialize results data frame
    all_significant_qtls <- data.frame()
    qtl_counts <- data.frame(
        significance_level = c("suggestive", "0.10", "0.05", "0.01"),
        qtl_count = 0,
        phenotypes_with_qtls = 0,
        stringsAsFactors = FALSE
    )

    cat("Identifying significant QTLs for each significance level...\\n")

    # Process each significance level
    for (sig_level in c("suggestive", "0.10", "0.05", "0.01")) {
        cat("Processing", sig_level, "significance level...\\n")

        # Get thresholds for this significance level
        level_thresholds <- thresholds[thresholds\$significance_level == sig_level, ]

        total_qtls <- 0
        phenotypes_with_qtls <- 0

        # Process each phenotype
        for (i in 1:nrow(level_thresholds)) {
            pheno_name <- level_thresholds\$phenotype[i]
            threshold <- level_thresholds\$lod_threshold[i]

            # Find peaks for this phenotype at this significance level
            if (pheno_name %in% colnames(scan_results)) {
                # Create single-phenotype scan result
                pheno_scan <- scan_results[, pheno_name, drop=FALSE]

                # Find peaks above threshold
                peaks <- find_peaks(pheno_scan, cross2\$gmap, threshold=threshold, drop=1.5)

                if (nrow(peaks) > 0) {
                    # Add metadata to peaks
                    peaks\$study <- "${prefix}"
                    peaks\$phenotype <- pheno_name
                    peaks\$significance_level <- sig_level
                    peaks\$alpha <- level_thresholds\$alpha[i]
                    peaks\$threshold_used <- threshold

                    # Reorder columns
                    peaks <- peaks[, c("study", "phenotype", "significance_level", "alpha",
                                      "threshold_used", "lodindex", "lodcolumn", "chr", "pos", "lod")]

                    all_significant_qtls <- rbind(all_significant_qtls, peaks)
                    total_qtls <- total_qtls + nrow(peaks)
                    phenotypes_with_qtls <- phenotypes_with_qtls + 1
                }
            }
        }

        # Update counts
        qtl_counts[qtl_counts\$significance_level == sig_level, "qtl_count"] <- total_qtls
        qtl_counts[qtl_counts\$significance_level == sig_level, "phenotypes_with_qtls"] <- phenotypes_with_qtls

        validation_log <- c(validation_log, paste("✓", sig_level, "level:", total_qtls, "QTLs in", phenotypes_with_qtls, "phenotypes"))
    }

    # Write significant QTLs to file
    if (nrow(all_significant_qtls) > 0) {
        write.csv(all_significant_qtls, file = "${prefix}_significant_qtls.csv", row.names = FALSE)
        validation_log <- c(validation_log, paste("✓ Total significant QTLs identified:", nrow(all_significant_qtls)))
    } else {
        # Create empty file
        empty_qtls <- data.frame(
            study=character(0), phenotype=character(0), significance_level=character(0),
            alpha=numeric(0), threshold_used=numeric(0), lodindex=numeric(0),
            lodcolumn=character(0), chr=character(0), pos=numeric(0), lod=numeric(0)
        )
        write.csv(empty_qtls, file = "${prefix}_significant_qtls.csv", row.names = FALSE)
        validation_log <- c(validation_log, "⚠ No significant QTLs identified at any significance level")
    }

    # Create summary report
    cat("Generating QTL summary report...\\n")

    summary_lines <- c()
    summary_lines <- c(summary_lines, paste("QTL Mapping Summary for Study:", "${prefix}"))
    summary_lines <- c(summary_lines, paste("Generated:", Sys.time()))
    summary_lines <- c(summary_lines, "")
    summary_lines <- c(summary_lines, "=== Dataset Overview ===")
    summary_lines <- c(summary_lines, paste("Total individuals:", nrow(cross2\$pheno)))
    summary_lines <- c(summary_lines, paste("Total phenotypes:", ncol(cross2\$pheno)))
    summary_lines <- c(summary_lines, paste("Chromosomes analyzed:", paste(names(cross2\$gmap), collapse=", ")))
    summary_lines <- c(summary_lines, paste("Total markers scanned:", nrow(scan_results)))
    summary_lines <- c(summary_lines, "")
    summary_lines <- c(summary_lines, "=== QTL Summary by Significance Level ===")

    for (i in 1:nrow(qtl_counts)) {
        level <- qtl_counts\$significance_level[i]
        count <- qtl_counts\$qtl_count[i]
        phenos <- qtl_counts\$phenotypes_with_qtls[i]
        summary_lines <- c(summary_lines, paste(level, ":", count, "QTLs in", phenos, "phenotypes"))
    }

    if (nrow(all_significant_qtls) > 0) {
        summary_lines <- c(summary_lines, "")
        summary_lines <- c(summary_lines, "=== Chromosome Distribution ===")
        chr_counts <- table(all_significant_qtls\$chr)
        for (chr in names(chr_counts)) {
            summary_lines <- c(summary_lines, paste("Chr", chr, ":", chr_counts[chr], "QTLs"))
        }

        summary_lines <- c(summary_lines, "")
        summary_lines <- c(summary_lines, "=== Top QTLs (LOD ≥ 10) ===")
        top_qtls <- all_significant_qtls[all_significant_qtls\$lod >= 10, ]
        if (nrow(top_qtls) > 0) {
            top_qtls <- top_qtls[order(top_qtls\$lod, decreasing=TRUE), ]
            for (i in 1:min(10, nrow(top_qtls))) {
                qtl <- top_qtls[i, ]
                summary_lines <- c(summary_lines,
                    paste("  ", qtl\$phenotype, "- Chr", qtl\$chr, "@", round(qtl\$pos, 1), "cM, LOD =", round(qtl\$lod, 2)))
            }
        } else {
            summary_lines <- c(summary_lines, "  No QTLs with LOD ≥ 10")
        }
    }

    # Write summary
    writeLines(summary_lines, "${prefix}_qtl_summary.txt")

    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== QTL Identification Complete ===")
    validation_log <- c(validation_log, paste("✓ Significant QTLs file:", "${prefix}_significant_qtls.csv"))
    validation_log <- c(validation_log, paste("✓ Summary report:", "${prefix}_qtl_summary.txt"))
    validation_log <- c(validation_log, "✓ Ready for results visualization and downstream analysis")

    # Write validation report
    writeLines(validation_log, "qtl_identification_report.txt")

    cat("QTL identification completed successfully\\n")
    """
}