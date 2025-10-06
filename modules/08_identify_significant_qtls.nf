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
    cross2 <- readRDS("${cross2_file}")  # This is the filtered cross2 from module 7
    scan_results_full <- readRDS("${scan_results_file}")
    perm_results <- readRDS("${perm_results_file}")
    permThresh <- readRDS("${thresholds_file}")

    # Filter scan_results to match the phenotypes in filtered cross2
    filtered_phenos <- colnames(cross2\$pheno)
    scan_results <- scan_results_full[, filtered_phenos, drop=FALSE]

    # Align the order of permThresh columns with scan_results columns
    permThresh <- permThresh[, match(colnames(scan_results), colnames(permThresh))]

    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Filtered cross2 loaded:", ncol(cross2\$pheno), "phenotypes,", nrow(cross2\$pheno), "individuals"))
    validation_log <- c(validation_log, paste("✓ Scan results filtered:", paste(dim(scan_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Permutation results loaded:", paste(dim(perm_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Significance thresholds loaded:", paste(dim(permThresh), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Column alignment verified:", sum(colnames(permThresh) == colnames(scan_results)), "of", ncol(scan_results), "phenotypes"))
    validation_log <- c(validation_log, "")

    # Initialize results data frame
    all_significant_qtls <- data.frame()
    qtl_counts <- data.frame(
        significance_level = rownames(permThresh),
        qtl_count = 0,
        phenotypes_with_qtls = 0,
        stringsAsFactors = FALSE
    )

    cat("Identifying significant QTLs for each significance level...\\n")

    # Process each significance level (each row of permThresh)
    for (i in 1:nrow(permThresh)) {
        sig_level <- rownames(permThresh)[i]
        cat("Processing", sig_level, "significance level...\\n")

        # Get threshold vector for this significance level
        thresh <- permThresh[i, ]

        # Find peaks above this threshold
        peaks <- find_peaks(scan_results, cross2\$gmap, threshold=thresh, drop=1.5)

        if (nrow(peaks) > 0) {
            # Annotate each peak with ALL four threshold values for that phenotype
            for (j in 1:nrow(permThresh)) {
                desig <- paste0("signif_", gsub("%", "", rownames(permThresh)[j]))
                peaks[[desig]] <- permThresh[j, ][match(peaks\$lodcolumn, colnames(permThresh))]
            }

            # Add study and significance level metadata
            peaks\$study <- "${prefix}"
            peaks\$significance_level <- sig_level

            all_significant_qtls <- rbind(all_significant_qtls, peaks)

            # Count unique phenotypes at this level
            phenos_at_level <- length(unique(peaks\$lodcolumn))
            qtl_counts[qtl_counts\$significance_level == sig_level, "qtl_count"] <- nrow(peaks)
            qtl_counts[qtl_counts\$significance_level == sig_level, "phenotypes_with_qtls"] <- phenos_at_level

            validation_log <- c(validation_log, paste("✓", sig_level, "level:", nrow(peaks), "QTLs in", phenos_at_level, "phenotypes"))
        } else {
            validation_log <- c(validation_log, paste("✓", sig_level, "level: 0 QTLs"))
        }
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