// ─── MODULE 8: Significant QTL Identification (scatter-gather) ───────────────
//
// Architecture:
//   FIND_PEAKS_CHR  — loads full scan matrix, subsets to one chromosome,
//                     runs find_peaks (20 SLURM jobs submitted in parallel
//                     via Channel.of() in the workflow)
//   GATHER_PEAKS    — collects all 20 chr peak files, deduplicates to best
//                     peak per gene, annotates all 4 thresholds, writes outputs
//
// Scatter is driven by Channel.of("1","2",...,"19","X") — no glob/flatten/join
// needed, which avoids the resume-mode caching ambiguity of the previous
// SPLIT_SCAN_BY_CHR → flatten → join pattern.

process FIND_PEAKS_CHR {
    tag "find_peaks chr ${chr}"

    input:
    val(chr)
    path(scan_results_file)
    path(cross2_file)
    path(thresholds_file)

    output:
    path("peaks_chr_${chr}.rds"), emit: chr_peaks

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(qtl2))

    cat("Chr ${chr}: loading scan results...\\n")
    scan_results <- readRDS("${scan_results_file}")
    cross2       <- readRDS("${cross2_file}")
    permThresh   <- readRDS("${thresholds_file}")

    # Restrict to phenotypes in the filtered cross2
    filtered_phenos <- colnames(cross2\$pheno)
    scan_results    <- scan_results[, intersect(filtered_phenos, colnames(scan_results)), drop=FALSE]
    class(scan_results) <- "scan1"

    # Subset to this chromosome's markers
    chr_markers     <- names(cross2\$gmap[["${chr}"]])
    markers_in_scan <- chr_markers[chr_markers %in% rownames(scan_results)]

    if (length(markers_in_scan) == 0) {
        cat("Chr ${chr}: no markers in scan results — writing empty peaks file\\n")
        saveRDS(data.frame(), "peaks_chr_${chr}.rds")
        quit(save="no", status=0)
    }

    chr_scan        <- scan_results[markers_in_scan, , drop=FALSE]
    class(chr_scan) <- "scan1"
    chr_gmap        <- cross2\$gmap["${chr}"]

    # Per-phenotype 63% threshold
    thresh_63        <- permThresh[1, match(colnames(chr_scan), colnames(permThresh))]
    names(thresh_63) <- colnames(chr_scan)

    peaks <- find_peaks(chr_scan, chr_gmap, threshold=thresh_63, drop=1.5, cores=1)
    cat("Chr ${chr}:", nrow(peaks), "peaks above 63% threshold\\n")

    saveRDS(peaks, "peaks_chr_${chr}.rds")
    """
}

process GATHER_PEAKS {
    tag "Gathering peaks for ${prefix}"
    publishDir "${params.outdir}/08_significant_qtls", mode: 'copy'

    input:
    path(chr_peaks_files)   // collected from FIND_PEAKS_CHR — all 20 files staged here
    path(cross2_file)
    path(perm_results_file)
    path(thresholds_file)
    val(prefix)

    output:
    path("${prefix}_significant_qtls.csv"), emit: significant_qtls
    path("${prefix}_qtl_summary.txt"),      emit: qtl_summary
    path("qtl_identification_report.txt"),   emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(qtl2)
        library(dplyr)
    })

    validation_log <- c(
        "=== QTL Identification Report ===",
        paste("Timestamp:", Sys.time()),
        paste("Study Prefix:", "${prefix}"),
        ""
    )

    cat("Loading per-chromosome peak files...\\n")
    peak_files     <- list.files(".", pattern="^peaks_chr_")
    all_peaks_list <- lapply(peak_files, readRDS)
    non_empty      <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, all_peaks_list)
    all_peaks      <- if (length(non_empty) > 0) do.call(rbind, non_empty) else data.frame()

    cat("Chromosomes processed:", length(peak_files), "— raw peaks:", nrow(all_peaks), "\\n")

    cross2       <- readRDS("${cross2_file}")
    perm_results <- readRDS("${perm_results_file}")
    permThresh   <- readRDS("${thresholds_file}")

    filtered_phenos <- colnames(cross2\$pheno)
    permThresh      <- permThresh[, match(filtered_phenos, colnames(permThresh)), drop=FALSE]

    validation_log <- c(validation_log,
        "=== Data Loading Complete ===",
        paste("✓ Filtered cross2:", ncol(cross2\$pheno), "phenotypes,", nrow(cross2\$pheno), "individuals"),
        paste("✓ Permutation results loaded:", paste(dim(perm_results), collapse=" x ")),
        paste("✓ Significance thresholds:", paste(dim(permThresh), collapse=" x ")),
        paste("✓ Chromosomes with peaks files:", length(peak_files)),
        paste("✓ Raw peaks across all chromosomes:", nrow(all_peaks)),
        ""
    )

    qtl_counts <- data.frame(
        significance_level   = rownames(permThresh),
        qtl_count            = 0L,
        phenotypes_with_qtls = 0L,
        stringsAsFactors     = FALSE
    )

    all_significant_qtls <- data.frame()

    if (nrow(all_peaks) > 0) {
        cat("Total peaks (one per chromosome per gene):", nrow(all_peaks), "\\n")

        # Annotate each peak with all four threshold values for that gene
        for (j in 1:nrow(permThresh)) {
            desig              <- paste0("signif_", gsub("%", "", rownames(permThresh)[j]))
            all_peaks[[desig]] <- permThresh[j, ][match(all_peaks\$lodcolumn, colnames(permThresh))]
        }

        all_peaks\$study <- "${prefix}"

        # Emit one row per (peak × qualifying threshold), matching pre-refactor behavior
        for (i in 1:nrow(permThresh)) {
            sig_level  <- rownames(permThresh)[i]
            thresh_col <- paste0("signif_", gsub("%", "", sig_level))
            peaks_at   <- all_peaks[all_peaks\$lod >= all_peaks[[thresh_col]], , drop = FALSE]
            if (nrow(peaks_at) > 0) {
                peaks_at\$significance_level <- sig_level
                all_significant_qtls <- rbind(all_significant_qtls, peaks_at)
            }
            phenos_at <- length(unique(peaks_at\$lodcolumn))
            qtl_counts[qtl_counts\$significance_level == sig_level, "qtl_count"]            <- nrow(peaks_at)
            qtl_counts[qtl_counts\$significance_level == sig_level, "phenotypes_with_qtls"] <- phenos_at
            validation_log <- c(validation_log,
                paste("✓", sig_level, "level:", nrow(peaks_at), "QTLs in", phenos_at, "phenotypes"))
        }
    } else {
        validation_log <- c(validation_log, "✓ No QTLs above 63% threshold found")
    }

    if (nrow(all_significant_qtls) > 0) {
        write.csv(all_significant_qtls, file="${prefix}_significant_qtls.csv", row.names=FALSE)
        validation_log <- c(validation_log,
            paste("✓ Total significant QTLs identified:", nrow(all_significant_qtls)))
    } else {
        empty_qtls <- data.frame(
            study=character(0), phenotype=character(0), significance_level=character(0),
            alpha=numeric(0), threshold_used=numeric(0), lodindex=numeric(0),
            lodcolumn=character(0), chr=character(0), pos=numeric(0), lod=numeric(0)
        )
        write.csv(empty_qtls, file="${prefix}_significant_qtls.csv", row.names=FALSE)
        validation_log <- c(validation_log, "⚠ No significant QTLs identified at any significance level")
    }

    cat("Generating QTL summary report...\\n")
    summary_lines <- c(
        paste("QTL Mapping Summary for Study:", "${prefix}"),
        paste("Generated:", Sys.time()),
        "",
        "=== Dataset Overview ===",
        paste("Total individuals:", nrow(cross2\$pheno)),
        paste("Total phenotypes:", ncol(cross2\$pheno)),
        paste("Chromosomes analyzed:", paste(names(cross2\$gmap), collapse=", ")),
        "",
        "=== QTL Summary by Significance Level ==="
    )

    for (i in 1:nrow(qtl_counts)) {
        summary_lines <- c(summary_lines,
            paste(qtl_counts\$significance_level[i], ":",
                  qtl_counts\$qtl_count[i], "QTLs in",
                  qtl_counts\$phenotypes_with_qtls[i], "phenotypes"))
    }

    if (nrow(all_significant_qtls) > 0) {
        summary_lines <- c(summary_lines, "", "=== Chromosome Distribution ===")
        chr_counts <- table(all_significant_qtls\$chr)
        for (chr in names(chr_counts)) {
            summary_lines <- c(summary_lines, paste("Chr", chr, ":", chr_counts[chr], "QTLs"))
        }

        summary_lines <- c(summary_lines, "", "=== Top QTLs (LOD >= 10) ===")
        top_qtls <- all_significant_qtls[all_significant_qtls\$lod >= 10, ]
        if (nrow(top_qtls) > 0) {
            top_qtls <- top_qtls[order(top_qtls\$lod, decreasing=TRUE), ]
            for (i in 1:min(10, nrow(top_qtls))) {
                q <- top_qtls[i, ]
                summary_lines <- c(summary_lines,
                    paste("  ", q\$lodcolumn, "- Chr", q\$chr, "@", round(q\$pos, 1), "cM, LOD =", round(q\$lod, 2)))
            }
        } else {
            summary_lines <- c(summary_lines, "  No QTLs with LOD >= 10")
        }
    }

    writeLines(summary_lines, "${prefix}_qtl_summary.txt")

    validation_log <- c(validation_log, "",
        "=== QTL Identification Complete ===",
        paste("✓ Significant QTLs file:", "${prefix}_significant_qtls.csv"),
        paste("✓ Summary report:", "${prefix}_qtl_summary.txt"),
        "✓ Ready for results visualization and downstream analysis"
    )
    writeLines(validation_log, "qtl_identification_report.txt")

    cat("QTL identification completed successfully\\n")
    """
}
