process VISUALIZE_QTLS {
    tag "Generating QTL visualizations for ${prefix}"
    publishDir "${params.outdir}/10_visualize", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '8h'

    input:
    path(alleleprob_file)
    path(genetic_map_file)
    path(scan_results_file)
    path(cross2_file)
    path(significant_qtls_file)
    val(prefix)

    output:
    path("chr*/*.png"), emit: qtl_plots, optional: true
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
    gmap <- readRDS("${genetic_map_file}")
    scan_results <- readRDS("${scan_results_file}")
    cross2 <- readRDS("${cross2_file}")

    # Load significant QTLs and filter to 99% significance level
    significant_qtls <- read.csv("${significant_qtls_file}", stringsAsFactors = FALSE)
    qtls_99 <- significant_qtls %>% filter(significance_level == "99%")

    validation_log <- c(validation_log, "=== Data Loading Complete ===")
    validation_log <- c(validation_log, paste("✓ Alleleprobs loaded:", length(alleleprobs), "chromosomes"))
    validation_log <- c(validation_log, paste("✓ Genetic map loaded:", length(gmap), "chromosomes"))
    validation_log <- c(validation_log, paste("✓ Scan results loaded:", paste(dim(scan_results), collapse=" x ")))
    validation_log <- c(validation_log, paste("✓ Cross2 object loaded:", nrow(cross2\$pheno), "individuals,", ncol(cross2\$pheno), "phenotypes"))
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
    cat("Generating QTL coefficient plots...\\n")
    cat("This will generate", nrow(qtls_99), "plots (one per significant QTL at 99% threshold)\\n")

    plots_generated <- 0
    plots_failed <- 0

    for (i in 1:nrow(qtls_99)) {
        qtl <- qtls_99[i, ]
        chr <- as.character(qtl\$chr)
        gene_id <- qtl\$lodcolumn
        pos <- qtl\$pos
        lod <- qtl\$lod

        # Create safe filename (remove any special characters)
        safe_gene_id <- gsub("[^A-Za-z0-9_-]", "_", gene_id)
        filename <- paste0("chr", chr, "/", safe_gene_id, "_chr", chr, "_",
                          round(pos, 1), "cM_LOD", round(lod, 2), ".png")

        tryCatch({
            # Check if phenotype exists in cross2 and scan_results
            if (!gene_id %in% colnames(cross2\$pheno)) {
                stop(paste("Phenotype", gene_id, "not found in cross2 object"))
            }
            if (!gene_id %in% colnames(scan_results)) {
                stop(paste("Phenotype", gene_id, "not found in scan_results"))
            }

            # Calculate coefficients for this chromosome and phenotype
            coef_chr <- scan1coef(alleleprobs[, chr], cross2\$pheno[, gene_id, drop=FALSE])

            # Extract chromosome-specific genetic map
            gmap_chr <- gmap[chr]

            # Extract scan results for this phenotype (all chromosomes for context)
            # but we'll focus on the specific chromosome in the plot
            scan_pheno <- scan_results[, gene_id, drop=FALSE]

            # Open PNG device
            png(filename, width=800, height=600)
            par(mar=c(4.1, 4.1, 0.6, 0.6))

            # Generate plot with allele effects and LOD score
            plot_coefCC(coef_chr, gmap_chr, scan1_output=scan_pheno,
                       bgcolor="gray95", legend="bottomleft")

            # Close device
            dev.off()

            plots_generated <- plots_generated + 1

            # Progress update every 100 plots
            if (plots_generated %% 100 == 0) {
                cat("  Generated", plots_generated, "of", nrow(qtls_99), "plots...\\n")
            }

        }, error = function(e) {
            cat("WARNING: Failed to generate plot for", gene_id, "on chr", chr, ":", e\$message, "\\n")
            plots_failed <- plots_failed + 1
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
    validation_log <- c(validation_log, paste("✓ Each plot shows founder allele effects and LOD scores"))
    validation_log <- c(validation_log, paste("✓ Plot dimensions: 800x600 pixels with mar=c(4.1, 4.1, 0.6, 0.6)"))

    # Write validation report
    writeLines(validation_log, "visualization_report.txt")

    cat("\\nVisualization report saved to: visualization_report.txt\\n")
    """
}
