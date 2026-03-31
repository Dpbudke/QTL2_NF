#!/usr/bin/env nextflow

// Broman Sample Mixup QC Module
// Detects sample mix-ups in Diversity Outbred mouse data using expression vs genotypes
// Based on methodology from Broman et al. (2015) G3 and Westra et al. (2011)
// https://kbroman.org/qtl2/assets/vignettes/do_mixups.html

process BROMAN_MIXUP_QC {
    tag "Broman mixup QC for ${study_prefix}"
    publishDir "${params.outdir}/mixup_qc", mode: 'copy'

    cpus 8
    memory '64 GB'
    time '4h'

    input:
    path(cross2_file)
    path(genoprob_file)
    path(expr_data_file)       // Expression data (RDS matrix OR CSV file with samples x genes)
    path(expr_peaks_file)      // QTL peaks CSV (from module 8 significant_qtls output)
    val(study_prefix)
    val(n_top_eqtl)            // Number of top eQTL to use (default: 100)

    output:
    path("${study_prefix}_mixup_report.html"), emit: mixup_report
    path("${study_prefix}_mixup_problems.csv"), emit: mixup_problems
    path("${study_prefix}_distance_matrix.csv"), emit: distances
    path("${study_prefix}_mixup_plots.jpg"), emit: plots
    path("mixup_qc_log.txt"), emit: qc_log

    script:
    """
    #!/usr/bin/env Rscript

    # Set up writable R library path in work directory
    work_lib <- file.path(getwd(), "R_libs")
    dir.create(work_lib, showWarnings = FALSE, recursive = TRUE)
    .libPaths(c(work_lib, .libPaths()))

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)

        # Install lineup2 from CRAN if not available
        if (!require("lineup2", quietly = TRUE)) {
            message("Installing lineup2 package to work directory...")
            install.packages("lineup2", repos = "https://cloud.r-project.org",
                           lib = work_lib, dependencies = TRUE)
            library(lineup2)
        }
    })

    # Initialize logging
    log_file <- "mixup_qc_log.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message("=== BROMAN SAMPLE MIXUP QC ANALYSIS ===")
    log_message("Study: ${study_prefix}")
    log_message("Top eQTL to analyze: ${n_top_eqtl}")

    # Load data
    log_message("Loading cross2 object...")
    cross2 <- readRDS("${cross2_file}")

    log_message("Loading genotype probabilities...")
    genoprobs <- readRDS("${genoprob_file}")

    log_message("Loading expression data...")
    # Handle both RDS and CSV input formats
    expr_file <- "${expr_data_file}"
    if (grepl("\\\\.rds\$", expr_file, ignore.case=TRUE)) {
        log_message("  Format: RDS file")
        expr_data <- readRDS(expr_file)
    } else if (grepl("\\\\.csv\$", expr_file, ignore.case=TRUE)) {
        log_message("  Format: CSV file")
        expr_df <- read.csv(expr_file, row.names=1, check.names=FALSE)
        expr_data <- as.matrix(expr_df)
        log_message("  Converted CSV to matrix")
    } else {
        stop("Expression data file must be either .rds or .csv format")
    }

    log_message("Loading eQTL peaks...")
    eqtl_peaks <- read.csv("${expr_peaks_file}", stringsAsFactors = FALSE)

    # Validate data
    log_message(paste("Cross2 samples:", nrow(cross2\$pheno)))
    log_message(paste("Expression data dimensions:", nrow(expr_data), "x", ncol(expr_data)))
    log_message(paste("Total eQTL peaks:", nrow(eqtl_peaks)))

    # Select top N eQTL by LOD score, keeping only unique genes
    n_traits <- as.integer(${n_top_eqtl})

    # Sort by LOD score and keep only the top eQTL per gene (unique lodcolumn)
    eqtl_peaks_sorted <- eqtl_peaks[order(eqtl_peaks\$lod, decreasing=TRUE), ]
    eqtl_peaks_unique <- eqtl_peaks_sorted[!duplicated(eqtl_peaks_sorted\$lodcolumn), ]

    log_message(paste("Total peaks:", nrow(eqtl_peaks), "| Unique genes:", nrow(eqtl_peaks_unique)))

    if (nrow(eqtl_peaks_unique) < n_traits) {
        log_message(paste("WARNING: Only", nrow(eqtl_peaks_unique), "unique genes available, using all"))
        n_traits <- nrow(eqtl_peaks_unique)
    }

    top_peaks <- eqtl_peaks_unique[1:n_traits, ]
    log_message(paste("Selected top", n_traits, "unique eQTL for analysis"))
    log_message(paste("LOD range:", round(min(top_peaks\$lod), 2), "to", round(max(top_peaks\$lod), 2)))

    # Extract observed expression for top eQTL genes
    log_message("Extracting observed expression values...")
    gene_ids <- top_peaks\$lodcolumn
    obs_expr <- expr_data[, gene_ids, drop=FALSE]
    log_message(paste("Observed expression matrix:", nrow(obs_expr), "samples x", ncol(obs_expr), "genes"))

    # Calculate predicted expression from genotype probabilities
    log_message("Calculating predicted expression from genotypes...")
    exp_expr <- NULL

    # Extract map from cross2 object for marker lookup
    map <- cross2\$gmap
    if (is.null(map)) {
        log_message("WARNING: gmap not found, trying pmap")
        map <- cross2\$pmap
    }

    for(i in 1:nrow(top_peaks)) {
        if (i %% 10 == 0) {
            log_message(paste("  Processing gene", i, "of", nrow(top_peaks)))
        }

        # Pull genotype probabilities at eQTL marker
        pr <- pull_genoprobpos(genoprobs, map, top_peaks\$chr[i], top_peaks\$pos[i])

        # Fit single-QTL model (zerosum=FALSE keeps original scale)
        out_fit1 <- fit1(pr, obs_expr[,i,drop=FALSE], zerosum=FALSE)

        # Calculate predicted expression
        fitted <- pr %*% out_fit1\$coef

        # Initialize matrix on first iteration
        if(is.null(exp_expr)) {
            exp_expr <- matrix(nrow=nrow(fitted), ncol=ncol(obs_expr))
            dimnames(exp_expr) <- list(rownames(pr), colnames(obs_expr))
        }

        exp_expr[rownames(fitted), i] <- fitted
    }

    log_message(paste("Predicted expression matrix:", nrow(exp_expr), "samples x", ncol(exp_expr), "genes"))

    # Calculate RMS distances
    log_message("Calculating RMS distances...")
    d_evg <- dist_betw_matrices(obs_expr, exp_expr)

    # Extract self and minimum distances
    self_dist <- get_self(d_evg)
    min_dist <- get_best(d_evg)

    # Identify problems
    problems <- self_dist > min_dist
    n_problems <- sum(problems)

    log_message(paste("Samples with potential mix-ups:", n_problems))

    # Create problem report
    problem_df <- data.frame(
        sample = names(self_dist),
        self_distance = self_dist,
        min_distance = min_dist,
        best_match = colnames(d_evg)[apply(d_evg, 1, which.min)],
        is_problem = problems
    )

    write.csv(problem_df, "${study_prefix}_mixup_problems.csv", row.names=FALSE)
    write.csv(d_evg, "${study_prefix}_distance_matrix.csv")

    # Create plots
    jpeg("${study_prefix}_mixup_plots.jpg", width=10, height=8, units="in", res=300)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    plot(self_dist, min_dist,
         xlab="Self distance", ylab="Minimum distance",
         main=paste("Mixup QC:", "${study_prefix}"),
         pch=19, col=ifelse(problems, "red", "gray60"))
    abline(0, 1, lty=2, col="blue")
    dev.off()

    # Create HTML report
    html <- paste0("<html><body><h1>Broman Mixup QC - ${study_prefix}</h1>",
                   "<p>Samples with problems: ", n_problems, "</p></body></html>")
    writeLines(html, "${study_prefix}_mixup_report.html")

    log_message("=== ANALYSIS COMPLETE ===")
    close(log_conn)
    """
}
