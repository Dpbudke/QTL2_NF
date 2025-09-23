process PREPARE_GENOME_SCAN {
    tag "Preparing genome scan for ${prefix}"
    publishDir "${params.outdir}/05_genome_scan_preparation", mode: 'copy'

    // Increased resources for better performance
    cpus 32
    memory '128 GB'
    time '8h'

    input:
    path(cross2_file)
    val(prefix)

    output:
    path("${prefix}_genoprob.rds"), emit: genoprob
    path("${prefix}_kinship_loco.rds"), emit: kinship
    path("${prefix}_genetic_map.rds"), emit: genetic_map
    path("genoprob_validation_report.txt"), emit: validation_report

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })

    # Initialize validation report
    validation_log <- c()
    validation_log <- c(validation_log, paste("=== Genome Scan Preparation Report ==="))
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")

    cat("Loading cross2 object...\\n")
    cross2 <- readRDS("${cross2_file}")

    validation_log <- c(validation_log, "=== Cross2 Object Loaded ===")
    validation_log <- c(validation_log, paste("✓ Number of individuals:", nrow(cross2\$pheno)))
    validation_log <- c(validation_log, paste("✓ Number of chromosomes:", length(cross2\$gmap)))
    validation_log <- c(validation_log, paste("✓ Total markers:", sum(sapply(cross2\$gmap, length))))
    validation_log <- c(validation_log, paste("✓ Number of phenotypes:", ncol(cross2\$pheno)))
    validation_log <- c(validation_log, "")

    # Get genetic map and insert pseudomarkers for high-resolution scanning
    cat("Preparing genetic map with pseudomarkers...\\n")
    gmap <- cross2\$gmap

    # Debug: check cross2 and genetic map structure
    cat("DEBUG: Cross2 class:", class(cross2), "\\n")
    cat("DEBUG: Cross2 names:", paste(names(cross2), collapse=", "), "\\n")
    cat("DEBUG: Genotypes class:", class(cross2\$geno), "\\n")
    cat("DEBUG: Genotypes names:", paste(names(cross2\$geno), collapse=", "), "\\n")
    if (length(cross2\$geno) > 0) {
        cat("DEBUG: First geno chr class:", class(cross2\$geno[[1]]), "\\n")
        cat("DEBUG: First geno chr dim:", paste(dim(cross2\$geno[[1]]), collapse=" x "), "\\n")
    }
    cat("DEBUG: Genetic map class:", class(gmap), "\\n")
    cat("DEBUG: Genetic map length:", length(gmap), "\\n")
    if (length(gmap) > 0) {
        cat("DEBUG: First chromosome name:", names(gmap)[1], "\\n")
        cat("DEBUG: First chromosome class:", class(gmap[[1]]), "\\n")
        cat("DEBUG: First chromosome length:", length(gmap[[1]]), "\\n")
        if (length(gmap[[1]]) > 0) {
            cat("DEBUG: First few markers:", paste(head(names(gmap[[1]]), 5), collapse=", "), "\\n")
        }
    }

    # Insert pseudomarkers every 0.5 cM for high resolution
    gmap_pseudo <- insert_pseudomarkers(gmap, step=0.5)

    # Save genetic map
    saveRDS(gmap_pseudo, file = "${prefix}_genetic_map.rds")

    validation_log <- c(validation_log, "=== Genetic Map Preparation ===")
    validation_log <- c(validation_log, paste("✓ Original markers:", sum(sapply(gmap, length))))
    validation_log <- c(validation_log, paste("✓ Markers with pseudomarkers (0.5 cM):", sum(sapply(gmap_pseudo, length))))
    validation_log <- c(validation_log, "")

    # Calculate genotype probabilities
    cat("Calculating genotype probabilities...\\n")
    cat("This may take several minutes for large datasets...\\n")

    validation_log <- c(validation_log, "=== Genotype Probability Calculation ===")
    start_time <- Sys.time()

    # Use error probability appropriate for high-density arrays with all available cores
    genoprob <- calc_genoprob(cross2, gmap_pseudo, error_prob=0.002, cores=0)  # cores=0 uses all available

    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)

    # Save genotype probabilities
    saveRDS(genoprob, file = "${prefix}_genoprob.rds")

    validation_log <- c(validation_log, paste("✓ Genotype probabilities calculated in", duration, "minutes"))
    validation_log <- c(validation_log, paste("✓ Genoprob dimensions:", paste(dim(genoprob[[1]]), collapse=" x ")))
    validation_log <- c(validation_log, "")

    # Calculate kinship matrix using LOCO (Leave One Chromosome Out) method
    cat("Calculating LOCO kinship matrices...\\n")
    cat("This step accounts for population structure and relatedness...\\n")

    validation_log <- c(validation_log, "=== LOCO Kinship Matrix Calculation ===")
    start_time <- Sys.time()

    # LOCO kinship for proper mixed model analysis with all available cores
    kinship_loco <- calc_kinship(genoprob, type="loco", cores=0)  # cores=0 uses all available

    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)

    # Save kinship matrices
    saveRDS(kinship_loco, file = "${prefix}_kinship_loco.rds")

    validation_log <- c(validation_log, paste("✓ LOCO kinship matrices calculated in", duration, "minutes"))
    validation_log <- c(validation_log, paste("✓ Number of kinship matrices (one per chromosome):", length(kinship_loco)))

    # Report kinship matrix properties
    if (length(kinship_loco) > 0) {
        first_kinship <- kinship_loco[[1]]
        validation_log <- c(validation_log, paste("✓ Kinship matrix dimensions:", paste(dim(first_kinship), collapse=" x ")))
        validation_log <- c(validation_log, paste("✓ Kinship diagonal range:", round(min(diag(first_kinship)), 3), "to", round(max(diag(first_kinship)), 3)))
    }

    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Genome Scan Preparation Complete ===")
    validation_log <- c(validation_log, "✓ Ready for genome scanning and permutation testing")

    # Write validation report
    writeLines(validation_log, "genoprob_validation_report.txt")

    cat("Genome scan preparation completed successfully\\n")
    """
}