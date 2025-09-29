process PREPARE_GENOME_SCAN_SETUP {
    tag "Setting up genome scan preparation for ${prefix}"
    publishDir "${params.outdir}/05_genome_scan_preparation", mode: 'copy'

    cpus 128
    memory '128 GB'
    time '4h'

    input:
    path(cross2_file)
    val(prefix)

    output:
    path("${prefix}_genetic_map.rds"), emit: genetic_map
    path("${prefix}_genoprob.rds"), emit: genoprob
    path("${prefix}_alleleprob.rds"), emit: alleleprob
    path("${prefix}_kinship_loco.rds"), emit: kinship_loco
    path("chromosome_list.txt"), emit: chromosome_list
    path("prep_setup_report.txt"), emit: setup_report

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })

    cat("=== GENOME SCAN PREPARATION SETUP ===\\n")
    cat("Study Prefix: ${prefix}\\n")
    cat("Setup Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # Load cross2 object
    cat("Loading cross2 object...\\n")
    cross2 <- readRDS("${cross2_file}")

    # Get genetic map and insert pseudomarkers
    cat("Preparing genetic map with pseudomarkers...\\n")
    gmap <- cross2\$gmap
    gmap_pseudo <- insert_pseudomarkers(gmap, step=0.5)

    # Save genetic map
    saveRDS(gmap_pseudo, file = "${prefix}_genetic_map.rds")

    # Calculate genotype probabilities for ALL chromosomes
    # This step cannot be easily parallelized due to cross2 object structure requirements
    cat("Calculating genotype probabilities for all chromosomes...\\n")
    cat("This may take some time but runs on multiple cores...\\n")

    start_time <- Sys.time()
    genoprob <- calc_genoprob(cross2, gmap_pseudo, error_prob=0.002, cores=0)
    end_time <- Sys.time()
    duration <- round(as.numeric(end_time - start_time, units="mins"), 2)

    # Save genotype probabilities
    saveRDS(genoprob, file = "${prefix}_genoprob.rds")

    # Convert genotype probabilities to allele probabilities for additive models
    cat("Converting genotype probabilities to allele probabilities...\n")
    start_time_allele <- Sys.time()
    alleleprob <- genoprob_to_alleleprob(genoprob)
    end_time_allele <- Sys.time()
    duration_allele <- round(as.numeric(end_time_allele - start_time_allele, units="mins"), 2)

    # Save allele probabilities
    saveRDS(alleleprob, file = "${prefix}_alleleprob.rds")

    # Get chromosome list
    chromosomes <- names(cross2\$geno)
    writeLines(chromosomes, "chromosome_list.txt")

    # Calculate LOCO kinship matrices using allele probabilities (all chromosomes in one step with multi-core)
    cat("Calculating LOCO kinship matrices for all chromosomes using allele probabilities...\\n")
    start_time_kinship <- Sys.time()
    kinship_loco <- calc_kinship(alleleprob, "loco", cores=0)
    end_time_kinship <- Sys.time()
    duration_kinship <- round(as.numeric(end_time_kinship - start_time_kinship, units="mins"), 2)

    # Save LOCO kinship
    saveRDS(kinship_loco, file = "${prefix}_kinship_loco.rds")

    # Create setup report
    setup_report <- c(
        "=== Genome Scan Preparation Setup ===",
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Study Prefix: ${prefix}"),
        "",
        "=== Genotype Probability Calculation ===",
        paste("Number of individuals:", nrow(cross2\$pheno)),
        paste("Number of chromosomes:", length(chromosomes)),
        paste("Chromosomes:", paste(chromosomes, collapse=", ")),
        paste("Total markers:", sum(sapply(gmap, length))),
        paste("Markers with pseudomarkers:", sum(sapply(gmap_pseudo, length))),
        paste("Genoprob calculation time:", duration, "minutes"),
        "",
        "=== Allele Probability Conversion ===",
        paste("Alleleprob conversion time:", duration_allele, "minutes"),
        "Converted genotype probabilities to allele probabilities for additive models",
        "",
        "=== LOCO Kinship Calculation ===",
        paste("Number of chromosomes:", length(chromosomes)),
        paste("Chromosomes:", paste(chromosomes, collapse=", ")),
        paste("Kinship calculation time:", duration_kinship, "minutes"),
        paste("Total processing time:", round(duration + duration_allele + duration_kinship, 2), "minutes"),
        "",
        "=== Genome Scan Preparation Complete ===",
        "All genotype probabilities, allele probabilities, and LOCO kinship matrices ready for scanning"
    )

    writeLines(setup_report, "prep_setup_report.txt")

    cat("\\n=== GENOME SCAN PREPARATION COMPLETED ===\\n")
    cat("Genoprob calculated in", duration, "minutes\\n")
    cat("Alleleprob converted in", duration_allele, "minutes\\n")
    cat("LOCO kinship calculated in", duration_kinship, "minutes\\n")
    cat("Total time:", round(duration + duration_allele + duration_kinship, 2), "minutes\\n")
    cat("Ready for genome scanning\\n")
    """
}

process PREPARE_GENOME_SCAN_CHR_KINSHIP {
    tag "Calculating LOCO kinship for chromosome ${chromosome} for ${prefix}"
    publishDir "${params.outdir}/05_genome_scan_preparation/kinship", mode: 'copy'

    cpus 4
    memory '16 GB'
    time '1h'

    input:
    path(cross2_file)
    path(genoprob_file)
    val(prefix)
    val(chromosome)

    output:
    path("${prefix}_kinship_chr${chromosome}.rds"), emit: kinship_chr
    path("chr${chromosome}_kinship_report.txt"), emit: chr_report

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(qtl2)
    })

    cat("=== LOCO KINSHIP FOR CHROMOSOME ${chromosome} ===\\n")
    cat("Study Prefix: ${prefix}\\n")
    cat("Chromosome: ${chromosome}\\n")
    cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # Load data
    cross2 <- readRDS("${cross2_file}")
    genoprob_full <- readRDS("${genoprob_file}")

    cat("Calculating LOCO kinship excluding chromosome ${chromosome}...\\n")
    start_time <- Sys.time()

    # For LOCO kinship, use all chromosomes EXCEPT this one
    other_chrs <- names(genoprob_full)[names(genoprob_full) != "${chromosome}"]
    genoprob_loco <- genoprob_full[other_chrs]

    # Calculate kinship from LOCO genoprob
    kinship_chr <- calc_kinship(genoprob_loco, cores=0)

    end_time <- Sys.time()
    duration_kinship <- round(as.numeric(end_time - start_time, units="mins"), 2)

    # Save results
    saveRDS(kinship_chr, file = "${prefix}_kinship_chr${chromosome}.rds")

    # Create chromosome report
    chr_report <- c(
        paste("=== Chromosome ${chromosome} LOCO Kinship Report ==="),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Study Prefix: ${prefix}"),
        paste("Chromosome: ${chromosome}"),
        "",
        "=== Processing Results ===",
        paste("LOCO kinship calculated in:", duration_kinship, "minutes"),
        paste("Chromosomes used for kinship:", paste(other_chrs, collapse=", ")),
        "",
        "=== Data Dimensions ===",
        paste("Kinship matrix dimensions:", paste(dim(kinship_chr), collapse=" x ")),
        paste("Kinship diagonal range:", round(min(diag(kinship_chr)), 3), "to", round(max(diag(kinship_chr)), 3)),
        "",
        "=== Chromosome ${chromosome} LOCO Kinship Complete ==="
    )

    writeLines(chr_report, "chr${chromosome}_kinship_report.txt")

    cat("\\nChromosome ${chromosome} LOCO kinship completed successfully\\n")
    cat("Kinship calculation time:", duration_kinship, "minutes\\n")
    """
}

process COMBINE_GENOME_SCAN_PREP {
    tag "Combining genome scan preparation for ${prefix}"
    publishDir "${params.outdir}/05_genome_scan_preparation", mode: 'copy'

    cpus 2
    memory '8 GB'
    time '30m'

    input:
    path(genetic_map_file)
    path(genoprob_file)
    path(kinship_files)   // All chromosome kinship files
    path(chr_reports)     // All chromosome reports
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

    cat("=== COMBINING GENOME SCAN PREPARATION ===\\n")
    cat("Study Prefix: ${prefix}\\n")
    cat("Combine Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # Copy pre-calculated files
    file.copy("${genetic_map_file}", "${prefix}_genetic_map.rds")
    file.copy("${genoprob_file}", "${prefix}_genoprob.rds")

    # Load genoprob to get chromosome list
    genoprob <- readRDS("${genoprob_file}")
    chromosomes <- names(genoprob)

    cat("Combining LOCO kinship matrices from", length(chromosomes), "chromosomes...\\n")

    # Initialize combined kinship structure
    combined_kinship_loco <- list()

    # Load and combine kinship files
    for (chr in chromosomes) {
        kinship_file <- paste0("${prefix}_kinship_chr", chr, ".rds")
        if (file.exists(kinship_file)) {
            kinship_chr <- readRDS(kinship_file)
            combined_kinship_loco[[chr]] <- kinship_chr
            cat("  Loaded LOCO kinship for chromosome", chr, "\\n")
        } else {
            cat("  WARNING: Missing kinship file for chromosome", chr, "\\n")
        }
    }

    # Set proper class attributes for qtl2
    names(combined_kinship_loco) <- chromosomes

    # Save combined kinship
    saveRDS(combined_kinship_loco, file = "${prefix}_kinship_loco.rds")

    # Create validation report
    validation_log <- c()
    validation_log <- c(validation_log, "=== Genome Scan Preparation Report ===")
    validation_log <- c(validation_log, paste("Timestamp:", Sys.time()))
    validation_log <- c(validation_log, paste("Study Prefix:", "${prefix}"))
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Processing Summary ===")
    validation_log <- c(validation_log, paste("Total chromosomes processed:", length(chromosomes)))
    validation_log <- c(validation_log, paste("Chromosomes:", paste(chromosomes, collapse=", ")))
    validation_log <- c(validation_log, paste("Genoprob components:", length(genoprob)))
    validation_log <- c(validation_log, paste("LOCO kinship matrices combined:", length(combined_kinship_loco)))
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Data Validation ===")

    # Report dimensions for each chromosome
    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Per-Chromosome Results ===")
    for (chr in chromosomes) {
        if (chr %in% names(genoprob)) {
            genoprob_dim <- dim(genoprob[[chr]])
            validation_log <- c(validation_log, paste("Chr", chr, "genoprob:", paste(genoprob_dim, collapse=" x ")))
        }
        if (chr %in% names(combined_kinship_loco)) {
            kinship_dim <- dim(combined_kinship_loco[[chr]])
            validation_log <- c(validation_log, paste("Chr", chr, "LOCO kinship:", paste(kinship_dim, collapse=" x ")))
        }
    }

    validation_log <- c(validation_log, "")
    validation_log <- c(validation_log, "=== Genome Scan Preparation Complete ===")
    validation_log <- c(validation_log, "Ready for genome scanning and permutation testing")

    # Write validation report
    writeLines(validation_log, "genoprob_validation_report.txt")

    cat("\\n=== COMBINATION COMPLETED ===\\n")
    cat("All results successfully combined\\n")
    cat("Ready for genome scanning\\n")
    """
}