#!/usr/bin/env nextflow

// Optional Regional QTL Filtering Module
// Filters QTL peaks to specific genomic regions before permutation testing
// Dramatically reduces computational burden for targeted analyses

process FILTER_PEAKS_BY_REGION {
    tag "Filtering QTL peaks by genomic region"
    publishDir "${params.outdir}/06b_regional_filtering", mode: 'copy'

    cpus 4
    memory '16 GB'
    time '30m'

    input:
    path(all_peaks_file)
    path(genetic_map_file)
    path(gtf_file)
    val(study_prefix)
    val(region_spec)  // Format: "chr12:81-91" or "chr12" or "chr12:81-91;chr2:50-60"

    output:
    path("${study_prefix}_regional_filtered_phenotypes.txt"), emit: filtered_phenotypes
    path("${study_prefix}_regional_filtering_report.txt"), emit: filter_report
    path("${study_prefix}_filtered_peaks.csv"), emit: filtered_peaks

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    log_file <- "${study_prefix}_regional_filtering_report.txt"
    log_conn <- file(log_file, "w")

    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- paste0("[", timestamp, "] ", msg)
        cat(full_msg, "\\n", file = stderr())
        cat(full_msg, "\\n", file = log_conn)
        flush(log_conn)
    }

    log_message("=== REGIONAL QTL PEAK FILTERING ===")
    log_message("Region specification: ${region_spec}")

    # Parse region specification
    parse_regions <- function(spec) {
        # Split multiple regions by semicolon
        regions <- strsplit(spec, ";")[[1]]
        region_list <- lapply(regions, function(r) {
            # Parse chr:start-end or just chr
            if (grepl(":", r)) {
                parts <- strsplit(r, ":")[[1]]
                chr <- gsub("^chr", "", parts[1])  # Remove chr prefix if present
                range_parts <- strsplit(parts[2], "-")[[1]]
                list(
                    chr = chr,
                    start_mbp = as.numeric(range_parts[1]),
                    end_mbp = as.numeric(range_parts[2])
                )
            } else {
                # Whole chromosome
                chr <- gsub("^chr", "", r)  # Remove chr prefix if present
                list(chr = chr, start_mbp = 0, end_mbp = Inf)
            }
        })
        return(region_list)
    }

    target_regions <- parse_regions("${region_spec}")
    log_message(sprintf("Parsed %d target region(s):", length(target_regions)))
    for (i in seq_along(target_regions)) {
        r <- target_regions[[i]]
        if (is.infinite(r\$end_mbp)) {
            log_message(sprintf("  %d. chr%s (entire chromosome)", i, r\$chr))
        } else {
            log_message(sprintf("  %d. chr%s: %.1f-%.1f Mbp", i, r\$chr, r\$start_mbp, r\$end_mbp))
        }
    }

    # Load peaks (has lodcolumn, chr, pos, lod columns from qtl2::find_peaks)
    log_message("Loading QTL peaks...")
    peaks <- fread("${all_peaks_file}")
    total_peaks <- nrow(peaks)
    total_phenos <- length(unique(peaks\$lodcolumn))
    log_message(sprintf("Total peaks: %d", total_peaks))
    log_message(sprintf("Unique phenotypes with peaks: %d", total_phenos))

    # Load GTF and extract gene locations
    log_message("Loading genome annotation (GTF)...")
    gtf <- fread("${gtf_file}", skip = 5, header = FALSE, sep = "\\t")
    colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

    # Filter to gene features only
    gtf_genes <- gtf[feature == "gene"]
    log_message(sprintf("Total genes in GTF: %d", nrow(gtf_genes)))

    # Extract gene_id from attributes column
    log_message("Extracting gene IDs from GTF attributes...")
    gtf_genes\$gene_id <- sapply(strsplit(as.character(gtf_genes\$attributes), ";"), function(attrs) {
        gene_id_field <- grep("gene_id", attrs, value = TRUE)
        if (length(gene_id_field) > 0) {
            gsub('.*gene_id "([^"]+)".*', "\\\\1", gene_id_field[1])
        } else {
            NA
        }
    })

    # Remove chr prefix if present (standardize to match peaks)
    gtf_genes\$chr <- gsub("^chr", "", gtf_genes\$chr)

    # Create gene location table (chr, start_mbp, end_mbp, gene_id)
    gene_locations <- data.table(
        chr = gtf_genes\$chr,
        start_mbp = gtf_genes\$start / 1e6,
        end_mbp = gtf_genes\$end / 1e6,
        gene_id = gtf_genes\$gene_id
    )
    gene_locations <- gene_locations[!is.na(gene_id)]
    log_message(sprintf("Extracted %d gene locations", nrow(gene_locations)))

    # Cross-reference peaks with gene locations
    log_message("Cross-referencing peak phenotypes with GTF gene locations...")

    # Match phenotype IDs (lodcolumn) to gene IDs
    peaks_with_loc <- merge(
        peaks,
        gene_locations,
        by.x = "lodcolumn",
        by.y = "gene_id",
        all.x = TRUE,
        suffixes = c("_peak", "_gene")
    )

    # Report matching statistics
    matched_peaks <- sum(!is.na(peaks_with_loc\$chr_gene))
    matched_phenos <- length(unique(peaks_with_loc\$lodcolumn[!is.na(peaks_with_loc\$chr_gene)]))

    log_message("")
    log_message("=== PHENOTYPE MATCHING STATISTICS ===")
    log_message(sprintf("Matched peaks: %d/%d (%.1f%%)",
                       matched_peaks, total_peaks,
                       100 * matched_peaks / total_peaks))
    log_message(sprintf("Matched phenotypes: %d/%d (%.1f%%)",
                       matched_phenos, total_phenos,
                       100 * matched_phenos / total_phenos))

    if (matched_phenos < total_phenos) {
        unmatched <- unique(peaks\$lodcolumn[!peaks\$lodcolumn %in% gene_locations\$gene_id])
        log_message(sprintf("WARNING: %d phenotypes not found in GTF annotation", length(unmatched)))
        log_message(sprintf("First 10 unmatched phenotypes: %s",
                           paste(head(unmatched, 10), collapse=", ")))
    }

    # Filter to target regions
    log_message("")
    log_message("=== FILTERING TO TARGET REGIONS ===")

    # Vectorized filtering instead of apply()
    peaks_with_loc\$in_region <- FALSE

    for (region in target_regions) {
        # Check if gene overlaps with target region
        # Gene overlaps if: gene_end >= region_start AND gene_start <= region_end
        in_this_region <- !is.na(peaks_with_loc\$chr_gene) &
                         !is.na(peaks_with_loc\$start_mbp) &
                         !is.na(peaks_with_loc\$end_mbp) &
                         peaks_with_loc\$chr_gene == region\$chr &
                         peaks_with_loc\$end_mbp >= region\$start_mbp &
                         peaks_with_loc\$start_mbp <= region\$end_mbp

        peaks_with_loc\$in_region <- peaks_with_loc\$in_region | in_this_region
    }

    filtered_peaks <- peaks_with_loc[peaks_with_loc\$in_region, ]

    log_message(sprintf("Peaks in target regions: %d/%d (%.1f%%)",
                       nrow(filtered_peaks), total_peaks,
                       100 * nrow(filtered_peaks) / total_peaks))

    # Get unique phenotypes
    filtered_phenotypes <- unique(filtered_peaks\$lodcolumn)
    log_message(sprintf("Phenotypes with peaks in target regions: %d/%d (%.1f%%)",
                       length(filtered_phenotypes), total_phenos,
                       100 * length(filtered_phenotypes) / total_phenos))

    # Summary by region
    log_message("")
    log_message("=== SUMMARY BY REGION ===")
    for (i in seq_along(target_regions)) {
        r <- target_regions[[i]]

        region_peaks <- filtered_peaks[
            filtered_peaks\$chr_gene == r\$chr &
            filtered_peaks\$end_mbp >= r\$start_mbp &
            filtered_peaks\$start_mbp <= r\$end_mbp
        ]
        region_phenos <- length(unique(region_peaks\$lodcolumn))

        if (is.infinite(r\$end_mbp)) {
            log_message(sprintf("chr%s: %d peaks, %d phenotypes",
                               r\$chr, nrow(region_peaks), region_phenos))
        } else {
            log_message(sprintf("chr%s:%.1f-%.1f Mbp: %d peaks, %d phenotypes",
                               r\$chr, r\$start_mbp, r\$end_mbp,
                               nrow(region_peaks), region_phenos))
        }
    }

    # Prepare output - keep only original peak columns plus gene location info
    output_peaks <- filtered_peaks[, .(lodcolumn, chr = chr_peak, pos, lod,
                                       gene_chr = chr_gene,
                                       gene_start_mbp = start_mbp, gene_end_mbp = end_mbp)]

    # Save outputs
    writeLines(filtered_phenotypes, "${study_prefix}_regional_filtered_phenotypes.txt")
    fwrite(output_peaks, "${study_prefix}_filtered_peaks.csv")

    log_message("")
    log_message("=== FILTERING COMPLETE ===")
    log_message(sprintf("Workload reduction: %.1f%%",
                       100 * (1 - length(filtered_phenotypes) / total_phenos)))
    log_message(sprintf("Output files:"))
    log_message(sprintf("  - Filtered phenotype list: %s",
                       "${study_prefix}_regional_filtered_phenotypes.txt"))
    log_message(sprintf("  - Filtered peaks: %s",
                       "${study_prefix}_filtered_peaks.csv"))

    close(log_conn)
    """
}

workflow REGIONAL_FILTER {
    take:
    all_peaks
    genetic_map
    gtf
    study_prefix
    region_spec

    main:
    FILTER_PEAKS_BY_REGION(
        all_peaks,
        genetic_map,
        gtf,
        study_prefix,
        region_spec
    )

    emit:
    filtered_phenotypes = FILTER_PEAKS_BY_REGION.out.filtered_phenotypes
    filter_report = FILTER_PEAKS_BY_REGION.out.filter_report
    filtered_peaks = FILTER_PEAKS_BY_REGION.out.filtered_peaks
}
