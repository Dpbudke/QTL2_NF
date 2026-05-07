process CLASSIFY_CIS_TRANS_EQTLS {
    tag "Classifying cis vs trans eQTLs"
    publishDir "${params.outdir}/00_analyses", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    input:
    path qtl_file
    path gtf_file
    path cross2_file

    output:
    path "eqtl_cis_trans_classification.csv", emit: classification
    path "eqtl_classification_summary.txt", emit: summary
    path "qtl2_position_map.rds",            emit: position_map
    path "eqtl_classification_report.html",  emit: report, optional: true

    script:
    def cis_window_mb = params.cis_window_mb ?: 4.0
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(qtl2)

    # Read significant QTLs
    cat("Reading QTL data...\\n")
    qtls <- fread("${qtl_file}")

    # Load cross2 to extract gmap/pmap, then save as reusable position map output
    cat("Loading cross2 object for coordinate conversion...\\n")
    cross2 <- readRDS("${cross2_file}")
    gmap <- cross2\$gmap
    pmap <- cross2\$pmap
    saveRDS(list(gmap = gmap, pmap = pmap), "qtl2_position_map.rds")
    cat("  Position map saved\\n")

    # Convert QTL positions from cM to Mb — vectorized per chromosome
    cat("Converting QTL positions from cM to Mb...\\n")
    qtls[, pos_mb := NA_real_]
    for (chr_name in intersect(unique(qtls\$chr), intersect(names(gmap), names(pmap)))) {
        idx <- qtls\$chr == chr_name
        qtls[idx, pos_mb := approx(gmap[[chr_name]], pmap[[chr_name]],
                                   xout = pos, rule = 2)\$y]
    }
    cat("  Converted", sum(!is.na(qtls\$pos_mb)), "of", nrow(qtls), "QTL positions\\n")

    # Read and process GTF file
    cat("Reading GTF annotation...\\n")
    gtf <- fread("gunzip -c ${gtf_file}", header = FALSE, sep = "\\t")

    # Filter for gene entries only
    genes <- gtf[V3 == "gene", .(chr = V1, start = V4, end = V5, strand = V7, attributes = V9)]

    # Extract gene_id and gene_name from attributes
    genes[, gene_id := gsub('.*gene_id "([^"]+)".*', '\\\\1', attributes)]
    genes[, gene_name := gsub('.*gene_name "([^"]+)".*', '\\\\1', attributes)]
    genes[, gene_name := ifelse(grepl('gene_name', attributes), gene_name, NA)]
    genes[, attributes := NULL]

    # Calculate TSS based on strand
    # Forward strand (+): TSS = start
    # Reverse strand (-): TSS = end
    genes[, tss := ifelse(strand == "+", start, end)]
    genes[, tss_mb := tss / 1e6]  # Convert to Mb for comparison with QTL positions

    cat("Found", nrow(genes), "genes in GTF annotation\\n")
    cat("Processing", nrow(qtls), "significant QTLs\\n")

    # Join QTLs with gene annotations using merge
    qtl_classification <- merge(qtls, genes, by.x = "lodcolumn", by.y = "gene_id", all.x = TRUE)

    # Calculate distance between QTL peak (in Mb) and TSS (in Mb)
    qtl_classification[, same_chr := (chr.x == chr.y)]
    qtl_classification[, distance_to_tss_mb := abs(pos_mb - tss_mb)]

    # Classify as cis or trans
    qtl_classification[, eqtl_type := "trans"]
    qtl_classification[is.na(tss_mb), eqtl_type := "unknown_gene"]
    qtl_classification[!is.na(tss_mb) & same_chr & distance_to_tss_mb <= ${cis_window_mb}, eqtl_type := "cis"]

    # Rename columns for clarity
    setnames(qtl_classification, "chr.x", "qtl_chr")
    setnames(qtl_classification, "chr.y", "gene_chr")
    setnames(qtl_classification, "pos", "qtl_pos_cM")
    setnames(qtl_classification, "pos_mb", "qtl_pos_mb")
    setnames(qtl_classification, "start", "gene_start")
    setnames(qtl_classification, "end", "gene_end")
    setnames(qtl_classification, "strand", "gene_strand")
    setnames(qtl_classification, "tss", "gene_tss")
    setnames(qtl_classification, "tss_mb", "gene_tss_mb")

    # Reorder columns
    col_order <- c("lodindex", "lodcolumn", "gene_name",
                   "qtl_chr", "qtl_pos_cM", "qtl_pos_mb",
                   "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_tss", "gene_tss_mb",
                   "distance_to_tss_mb", "eqtl_type",
                   "lod", "ci_lo", "ci_hi",
                   "signif_63", "signif_90", "signif_95", "signif_99",
                   "study", "significance_level")

    # Only keep columns that exist
    col_order <- col_order[col_order %in% names(qtl_classification)]
    other_cols <- setdiff(names(qtl_classification), col_order)
    setcolorder(qtl_classification, c(col_order, other_cols))

    # Write classification results
    cat("Writing classification results...\\n")
    fwrite(qtl_classification, "eqtl_cis_trans_classification.csv")

    # Generate summary statistics
    summary_stats <- qtl_classification[, .(
        count = .N,
        mean_lod = mean(lod, na.rm = TRUE),
        median_distance_mb = median(distance_to_tss_mb, na.rm = TRUE)
    ), by = eqtl_type]

    summary_stats[, percentage := count / sum(count) * 100]

    # Write summary to text file
    sink("eqtl_classification_summary.txt")
    cat("===========================================\\n")
    cat("eQTL CIS vs TRANS CLASSIFICATION SUMMARY\\n")
    cat("===========================================\\n")
    cat("\\nAnalysis Parameters:\\n")
    cat("  Cis window: +/-", ${cis_window_mb}, "Mb from TSS (", ${cis_window_mb} * 2, "Mb total)\\n")
    cat("  Total QTLs analyzed:", nrow(qtl_classification), "\\n")
    cat("\\nClassification Results:\\n")
    print(summary_stats)
    cat("\\nBy Significance Level:\\n")
    print(table(qtl_classification\$eqtl_type, qtl_classification\$significance_level))
    sink()

    # Print summary to console
    cat("\\n===========================================\\n")
    cat("CLASSIFICATION SUMMARY\\n")
    cat("===========================================\\n")
    print(summary_stats)
    cat("\\nResults written to:\\n")
    cat("  - eqtl_cis_trans_classification.csv\\n")
    cat("  - eqtl_classification_summary.txt\\n")
    """
}
