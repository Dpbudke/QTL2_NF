process PLOT_EQTL_MAP {
    tag "Plotting genome-wide eQTL map"
    publishDir "${params.outdir}/00_analyses", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    cpus 2
    memory '16 GB'
    time '30m'

    input:
    path classification_csv
    path position_map_rds
    val study_prefix

    output:
    path "${study_prefix}_eqtl_map.png",     emit: map
    path "${study_prefix}_eqtl_map_log.txt", emit: log

    script:
    def analysis_type = params.analysis_type ?: (params.study_type ?: 'eQTL')
    """
    #!/usr/bin/env Rscript

    # Genome-wide eQTL map: x = eQTL peak position, y = gene TSS position (both
    # concatenated across chromosomes). cis = blue diagonal, trans = gray off-diagonal.
    # Filters to eQTLs significant at the 95% genome-wide threshold or stronger.

    suppressMessages({ library(data.table); library(ggplot2) })

    out_png <- "${study_prefix}_eqtl_map.png"
    out_log <- "${study_prefix}_eqtl_map_log.txt"

    log_con <- file(out_log, "w")
    logm <- function(...) { msg <- paste0(...); cat(msg, "\\n"); cat(msg, "\\n", file = log_con); flush(log_con) }

    logm("=== GENOME-WIDE eQTL MAP ===")
    logm("Classification CSV: ", "${classification_csv}")

    d <- fread("${classification_csv}")
    logm("Total rows in classification: ", nrow(d))

    # Keep eQTLs significant at the 95% genome-wide threshold. The classification
    # table nests significance levels (a QTL passing 95% has rows at 63/90/95, plus
    # 99 if it passes that), so significance_level == "95%" is the distinct, complete
    # 95%-significant set (the stronger 99% eQTLs carry a "95%" row too).
    d <- d[significance_level == "95%" & eqtl_type %in% c("cis", "trans") & !is.na(gene_tss_mb)]
    logm("Rows after 95% filter (significance_level == 95%, cis/trans, gene pos known): ", nrow(d))

    # Chromosome offsets from the marker physical map (shared across both axes)
    pm   <- readRDS("${position_map_rds}")
    pmap <- pm\$pmap
    chr_order <- c(as.character(1:19), "X")
    chr_order <- chr_order[chr_order %in% names(pmap)]
    chr_len   <- sapply(chr_order, function(c) max(pmap[[c]], na.rm = TRUE))
    offset    <- setNames(cumsum(c(0, head(chr_len, -1))), chr_order)
    chr_mid   <- offset + chr_len / 2
    boundaries <- c(offset, sum(chr_len))                # gridline positions
    logm("Chromosomes mapped: ", paste(chr_order, collapse = ","))

    # Drop genes/QTLs on chromosomes not in the marker map (e.g. MT)
    d <- d[as.character(qtl_chr) %in% chr_order & as.character(gene_chr) %in% chr_order]
    logm("Rows after restricting to mapped chromosomes: ", nrow(d))

    d[, qtl_gw  := offset[as.character(qtl_chr)]  + qtl_pos_mb]
    d[, gene_gw := offset[as.character(gene_chr)] + gene_tss_mb]

    n_cis   <- sum(d\$eqtl_type == "cis")
    n_trans <- sum(d\$eqtl_type == "trans")
    logm("Plotted cis:   ", n_cis)
    logm("Plotted trans: ", n_trans)

    # Draw trans first so the blue cis diagonal sits on top
    d[, eqtl_type := factor(eqtl_type, levels = c("trans", "cis"))]
    setorder(d, eqtl_type)

    # Alternating chromosome bands (shade odd chromosomes); applied on both axes so
    # their overlap forms a checkerboard that makes per-chromosome blocks easy to read.
    band_df <- data.frame(lo = offset, hi = offset + chr_len,
                          shade = seq_along(chr_order) %% 2 == 1)
    band_df <- band_df[band_df\$shade, ]

    p <- ggplot(d, aes(x = qtl_gw, y = gene_gw, color = eqtl_type)) +
      geom_rect(data = band_df, inherit.aes = FALSE,
                aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
                fill = "grey85", alpha = 0.45) +
      geom_rect(data = band_df, inherit.aes = FALSE,
                aes(ymin = lo, ymax = hi, xmin = -Inf, xmax = Inf),
                fill = "grey85", alpha = 0.45) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey60", linewidth = 0.3) +
      geom_point(size = 0.6, alpha = 0.5) +
      scale_color_manual(values = c(cis = "#1f5fbf", trans = "black"),
                         breaks = c("cis", "trans"),
                         labels = c(paste0("cis (", n_cis, ")"),
                                    paste0("trans (", n_trans, ")")),
                         name = "eQTL type") +
      scale_x_continuous(breaks = chr_mid, labels = chr_order, expand = c(0.01, 0.01)) +
      scale_y_continuous(breaks = chr_mid, labels = chr_order, expand = c(0.01, 0.01)) +
      labs(title = paste0("${study_prefix} — ${analysis_type} eQTL map (95% GW)"),
           x = "eQTL peak position (chromosome)",
           y = "Gene position (chromosome)") +
      coord_fixed() +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            legend.position = "right") +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

    ggsave(out_png, p, width = 10, height = 10, dpi = 300)
    logm("Saved: ", out_png)
    logm("=== DONE ===")
    close(log_con)
    """
}
