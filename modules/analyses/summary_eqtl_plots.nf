process SUMMARY_EQTL_PLOTS {
    tag "Cross-model eQTL summary figures"
    publishDir "${params.summary_outdir}/Summary", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    cpus 2
    memory '16 GB'
    time '30m'

    input:
    val study_prefix
    val results_base           // e.g. Results_final/eQTL
    val models_str             // comma-separated model names, in plot order

    output:
    path "${study_prefix}_eqtl_cis_trans_barchart.png", emit: barchart
    path "${study_prefix}_eqtl_counts_by_model.csv",    emit: counts
    path "${study_prefix}_eqtl_summary_log.txt",        emit: log

    script:
    """
    #!/usr/bin/env Rscript

    # Cross-model eQTL yield summary: stacked bar chart of distinct cis/trans eQTL
    # counts for each analysis model, so diet- and model-dependent shifts in eQTL
    # yield are directly comparable. Counts use the 95% genome-wide threshold; the
    # classification table nests significance levels, so significance_level == "95%"
    # is the distinct, complete 95%-significant set per model.

    suppressMessages({ library(data.table); library(ggplot2) })

    out_png <- "${study_prefix}_eqtl_cis_trans_barchart.png"
    out_csv <- "${study_prefix}_eqtl_counts_by_model.csv"
    out_log <- "${study_prefix}_eqtl_summary_log.txt"

    log_con <- file(out_log, "w")
    logm <- function(...) { msg <- paste0(...); cat(msg, "\\n"); cat(msg, "\\n", file = log_con); flush(log_con) }

    logm("=== CROSS-MODEL eQTL SUMMARY ===")

    # Plot order is the order models are passed in (single-diet -> additive -> interactive)
    model_levels <- strsplit("${models_str}", ",")[[1]]
    base <- "${results_base}"
    logm("Models: ", paste(model_levels, collapse = ", "))
    logm("Results base: ", base)

    rows <- rbindlist(lapply(model_levels, function(model) {
      f <- file.path(base, model, "00_analyses", "eqtl_cis_trans_classification.csv")
      if (!file.exists(f)) stop("Missing classification CSV for model '", model, "': ", f)
      d <- fread(f)
      d <- d[significance_level == "95%" & eqtl_type %in% c("cis", "trans")]
      logm("  ", model, ": ", nrow(d), " distinct 95% cis/trans eQTLs")
      data.table(model = model,
                 eqtl_type = d\$eqtl_type)[, .(n = .N), by = .(model, eqtl_type)]
    }))

    # Apply fixed order; trans drawn on top of cis
    rows[, model := factor(model, levels = model_levels)]
    rows[, eqtl_type := factor(eqtl_type, levels = c("trans", "cis"))]
    setorder(rows, model, eqtl_type)

    counts <- dcast(rows, model ~ eqtl_type, value.var = "n", fill = 0)
    fwrite(counts, out_csv)
    logm("Wrote count table: ", out_csv)

    totals <- rows[, .(total = sum(n)), by = model]

    p <- ggplot(rows, aes(x = model, y = n, fill = eqtl_type)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5),
                color = "white", size = 3.2) +
      geom_text(data = totals, inherit.aes = FALSE,
                aes(x = model, y = total, label = total), vjust = -0.4, size = 3.5) +
      scale_fill_manual(values = c(cis = "#1f5fbf", trans = "black"),
                        breaks = c("cis", "trans"), name = "eQTL type") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(title = paste0("${study_prefix} — eQTL yield by model (95% GW)"),
           subtitle = "Distinct cis/trans eQTLs per analysis model",
           x = NULL, y = "eQTL count") +
      theme_bw(base_size = 12) +
      theme(panel.grid.major.x = element_blank(),
            legend.position = "right")

    ggsave(out_png, p, width = 8, height = 6, dpi = 300)
    logm("Saved: ", out_png)
    logm("=== DONE ===")
    close(log_con)
    """
}
