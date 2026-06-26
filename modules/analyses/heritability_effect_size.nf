// eQTL heritability & effect-size analysis (cross-model summary).
// COMPUTE_HERIT_VAREXP: per-gene narrow-sense h2 (est_herit, single overall kinship) and
// per-eQTL Haley-Knott variance explained (dR2). PLOT_HERIT_EFFECTSIZE: 4 figures.

process COMPUTE_HERIT_VAREXP {
    tag "Heritability + variance explained (${study_prefix})"
    publishDir "${params.summary_outdir}/Summary/Heritability", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    cpus 32
    memory '64 GB'
    time '2h'

    input:
    val study_prefix
    path cross2_rds
    path alleleprob_rds
    path genetic_map_rds
    path classification_csv

    output:
    path "${study_prefix}_herit_varexp.csv", emit: results
    path "${study_prefix}_herit_compute_log.txt", emit: log

    script:
    """
    #!/usr/bin/env Rscript
    suppressMessages({ library(qtl2); library(data.table); library(parallel) })
    ncores <- max(1L, ${task.cpus})

    log_con <- file("${study_prefix}_herit_compute_log.txt", "w")
    logm <- function(...) { m <- paste0(...); cat(m, "\\n"); cat(m, "\\n", file = log_con); flush(log_con) }
    logm("=== HERITABILITY + VARIANCE EXPLAINED ===")

    cross2 <- readRDS("${cross2_rds}")
    aprob  <- readRDS("${alleleprob_rds}")
    gmap   <- readRDS("${genetic_map_rds}")          # pseudomarker map (matches alleleprob)
    clsf   <- fread("${classification_csv}")

    # Covariates (mirror modules/06_qtl_analysis.nf)
    covar_df <- as.data.frame(cross2\$covar)
    if ("coat_color" %in% colnames(covar_df)) covar_df\$coat_color <- NULL
    covar_df <- covar_df[, !sapply(covar_df, function(x) length(unique(na.omit(x))) < 2), drop = FALSE]
    addcovar <- model.matrix(as.formula(paste("~", paste(colnames(covar_df), collapse = " + "))),
                             data = covar_df)[, -1, drop = FALSE]
    logm("addcovar columns: ", ncol(addcovar))

    logm("Building single overall kinship matrix...")
    K <- calc_kinship(aprob, "overall", cores = 0)
    logm("  kinship: ", nrow(K), " x ", ncol(K))

    eq <- clsf[significance_level == "95%" & eqtl_type %in% c("cis","trans")]
    eq <- eq[lodcolumn %in% colnames(cross2\$pheno)]
    genes <- unique(eq\$lodcolumn)
    logm("95% cis/trans eQTLs: ", nrow(eq), "  unique genes: ", length(genes))

    # Narrow-sense heritability per gene
    logm("est_herit over ", length(genes), " genes on ", ncores, " cores...")
    h2v <- unlist(mclapply(genes, function(g)
      tryCatch(as.numeric(est_herit(cross2\$pheno[, g, drop=FALSE], K, addcovar)),
               error = function(e) NA_real_), mc.cores = ncores))
    names(h2v) <- genes

    # Haley-Knott variance explained per eQTL peak (dR2 between full and null)
    r2 <- function(X, y) { f <- lm.fit(X, y); 1 - sum(f\$residuals^2)/sum((y - mean(y))^2) }
    varexp_one <- function(i) {
      chr <- as.character(eq\$qtl_chr[i]); g <- eq\$lodcolumn[i]
      mk  <- names(which.min(abs(gmap[[chr]] - eq\$qtl_pos_cM[i])))
      AP  <- aprob[[chr]][, , mk]; y <- cross2\$pheno[, g]
      ids <- rownames(AP); ids <- ids[!is.na(y[ids])]
      Xn  <- cbind(1, addcovar[ids, , drop=FALSE])
      Xf  <- cbind(Xn, AP[ids, -1, drop=FALSE])     # drop 1 founder as reference
      tryCatch(r2(Xf, y[ids]) - r2(Xn, y[ids]), error = function(e) NA_real_)
    }
    logm("Haley-Knott varexp over ", nrow(eq), " eQTLs...")
    eq[, varexp := unlist(mclapply(seq_len(.N), varexp_one, mc.cores = ncores))]
    eq[, h2 := h2v[lodcolumn]]

    res <- eq[, .(lodcolumn, gene_name, eqtl_type, qtl_chr, qtl_pos_cM, qtl_pos_mb, lod, h2, varexp)]
    fwrite(res, "${study_prefix}_herit_varexp.csv")
    logm("Wrote ${study_prefix}_herit_varexp.csv  rows: ", nrow(res))
    logm("  h2 median: ", round(median(res\$h2, na.rm=TRUE), 3),
         "  varexp median: ", round(median(res\$varexp, na.rm=TRUE), 3))
    logm("=== DONE ===")
    close(log_con)
    """
}

process PLOT_HERIT_EFFECTSIZE {
    tag "Heritability/effect-size figures (${study_prefix})"
    publishDir "${params.summary_outdir}/Summary/Heritability", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    cpus 2
    memory '16 GB'
    time '30m'

    input:
    val study_prefix
    path herit_csv
    val results_base          // absolute path to Results_final/eQTL
    val diet_models_str       // comma-separated: AIN76_only,HC_only,Diet_interactive

    output:
    path "${study_prefix}_heritability_distribution.png",   emit: fig_a
    path "${study_prefix}_herit_vs_varexp.png",             emit: fig_b
    path "${study_prefix}_varexp_cis_trans.png",            emit: fig_c
    path "${study_prefix}_varexp_by_diet_category.png",     emit: fig_d
    path "${study_prefix}_diet_categories.csv",             emit: categories
    path "${study_prefix}_herit_summary.txt",               emit: summary
    path "${study_prefix}_herit_figures_log.txt",           emit: log

    script:
    """
    #!/usr/bin/env Rscript
    suppressMessages({ library(data.table); library(ggplot2) })

    log_con <- file("${study_prefix}_herit_figures_log.txt", "w")
    logm <- function(...) { m <- paste0(...); cat(m, "\\n"); cat(m, "\\n", file = log_con); flush(log_con) }

    herit  <- fread("${herit_csv}")
    base   <- "${results_base}"
    models <- strsplit("${diet_models_str}", ",")[[1]]

    # Detection-based, mutually exclusive diet-dependence categories
    geneset <- function(model) {
      d <- fread(file.path(base, model, "00_analyses/eqtl_cis_trans_classification.csv"))
      unique(d[significance_level == "95%" & eqtl_type %in% c("cis","trans"), lodcolumn])
    }
    AIN <- geneset("AIN76_only"); HC <- geneset("HC_only"); INT <- geneset("Diet_interactive")
    categorize <- function(g) {
      inA <- g %in% AIN; inH <- g %in% HC; inI <- g %in% INT
      ifelse(inA & inH, "constitutive",
      ifelse(inA & !inH, "AIN-specific",
      ifelse(inH & !inA, "HC-specific",
      ifelse(inI, "interactive", "uncategorized"))))
    }
    herit[, category := categorize(lodcolumn)]
    logm("category counts (per eQTL):")
    for (cc in herit[, .N, by=category][order(-N)][, paste0(category, "=", N)]) logm("  ", cc)
    cat_levels <- c("constitutive","AIN-specific","HC-specific","interactive")
    fwrite(herit, "${study_prefix}_diet_categories.csv")

    # Figure A: heritability distribution (unique genes)
    genes <- unique(herit[, .(lodcolumn, h2)])
    med_h2 <- median(genes\$h2, na.rm = TRUE)
    pA <- ggplot(genes, aes(x = h2)) +
      geom_histogram(bins = 40, fill = "#4575b4", color = "white") +
      geom_vline(xintercept = med_h2, linetype = "dashed", color = "#d62728", linewidth = 0.7) +
      annotate("text", x = med_h2, y = Inf, label = sprintf("  median = %.2f", med_h2),
               hjust = 0, vjust = 1.5, color = "#d62728") +
      labs(title = paste0("${study_prefix} — narrow-sense heritability of eQTL genes (additive model)"),
           subtitle = sprintf("n = %d unique eQTL genes", nrow(genes)),
           x = expression(paste("Narrow-sense heritability (", h^2, ")")), y = "Number of genes") +
      theme_bw(base_size = 12)
    ggsave("${study_prefix}_heritability_distribution.png", pA, width = 8, height = 6, dpi = 300)

    # Figure B: h2 vs varexp, colored by LOD
    rho <- cor(herit\$h2, herit\$varexp, method = "spearman", use = "complete.obs")
    pB <- ggplot(herit, aes(x = h2, y = varexp, color = lod)) +
      geom_point(size = 0.8, alpha = 0.5) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
      scale_color_viridis_c(name = "LOD", trans = "log10") +
      annotate("text", x = -Inf, y = Inf, label = sprintf("  Spearman rho = %.2f", rho),
               hjust = 0, vjust = 1.5) +
      labs(title = paste0("${study_prefix} — heritability vs. variance explained by peak"),
           subtitle = sprintf("n = %d eQTLs", nrow(herit)),
           x = expression(paste("Narrow-sense heritability (", h^2, ")")),
           y = "Variance explained by peak marker (dR2)") +
      theme_bw(base_size = 12)
    ggsave("${study_prefix}_herit_vs_varexp.png", pB, width = 8, height = 6, dpi = 300)

    # Figure C: varexp cis vs trans (Wilcoxon rank-sum)
    wc <- wilcox.test(varexp ~ eqtl_type, data = herit)
    nC <- herit[, .N, by = eqtl_type]
    pC <- ggplot(herit, aes(x = eqtl_type, y = varexp, fill = eqtl_type)) +
      geom_boxplot(outlier.size = 0.4, width = 0.6) +
      scale_fill_manual(values = c(cis = "#1f5fbf", trans = "black"), guide = "none") +
      annotate("text", x = 1.5, y = max(herit\$varexp, na.rm = TRUE),
               label = sprintf("Wilcoxon p = %s", format.pval(wc\$p.value, digits = 2)), vjust = 1) +
      labs(title = paste0("${study_prefix} — variance explained: cis vs trans"),
           subtitle = sprintf("cis n=%d, trans n=%d", nC[eqtl_type=="cis", N], nC[eqtl_type=="trans", N]),
           x = NULL, y = "Variance explained by peak marker (dR2)") +
      theme_bw(base_size = 12)
    ggsave("${study_prefix}_varexp_cis_trans.png", pC, width = 6, height = 6, dpi = 300)

    # Figure D: varexp by diet-dependence category (Kruskal-Wallis)
    dd <- herit[category %in% cat_levels]
    dd[, category := factor(category, levels = cat_levels)]
    kw <- kruskal.test(varexp ~ category, data = dd)
    nD <- dd[, .N, by = category]
    pD <- ggplot(dd, aes(x = category, y = varexp, fill = category)) +
      geom_boxplot(outlier.size = 0.4, width = 0.6) +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      geom_text(data = nD, aes(x = category, y = 0, label = paste0("n=", N)), vjust = 1.6, size = 3, inherit.aes = FALSE) +
      annotate("text", x = 2.5, y = max(dd\$varexp, na.rm = TRUE),
               label = sprintf("Kruskal-Wallis p = %s", format.pval(kw\$p.value, digits = 2)), vjust = 1) +
      labs(title = paste0("${study_prefix} — variance explained by diet-dependence category"),
           x = NULL, y = "Variance explained by peak marker (dR2)") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave("${study_prefix}_varexp_by_diet_category.png", pD, width = 7, height = 6, dpi = 300)

    logm("Saved 4 figures. Fig C medians cis=", round(herit[eqtl_type=='cis', median(varexp, na.rm=TRUE)],3),
         " trans=", round(herit[eqtl_type=='trans', median(varexp, na.rm=TRUE)],3))

    # ---- Summary stats text file (the statistics backing Figures A–D) ----
    sm <- file("${study_prefix}_herit_summary.txt", "w")
    w   <- function(...) cat(paste0(..., "\\n"), file = sm)
    f3  <- function(x) sprintf("%.3f", as.numeric(x))
    ct  <- suppressWarnings(cor.test(herit\$h2, herit\$varexp, method = "spearman"))

    w("================================================================")
    w(" ${study_prefix} - eQTL Heritability & Effect-Size Summary")
    w(" Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    w("================================================================")
    w("Source: ${herit_csv}")
    w("Total eQTLs (95% GW, cis/trans): ", nrow(herit))
    w("Unique eQTL genes               : ", nrow(genes))
    w("")

    h2q <- quantile(genes\$h2, c(0, .25, .5, .75, 1), na.rm = TRUE)
    w("--- Figure A: narrow-sense heritability (h2) distribution, unique genes ---")
    w("  n genes (non-NA h2): ", sum(!is.na(genes\$h2)))
    w("  mean h2  : ", f3(mean(genes\$h2, na.rm = TRUE)))
    w("  median h2: ", f3(h2q[3]))
    w("  IQR      : ", f3(h2q[2]), " - ", f3(h2q[4]))
    w("  range    : ", f3(h2q[1]), " - ", f3(h2q[5]))
    w("")

    w("--- Figure B: heritability vs. variance explained by peak (dR2) ---")
    w("  n eQTLs (complete pairs): ", sum(complete.cases(herit\$h2, herit\$varexp)))
    w("  Spearman rho           : ", f3(rho))
    w("  Spearman p             : ", format.pval(ct\$p.value, digits = 3))
    w("")

    cis_v <- herit[eqtl_type == "cis", varexp]; trans_v <- herit[eqtl_type == "trans", varexp]
    w("--- Figure C: variance explained (dR2), cis vs trans ---")
    w("  cis  : n=", nC[eqtl_type == "cis", N],   "  median dR2=", f3(median(cis_v, na.rm = TRUE)),   "  mean=", f3(mean(cis_v, na.rm = TRUE)))
    w("  trans: n=", nC[eqtl_type == "trans", N], "  median dR2=", f3(median(trans_v, na.rm = TRUE)), "  mean=", f3(mean(trans_v, na.rm = TRUE)))
    w("  Wilcoxon rank-sum: W=", as.numeric(wc\$statistic), "  p=", format.pval(wc\$p.value, digits = 3))
    w("")

    w("--- Figure D: variance explained (dR2) by diet-dependence category ---")
    catstat <- dd[, .(n = .N, median_dR2 = median(varexp, na.rm = TRUE), mean_dR2 = mean(varexp, na.rm = TRUE)), by = category]
    catstat <- catstat[order(factor(category, levels = cat_levels))]
    for (i in seq_len(nrow(catstat)))
      w(sprintf("  %-14s n=%-4d median dR2=%.3f  mean=%.3f",
                catstat\$category[i], catstat\$n[i], catstat\$median_dR2[i], catstat\$mean_dR2[i]))
    w("  Kruskal-Wallis: chi-sq=", f3(kw\$statistic), "  df=", as.numeric(kw\$parameter), "  p=", format.pval(kw\$p.value, digits = 3))
    w("")
    w("Diet-dependence category counts across all eQTLs (incl. uncategorized):")
    for (cc in herit[, .N, by = category][order(-N)][, paste0("  ", category, " = ", N)]) w(cc)
    close(sm)
    logm("Wrote ${study_prefix}_herit_summary.txt")

    close(log_con)
    """
}
