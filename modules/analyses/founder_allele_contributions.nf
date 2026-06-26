// Founder allele contributions to eQTL (additive-model, 99% GW catalog).
// COMPUTE_FOUNDER_EFFECTS: per-eQTL qtl2 BLUP founder effects at the peak marker + TIMBR
//   allele-collapsed posterior effects. PLOT_FOUNDER_CONTRIBUTIONS: 7 figures (A, B, C1, C2, D1-D3).

process COMPUTE_FOUNDER_EFFECTS {
    tag "Founder effects: qtl2 BLUP + TIMBR (${study_prefix})"
    publishDir "${params.summary_outdir}/Summary/FounderAlleles", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    cpus 32
    memory '64 GB'
    time '4h'

    input:
    val study_prefix
    path cross2_rds
    path alleleprob_rds
    path genetic_map_rds
    path kinship_loco_rds
    path classification_csv
    path timbr_master_csv
    val timbr_dir            // absolute path to <additive_model>/10_timbr
    val sig_level            // e.g. "99%"

    output:
    path "${study_prefix}_founder_blup.csv",     emit: blup
    path "${study_prefix}_timbr_effects.csv",    emit: timbr
    path "${study_prefix}_founder_effects.csv",  emit: merged
    path "${study_prefix}_founder_compute_log.txt", emit: log

    script:
    """
    #!/usr/bin/env Rscript
    suppressMessages({ library(qtl2); library(data.table); library(parallel) })
    ncores <- max(1L, ${task.cpus})
    FOUNDERS <- c("A","B","C","D","E","F","G","H")

    log_con <- file("${study_prefix}_founder_compute_log.txt", "w")
    logm <- function(...) { m <- paste0(...); cat(m, "\\n"); cat(m, "\\n", file = log_con); flush(log_con) }
    logm("=== FOUNDER ALLELE CONTRIBUTIONS: COMPUTE ===")

    cross2  <- readRDS("${cross2_rds}")
    aprob   <- readRDS("${alleleprob_rds}")
    gmap    <- readRDS("${genetic_map_rds}")
    kloco   <- readRDS("${kinship_loco_rds}")
    clsf    <- fread("${classification_csv}")
    tmaster <- fread("${timbr_master_csv}")
    timbr_dir <- "${timbr_dir}"
    sig <- "${sig_level}"

    # Covariates (mirror modules/06_qtl_analysis.nf); Diet binary: ain76a is the same diet as ain76
    covar_df <- as.data.frame(cross2\$covar)
    if ("coat_color" %in% colnames(covar_df)) covar_df\$coat_color <- NULL
    covar_df <- covar_df[, !sapply(covar_df, function(x) length(unique(na.omit(x))) < 2), drop = FALSE]
    diet_col <- grep("diet", colnames(covar_df), ignore.case = TRUE, value = TRUE)[1]
    dv <- tolower(as.character(covar_df[[diet_col]])); dv[dv %in% c("ain76","ain76a")] <- "ain76"
    covar_df[[diet_col]] <- factor(dv, levels = c("ain76","hc"))
    addcovar <- model.matrix(as.formula(paste("~", paste(colnames(covar_df), collapse = " + "))),
                             data = covar_df)[, -1, drop = FALSE]
    logm("addcovar columns: ", ncol(addcovar))

    eq <- clsf[significance_level == sig & eqtl_type %in% c("cis","trans")]
    eq <- eq[lodcolumn %in% colnames(cross2\$pheno)]
    logm(sig, " cis/trans eQTLs: ", nrow(eq), "  unique genes: ", length(unique(eq\$lodcolumn)))

    ## (1) qtl2 BLUP founder effects at the peak marker, per chromosome (single-marker scan1blup)
    logm("scan1blup at peak markers on ", ncores, " cores (per-chromosome)...")
    blup_mat <- matrix(NA_real_, nrow(eq), 8, dimnames = list(NULL, FOUNDERS))
    chrs <- intersect(unique(as.character(eq\$qtl_chr)), names(gmap))
    for (chr in chrs) {
      ap_c <- aprob[, chr]; K_c <- kloco[[chr]]
      idx  <- which(as.character(eq\$qtl_chr) == chr)
      M <- mclapply(idx, function(i) {
        mk  <- names(which.min(abs(gmap[[chr]] - eq\$qtl_pos_cM[i])))
        ap1 <- ap_c; ap1[[chr]] <- ap_c[[chr]][, , mk, drop = FALSE]
        ph  <- cross2\$pheno[, eq\$lodcolumn[i], drop = FALSE]
        bl  <- tryCatch(scan1blup(ap1, ph, kinship = K_c, addcovar = addcovar), error = function(e) NULL)
        if (is.null(bl)) rep(NA_real_, 8) else as.numeric(bl[1, FOUNDERS])
      }, mc.cores = ncores)
      blup_mat[idx, ] <- do.call(rbind, M)
      logm("  chr ", chr, ": ", length(idx), " loci")
    }

    blup_dt <- data.table(lodcolumn = eq\$lodcolumn, gene_name = eq\$gene_name, eqtl_type = eq\$eqtl_type,
                          chr = eq\$qtl_chr, pos = eq\$qtl_pos_cM, peak_Mb = eq\$qtl_pos_mb, lod = eq\$lod)
    for (f in FOUNDERS) blup_dt[[paste0("blup_", f)]] <- blup_mat[, f]
    # Center-scale per locus: z-score across the 8 founders
    Z <- t(apply(blup_mat, 1, function(v) { s <- sd(v); if (is.na(s) || s == 0) rep(NA_real_, 8) else (v - mean(v)) / s }))
    for (j in seq_along(FOUNDERS)) blup_dt[[paste0("z_", FOUNDERS[j])]] <- Z[, j]
    fwrite(blup_dt, "${study_prefix}_founder_blup.csv")
    logm("Wrote founder_blup.csv rows: ", nrow(blup_dt), "  (non-NA effects: ", sum(!is.na(blup_mat[, 1])), ")")

    ## (2) TIMBR allele-collapsed posterior effects (colMeans of post.hap.effects), per locus
    logm("Reading TIMBR per-locus rds for ", nrow(tmaster), " loci...")
    timbr_one <- function(i) {
      gene <- tmaster\$lodcolumn[i]; chr <- as.character(tmaster\$chr[i])
      dd  <- file.path(timbr_dir, paste0("chr", chr))
      hit <- list.files(dd, pattern = paste0("^", gene, "_chr", chr, "_"), full.names = TRUE)
      if (!length(hit)) return(rep(NA_real_, 8))
      if (length(hit) > 1) {                       # pick dir whose Mb is closest to peak_Mb
        mb  <- suppressWarnings(as.numeric(sub("Mb\$", "", sub(".*_", "", basename(hit)))))
        hit <- hit[which.min(abs(mb - tmaster\$peak_Mb[i]))]
      }
      rds <- list.files(hit[1], pattern = "rds\$", full.names = TRUE)
      if (!length(rds)) return(rep(NA_real_, 8))
      tr <- tryCatch(readRDS(rds[1]), error = function(e) NULL)
      if (is.null(tr) || is.null(tr\$post.hap.effects)) return(rep(NA_real_, 8))
      colMeans(tr\$post.hap.effects)
    }
    TE <- do.call(rbind, mclapply(seq_len(nrow(tmaster)), timbr_one, mc.cores = ncores))
    timbr_dt <- data.table(lodcolumn = tmaster\$lodcolumn, chr = tmaster\$chr, pos = tmaster\$pos,
                          peak_Mb = tmaster\$peak_Mb, lod = tmaster\$lod,
                          n_functional_alleles = tmaster\$n_functional_alleles,
                          top_series_pattern = tmaster\$top_series_pattern)
    for (j in seq_along(FOUNDERS)) timbr_dt[[paste0("timbr_", FOUNDERS[j])]] <- TE[, j]
    fwrite(timbr_dt, "${study_prefix}_timbr_effects.csv")
    logm("Wrote timbr_effects.csv rows: ", nrow(timbr_dt), "  (non-NA effects: ", sum(!is.na(TE[, 1])), ")")

    ## (3) Merge BLUP + TIMBR on lodcolumn + chr (one TIMBR locus per gene/chr)
    tsub <- timbr_dt[, .SD[1], by = .(lodcolumn, chr),
                     .SDcols = c("n_functional_alleles","top_series_pattern", paste0("timbr_", FOUNDERS))]
    merged <- merge(blup_dt, tsub, by = c("lodcolumn","chr"), all.x = TRUE)
    fwrite(merged, "${study_prefix}_founder_effects.csv")
    logm("Merged founder_effects.csv rows: ", nrow(merged), "  with TIMBR: ", sum(!is.na(merged\$timbr_A)))
    logm("=== DONE ===")
    close(log_con)
    """
}

process PLOT_FOUNDER_CONTRIBUTIONS {
    tag "Founder contribution figures (${study_prefix})"
    publishDir "${params.summary_outdir}/Summary/FounderAlleles", mode: 'copy'

    container "file:///${projectDir}/singularity_cache/dpbudke-qtl2-pipeline-latest.img"

    cpus 2
    memory '16 GB'
    time '30m'

    input:
    val study_prefix
    path founder_blup_csv
    path timbr_effects_csv
    path founder_effects_csv

    // All figures (A, B, C1, C2, D1-D3, E) are produced. The TIMBR per-locus
    // post.hap.effects founder columns are now correctly aligned to the qtl2
    // canonical A-H order (P-matrix encoding bug fixed in modules/10_timbr.nf,
    // TIMBR rerun at 95% GW).
    output:
    path "${study_prefix}_founder_effect_boxplots.png",     emit: fig_a
    path "${study_prefix}_founder_effects_circos.png",      emit: fig_b
    path "${study_prefix}_timbr_singleton_frequency.png",   emit: fig_c1
    path "${study_prefix}_timbr_higheffect_frequency.png",  emit: fig_c2
    path "${study_prefix}_blup_vs_timbr_topLOD.png",        emit: fig_d1
    path "${study_prefix}_blup_vs_timbr_wildderived.png",   emit: fig_d2
    path "${study_prefix}_blup_vs_timbr_spanalleles.png",   emit: fig_d3
    path "${study_prefix}_blup_vs_timbr_concordance.png",   emit: fig_e
    path "${study_prefix}_founder_summary.txt",             emit: summary
    path "${study_prefix}_founder_figures_log.txt",         emit: log

    script:
    """
    #!/usr/bin/env Rscript
    suppressMessages({ library(data.table); library(ggplot2); library(circlize) })
    FOUNDERS <- c("A","B","C","D","E","F","G","H")
    SNAME <- c(A="A/J", B="C57BL/6J", C="129S1", D="NOD", E="NZO", F="CAST", G="PWK", H="WSB")
    # CC/DO consensus founder color palette (keyed by founder code)
    PAL <- c(A="#F0F000", B="#808080", C="#F08080", D="#1010F0",
             E="#00A0F0", F="#00A000", G="#F00000", H="#9000E0")
    PAL_S <- setNames(PAL, SNAME[names(PAL)])
    WILD <- c("F","G")  # CAST, PWK

    log_con <- file("${study_prefix}_founder_figures_log.txt", "w")
    logm <- function(...) { m <- paste0(...); cat(m, "\\n"); cat(m, "\\n", file = log_con); flush(log_con) }
    sm_con <- file("${study_prefix}_founder_summary.txt", "w")
    w <- function(...) cat(paste0(..., "\\n"), file = sm_con)

    blup   <- fread("${founder_blup_csv}")
    timbr  <- fread("${timbr_effects_csv}")
    merged <- fread("${founder_effects_csv}")

    parse_pat <- function(p) as.integer(strsplit(p, ",", fixed = TRUE)[[1]])

    # Manual BH-corrected Dunn post-hoc (dunn.test/FSA not in container)
    dunn_bh <- function(values, groups) {
      ok <- !is.na(values); values <- values[ok]; groups <- factor(groups[ok])
      N <- length(values); r <- rank(values)
      Rbar <- tapply(r, groups, mean); n <- tapply(r, groups, length)
      tv <- as.numeric(table(values)); tie <- sum(tv^3 - tv)
      sig2 <- (N * (N + 1) / 12) - tie / (12 * (N - 1))
      prs <- combn(levels(groups), 2)
      z <- apply(prs, 2, function(p) (Rbar[p[1]] - Rbar[p[2]]) / sqrt(sig2 * (1/n[p[1]] + 1/n[p[2]])))
      data.frame(pair = apply(prs, 2, paste, collapse = " vs "),
                 z = z, p_BH = p.adjust(2 * pnorm(-abs(z)), "BH"), row.names = NULL)
    }

    w("================================================================")
    w(" ${study_prefix} - Founder Allele Contributions to eQTL")
    w(" Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    w("================================================================")
    w("qtl2 BLUP eQTL (95% GW): ", nrow(blup), "  TIMBR loci: ", nrow(timbr),
      "  merged with TIMBR: ", sum(!is.na(merged\$timbr_A)))
    w("")

    ## ---- Figure A: per-founder mean-centered BLUP effect boxplots (conventional A/J..WSB order) ----
    bcols <- paste0("blup_", FOUNDERS)
    Bm <- as.matrix(blup[, ..bcols])
    MC <- Bm - rowMeans(Bm, na.rm = TRUE)               # mean-center per locus across the 8 founders
    mcdt <- as.data.table(MC); setnames(mcdt, FOUNDERS); mcdt[, .row := .I]
    long <- melt(mcdt, id.vars = ".row", measure.vars = FOUNDERS, variable.name = "founder", value.name = "mc")
    long <- long[!is.na(mc)]; long[, founder := as.character(founder)]
    long[, strain := factor(SNAME[founder], levels = SNAME[FOUNDERS])]   # fixed conventional founder order
    kw <- kruskal.test(mc ~ founder, data = long)
    dn <- dunn_bh(long\$mc, long\$founder)
    medv <- long[, .(m = median(mc)), by = founder]; medvec <- setNames(medv\$m, medv\$founder)
    pA <- ggplot(long, aes(x = strain, y = mc, fill = strain)) +
      geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
      geom_boxplot(outlier.size = 0.2, width = 0.65) +
      scale_fill_manual(values = PAL_S, guide = "none") +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
               label = sprintf("Kruskal-Wallis p = %s", format.pval(kw\$p.value, digits = 2))) +
      labs(title = "${study_prefix} - mean-centered founder BLUP effects at eQTL peaks (95% GW)",
           subtitle = sprintf("n = %d eQTL; wild-derived CAST/PWK in green/red", nrow(blup)),
           x = NULL, y = "Mean-centered founder BLUP effect") +
      theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave("${study_prefix}_founder_effect_boxplots.png", pA, width = 9, height = 6, dpi = 300)
    w("--- Figure A: mean-centered founder BLUP effects (median, A/J..WSB order) ---")
    for (f in FOUNDERS) w(sprintf("  %-9s median = %+.3f", SNAME[f], medvec[f]))
    w("  Kruskal-Wallis: chi-sq = ", round(kw\$statistic, 1), "  df = ", kw\$parameter,
      "  p = ", format.pval(kw\$p.value, digits = 3))
    w("  Dunn post-hoc (BH) pairs involving CAST(F)/PWK(G):")
    for (i in which(grepl("F|G", dn\$pair))) w(sprintf("    %-9s z = %+6.2f  p_BH = %s", dn\$pair[i], dn\$z[i], format.pval(dn\$p_BH[i], digits = 2)))
    w("")
    logm("Figure A done (KW p=", format.pval(kw\$p.value, digits = 3), ")")

    ## ---- Figure B: circular genome plot, 8 founder tracks of center-scaled effects ----
    chr_order <- c(as.character(1:19), "X")
    cb <- blup[as.character(chr) %in% chr_order & !is.na(z_A)]
    cb[, chrf := factor(as.character(chr), levels = chr_order)]
    xr <- cb[, .(xmax = max(peak_Mb, na.rm = TRUE)), by = chrf][order(chrf)]
    png("${study_prefix}_founder_effects_circos.png", width = 2400, height = 2400, res = 300)
    circos.clear()
    circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1, start.degree = 90, track.margin = c(0, 0.004))
    circos.initialize(factors = xr\$chrf, xlim = cbind(0, xr\$xmax))
    # chromosome labels on an outer thin track
    circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
      panel.fun = function(x, y) circos.text(CELL_META\$xcenter, 0.5, CELL_META\$sector.index, cex = 0.6))
    zlim <- quantile(unlist(cb[, paste0("z_", FOUNDERS), with = FALSE]), c(0.01, 0.99), na.rm = TRUE)
    for (f in FOUNDERS) {
      circos.track(ylim = zlim, track.height = 0.085, bg.border = "grey85",
        panel.fun = function(x, y) {
          ci <- CELL_META\$sector.index
          dd <- cb[chrf == ci]
          circos.lines(c(0, CELL_META\$xlim[2]), c(0, 0), col = "grey70", lwd = 0.4)
          circos.points(dd\$peak_Mb, pmin(pmax(dd[[paste0("z_", f)]], zlim[1]), zlim[2]),
                        pch = 16, cex = 0.12, col = PAL[f])
        })
      circos.text(0, mean(zlim), SNAME[f], sector.index = xr\$chrf[1], facing = "downward",
                  niceFacing = TRUE, adj = c(1.4, 0.5), cex = 0.5, col = PAL[f])
    }
    circos.clear()
    title("${study_prefix} - center-scaled founder effects across the genome (95% eQTL)", cex.main = 0.9)
    dev.off()
    logm("Figure B (circos) done")

    ## ---- Figure C1: TIMBR founder resolution (singleton frequency) ----
    ## Uses top_series_pattern (allele grouping), validated A-H ordered after the
    ## 10_timbr.nf founder-order fix + 95% rerun.
    ti <- timbr[!is.na(top_series_pattern) & top_series_pattern != "" & !is.na(timbr_A)]
    pat <- lapply(ti\$top_series_pattern, parse_pat)
    keep <- vapply(pat, function(v) length(v) == 8, logical(1))
    ti <- ti[keep]; pat <- pat[keep]
    singleton <- t(vapply(pat, function(v) { tb <- table(v); as.numeric(tb[as.character(v)] == 1) }, numeric(8)))
    sing_freq <- colMeans(singleton); names(sing_freq) <- FOUNDERS

    bar_fig <- function(freq, ttl, ylab, file) {
      df <- data.table(founder = FOUNDERS, strain = SNAME[FOUNDERS], freq = freq)
      df[, strain := factor(strain, levels = SNAME[FOUNDERS])]   # conventional A/J..WSB order
      p <- ggplot(df, aes(x = strain, y = freq, fill = strain)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = PAL_S, guide = "none") +
        labs(title = ttl, subtitle = sprintf("n = %d TIMBR loci (95%% GW)", nrow(ti)), x = NULL, y = ylab) +
        theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
      ggsave(file, p, width = 8, height = 6, dpi = 300)
    }
    bar_fig(sing_freq, "${study_prefix} - TIMBR: founder resolved as its own allele",
            "Fraction of loci where founder is a singleton allele", "${study_prefix}_timbr_singleton_frequency.png")
    w("--- Figure C1: TIMBR singleton frequency (own allele, A/J..WSB order) ---")
    for (f in FOUNDERS) w(sprintf("  %-9s %.3f", SNAME[f], sing_freq[f]))
    w("")
    logm("Figure C1 done")

    ## =====================================================================
    ## Figures C2, D1-D3, E use TIMBR per-locus post.hap.effects, whose founder
    ## columns are aligned with the qtl2 canonical A-H order (P-matrix encoding
    ## fixed in modules/10_timbr.nf; TIMBR rerun at 95% GW).
    ## =====================================================================
    {
      teff <- as.matrix(ti[, paste0("timbr_", FOUNDERS), with = FALSE])
      highgrp <- t(vapply(seq_along(pat), function(i) {
        v <- pat[[i]]; ge <- tapply(abs(teff[i, ]), v, mean); mg <- names(which.max(ge))
        as.numeric(as.character(v) == mg)
      }, numeric(8)))
      high_freq <- colMeans(highgrp); names(high_freq) <- FOUNDERS
      bar_fig(high_freq, "${study_prefix} - TIMBR: founder in the highest-effect allele group",
              "Fraction of loci in the max-|effect| allele group", "${study_prefix}_timbr_higheffect_frequency.png")

      md <- merged[!is.na(timbr_A) & !is.na(blup_A)]
      md[, sing_wild := vapply(top_series_pattern, function(p) {
        v <- tryCatch(parse_pat(p), error = function(e) rep(NA_integer_, 8))
        if (length(v) != 8 || any(is.na(v))) return(FALSE)
        tb <- table(v); any(tb[as.character(v[match(WILD, FOUNDERS)])] == 1)
      }, logical(1))]
      panel_df <- function(rows) {
        out <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
          bl <- as.numeric(rows[i, paste0("blup_", FOUNDERS), with = FALSE]); bl <- bl - mean(bl)
          tm <- as.numeric(rows[i, paste0("timbr_", FOUNDERS), with = FALSE]); tm <- tm - mean(tm)
          lab <- sprintf("%s (chr%s, %d alleles)", rows\$gene_name[i], rows\$chr[i], rows\$n_functional_alleles[i])
          rbind(data.table(locus = lab, strain = SNAME[FOUNDERS], model = "qtl2 BLUP", eff = bl),
                data.table(locus = lab, strain = SNAME[FOUNDERS], model = "TIMBR collapsed", eff = tm))
        }))
        out[, strain := factor(strain, levels = SNAME[FOUNDERS])]
        out[, locus := factor(locus, levels = unique(locus))]
        out
      }
      d_fig <- function(rows, ttl, file) {
        df <- panel_df(rows)
        p <- ggplot(df, aes(x = strain, y = eff, fill = model)) +
          geom_col(position = position_dodge(width = 0.8), width = 0.75) +
          geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
          scale_fill_manual(values = c("qtl2 BLUP" = "#4575b4", "TIMBR collapsed" = "#d73027"), name = NULL) +
          facet_wrap(~ locus, scales = "free_y") +
          labs(title = ttl, x = NULL, y = "Centered founder effect") +
          theme_bw(base_size = 11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "top")
        ggsave(file, p, width = 12, height = 8, dpi = 300); df
      }
      d1 <- md[order(-lod)][1:min(6, .N)]
      d2 <- md[sing_wild == TRUE][order(-lod)][1:min(6, .N)]
      d3 <- md[!is.na(n_functional_alleles)][, .SD[which.max(lod)], by = n_functional_alleles][order(n_functional_alleles)][1:min(8, .N)]
      d_fig(d1, "${study_prefix} - qtl2 BLUP vs TIMBR: highest-LOD eQTL", "${study_prefix}_blup_vs_timbr_topLOD.png")
      d_fig(d2, "${study_prefix} - qtl2 BLUP vs TIMBR: CAST/PWK resolved as singleton", "${study_prefix}_blup_vs_timbr_wildderived.png")
      d_fig(d3, "${study_prefix} - qtl2 BLUP vs TIMBR: spanning allele-group counts", "${study_prefix}_blup_vs_timbr_spanalleles.png")

      me <- merged[!is.na(timbr_A) & !is.na(blup_A)]
      Tcols <- paste0("timbr_", FOUNDERS)
      Bc <- as.matrix(me[, ..bcols]); Bc <- Bc - rowMeans(Bc)
      Tc <- as.matrix(me[, ..Tcols]); Tc <- Tc - rowMeans(Tc)
      econc <- rbindlist(lapply(seq_along(FOUNDERS), function(j)
        data.table(strain = SNAME[FOUNDERS[j]], blup = Bc[, j], timbr = Tc[, j])))
      econc[, strain := factor(strain, levels = SNAME[FOUNDERS])]
      rlab <- econc[, .(r = cor(blup, timbr, use = "complete.obs")), by = strain]
      overall_r <- econc[, cor(blup, timbr, use = "complete.obs")]
      rng <- range(c(econc\$blup, econc\$timbr), na.rm = TRUE)
      pE <- ggplot(econc, aes(x = blup, y = timbr, color = strain)) +
        geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.3) +
        geom_point(size = 0.25, alpha = 0.18) +
        geom_text(data = rlab, aes(x = rng[1], y = rng[2], label = sprintf("r = %.2f", r)),
                  hjust = 0, vjust = 1, size = 3, color = "black", inherit.aes = FALSE) +
        scale_color_manual(values = PAL_S, guide = "none") + facet_wrap(~ strain, nrow = 2) +
        labs(title = sprintf("${study_prefix} - qtl2 BLUP vs TIMBR concordance (overall r = %.2f)", overall_r),
             x = "qtl2 BLUP effect (mean-centered)", y = "TIMBR collapsed effect (mean-centered)") +
        theme_bw(base_size = 12)
      ggsave("${study_prefix}_blup_vs_timbr_concordance.png", pE, width = 10, height = 6, dpi = 300)

      w("--- Figure E: qtl2 BLUP vs TIMBR founder concordance (n = ", nrow(me), " loci) ---")
      w(sprintf("  overall r = %+.3f", overall_r))
      for (i in seq_len(nrow(rlab))) w(sprintf("  %-9s r = %+.3f", as.character(rlab\$strain[i]), rlab\$r[i]))
      w("")
      logm("Figure E (concordance) done: overall r = ", sprintf("%.3f", overall_r))
    }

    logm("All figures produced: A, B, C1, C2, D1-D3, E")

    close(sm_con); close(log_con)
    """
}
