#!/usr/bin/env Rscript
# Model-level significance thresholds from the screen permutations, instead of each phenotype's own.
#
# Module 7 permutes every phenotype that passes the LOD pre-filter on its own: 50 screen permutations, then 1,000
# for the screen-passers. For a whole transcriptome that takes weeks of cluster time. When every phenotype is in
# normal scores (rank-based inverse normal, as Metadata_phenotype_QTL2_NF.Rmd makes the expression values), all
# phenotypes have nearly the same null distribution of the genome-wide maximum LOD, so the screen permutations already
# finished can be pooled into one set of thresholds for the model. This script:
#   1. pools every finished screen batch of a run (<prefix>_*_batch_1.rds: 50 permutations of one phenotype);
#   2. checks the assumption: the spread of the phenotypes' own 95% quantiles against the spread expected from
#      drawing the same number of values from the pooled null. A ratio near 1 means the phenotypes do not differ
#      beyond sampling noise; a ratio well above 1 means per-phenotype thresholds are needed;
#   3. writes a new input folder for main_resume.nf --resume_from significant_qtls: links to the run's steps 01 to 06
#      and its filtered cross2 object and phenotype list, <prefix>_permThresh.rds with the pooled 63, 90, 95 and 99%
#      thresholds for every phenotype that passed the pre-filter, <prefix>_fullPerm.rds (the pooled values),
#      pooled_thresholds.tsv, pooled_threshold_spread.tsv and pooled_screen_batches.txt.
# Stop the permutation run first (see "Thresholds for a whole model" in README.md); running batch tasks finish and
# write their files by rename, so a half-written batch is never read.
#
#   Rscript bin/pooled_thresholds.R RUN_OUTDIR NEW_INPUT_DIR STUDY_PREFIX [MIN_PHENOTYPES, default 300]
# then
#   nextflow run main_resume.nf <the run's settings> --resume_from significant_qtls --input_dir NEW_INPUT_DIR \
#     --outdir <a new results folder>
args <- commandArgs(TRUE)
if (length(args) < 3) stop("usage: pooled_thresholds.R RUN_OUTDIR NEW_INPUT_DIR STUDY_PREFIX [MIN_PHENOTYPES]")
run <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
prefix <- args[3]
min_n <- if (length(args) >= 4) as.integer(args[4]) else 300L
if (file.exists(out)) stop(out, " exists; give a new folder")
pd <- file.path(run, "07_permutation_testing")
files <- sort(list.files(file.path(pd, "batches_screen"), pattern = paste0("^", prefix, "_.*_batch_1\\.rds$"),
                         full.names = TRUE))
if (length(files) < min_n) stop(length(files), " finished screen batches; at least ", min_n, " are needed")
per <- lapply(files, function(f) as.numeric(readRDS(f)[, 1]))
perms <- unlist(per)
if (!all(is.finite(perms))) stop("non-finite permutation LODs in the screen batches")
levels <- c("63%", "90%", "95%", "99%")
q <- quantile(perms, c(0.63, 0.90, 0.95, 0.99), names = FALSE)

set.seed(1)
k <- as.integer(stats::median(lengths(per)))
obs <- vapply(per, function(x) quantile(x, 0.95, names = FALSE), 0)
sim <- replicate(length(per), quantile(sample(perms, k), 0.95, names = FALSE))
ratio <- stats::sd(obs) / stats::sd(sim)

genes <- readLines(file.path(pd, paste0(prefix, "_phenotype_list.txt")))
genes <- genes[nzchar(trimws(genes))]
od <- file.path(out, "07_permutation_testing")
dir.create(od, recursive = TRUE)
for (s in c("01_phenotype_processing", "02_genotype_processing", "03_control_file_generation", "04_cross2_creation",
            "05_genome_scan_preparation", "06_qtl_analysis")) {
  if (!dir.exists(file.path(run, s))) stop(run, " has no ", s)
  file.symlink(file.path(run, s), file.path(out, s))
}
for (f in paste0(prefix, c("_filtered_cross2.rds", "_phenotype_list.txt"))) {
  if (!file.exists(file.path(pd, f))) stop(pd, " has no ", f)
  file.symlink(file.path(pd, f), file.path(od, f))
}
th <- matrix(q, nrow = 4, ncol = length(genes), dimnames = list(levels, genes))
saveRDS(th, file.path(od, paste0(prefix, "_permThresh.rds")))
saveRDS(matrix(perms, ncol = 1, dimnames = list(NULL, "pooled_screen_permutations")),
        file.path(od, paste0(prefix, "_fullPerm.rds")))
writeLines(basename(files), file.path(od, "pooled_screen_batches.txt"))
write.table(data.frame(phenotypes_pooled = length(files), permutations = length(perms), level = levels,
                       lod = round(q, 3)),
            file.path(od, "pooled_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(phenotypes = length(per), permutations_each = k, sd_own_q95 = round(stats::sd(obs), 4),
                       sd_from_sampling = round(stats::sd(sim), 4), ratio = round(ratio, 3),
                       own_q95_1st_pct = round(quantile(obs, 0.01, names = FALSE), 3),
                       own_q95_99th_pct = round(quantile(obs, 0.99, names = FALSE), 3)),
            file.path(od, "pooled_threshold_spread.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d phenotypes pooled, %d permutations; thresholds 63/90/95/99%% = %s; for %d phenotypes\n",
            length(files), length(perms), paste(round(q, 2), collapse = " / "), length(genes)))
cat(sprintf("own 95%% quantiles: sd %.3f against %.3f from sampling alone (ratio %.2f)%s\n", stats::sd(obs),
            stats::sd(sim), ratio, if (ratio > 1.15) "; phenotypes differ beyond sampling, consider per-phenotype thresholds" else ""))
