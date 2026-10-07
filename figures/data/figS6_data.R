# Data for Fig. S6 (sleep->MVPA 30 min/day by accelerometer-measured sleep-duration group:
# full-cohort model vs separate models per sleep-duration group).
# Reads, read-only:
#   full-cohort model : results/primary_range_thr, results/outcomes/<o>/primary_range_thr
#                       (five diseases, 500 replicates; two negative controls, 200 replicates)
#   separate models   : results/robust_strat, results/outcomes/<o>/robust_strat (200 replicates)
#   main-effects model: results/robust_reduced (gallstone only, 200 replicates)
# Writes figures/data/figS6_groups.csv, figS6_ratios.csv and figS6_heterogeneity.csv.
# Run: Rscript figures/data/fig4_data.R

PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))

diseases <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs")
ncos     <- c("nco_cyst", "nco_earwax")
groups   <- c("short", "normal", "long")

runs <- rbind(
  data.table(disease = c(diseases, ncos), model = "full_cohort", analysis = "primary_range_thr"),
  data.table(disease = diseases,          model = "separate",    analysis = "robust_strat"),
  data.table(disease = "gallstone",       model = "main_effects", analysis = "robust_reduced")
)

# ---- group-specific RR (95% percentile CI), sleep -> MVPA, 30 min ----------------------------
grp <- rbindlist(lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i]
  x <- summarise_dir(res_dir(r$disease, r$analysis))
  x <- x[cohort %in% groups & from == "mvpa" & duration == -30]
  x[, `:=`(disease = r$disease, model = r$model)]
  x[, .(model, disease, cohort, direction, minutes, rr, lower, upper, prop_sub, reported, n_boot)]
}))

# ---- between-group ratios, computed within each bootstrap replicate -----------------------
p_boot <- function(r) min(1, 2 * min(mean(r <= 1), mean(r >= 1)))
ratio_one <- function(dir, num, den) {
  p <- readRDS(file.path(dir, "point_estimate.rds"))[from == "mvpa" & duration == -30]
  b <- rbindlist(lapply(list.files(file.path(dir, "boots"), full.names = TRUE), readRDS))
  b <- b[from == "mvpa" & duration == -30 & cohort %in% c(num, den)]
  w <- dcast(b, boot ~ cohort, value.var = "rr")
  stopifnot(!anyNA(w[[num]]), !anyNA(w[[den]]))          # every replicate has both groups (paired)
  r <- w[[num]] / w[[den]]
  data.table(ratio = p[cohort == num, rr] / p[cohort == den, rr],
             lower = unname(quantile(r, 0.025)), upper = unname(quantile(r, 0.975)),
             p_boot = p_boot(r), se_log = sd(log(r)), n_boot = length(r))
}

contrasts <- list(long_vs_normal = c("long", "normal"), long_vs_short = c("long", "short"))
rat <- rbindlist(lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i]
  rbindlist(lapply(names(contrasts), function(k) {
    cbind(data.table(model = r$model, disease = r$disease, contrast = k),
          ratio_one(res_dir(r$disease, r$analysis), contrasts[[k]][1], contrasts[[k]][2]))
  }))
}))

# ---- Cochran Q across the five diseases, full-cohort long-vs-normal ratio (log scale) -----
h <- rat[model == "full_cohort" & contrast == "long_vs_normal" & disease %in% diseases]
w <- 1 / h$se_log^2; y <- log(h$ratio)
Q <- sum(w * (y - sum(w * y) / sum(w))^2)
het <- data.table(model = "full_cohort", contrast = "long_vs_normal", direction = "sleep_to_mvpa",
                  minutes = 30, Q = Q, df = nrow(h) - 1, p = pchisq(Q, nrow(h) - 1, lower.tail = FALSE))

out <- file.path(PROJ, "figures/data")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fwrite(grp, file.path(out, "figS6_groups.csv"))
fwrite(rat, file.path(out, "figS6_ratios.csv"))
fwrite(het, file.path(out, "figS6_heterogeneity.csv"))

# ---- console checks ------------------------------------------------------------------------
f <- function(a, l, u) sprintf("%.2f (%.2f-%.2f)", a, l, u)
print(grp[, .(model, disease, cohort, RR = f(rr, lower, upper), prop_sub = round(prop_sub, 3), n_boot)], nrows = 200)
print(rat[, .(model, disease, contrast, ratio = f(ratio, lower, upper), p_boot, n_boot)], nrows = 200)
print(het)

# agreement with the stored long_vs_normal / long_vs_short tables of the full-cohort runs
chk <- rbindlist(lapply(c(diseases, ncos), function(o) {
  s <- fread(file.path(res_dir(o, "primary_range_thr"), "table_ratio_between_groups.csv"))
  s <- s[direction == "Sleep -> MVPA" & minutes == 30 & contrast %in% names(contrasts)]
  m <- merge(s, rat[model == "full_cohort" & disease == o], by = "contrast", suffixes = c(".file", ".new"))
  m[, .(disease = o, contrast, d_ratio = abs(ratio.file - ratio.new), d_lo = abs(lower.file - lower.new),
        d_up = abs(upper.file - upper.new), d_p = abs(p_boot.file - p_boot.new))]
}))
print(chk)
cat("max abs diff vs table_ratio_between_groups.csv:", max(unlist(chk[, .(d_ratio, d_lo, d_up, d_p)])), "\n")

