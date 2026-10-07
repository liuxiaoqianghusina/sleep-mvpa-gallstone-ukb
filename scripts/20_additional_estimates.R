# Additional estimates reported in the paper, computed from the saved point estimates and bootstrap replicates
# (no refitting); sleep -> MVPA, 30 min/day, full cohort:
#   negative_control_ratios.csv  RR of each disease / RR of each negative control outcome. The two sets of replicates
#                                are independent and are paired by replicate number, so the CIs are conservative.
#   mvpa_vs_light.csv            RR(sleep -> MVPA) / RR(sleep -> light activity) within the same replicate
#   attenuation.csv              relative change in log RR after each core sensitivity analysis
#   ci99.csv                     99% bootstrap percentile CIs of the primary estimates
#   feasibility.csv              range of the proportion of participants able to take the reallocations
# Output: results/additional_estimates/
# Usage:  Rscript scripts/20_additional_estimates.R

suppressPackageStartupMessages(library(data.table))
PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))
OUT <- file.path(PROJ, "results", "additional_estimates")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DIS <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs")

boots <- function(d) rbindlist(lapply(list.files(file.path(d, "boots"), full.names = TRUE), readRDS))
point <- function(d) readRDS(file.path(d, "point_estimate.rds"))
p_boot <- function(r) min(1, 2 * min(mean(r <= 1), mean(r >= 1)))
get30 <- function(o, a, fr = "mvpa", dur = -30) {
  d <- res_dir(o, a)
  list(p = point(d)[cohort == "full" & from == fr & duration == dur, rr],
       b = boots(d)[cohort == "full" & from == fr & duration == dur, .(boot, rr)])
}

# 1. disease / negative control
nco <- rbindlist(lapply(c("nco_cyst", "nco_earwax"), function(nc) {
  n <- get30(nc, "primary_range_thr")
  rbindlist(lapply(DIS, function(o) {
    x <- get30(o, "primary_range_thr")
    m <- merge(x$b, n$b, by = "boot")
    r <- m$rr.x / m$rr.y
    data.table(disease = o, negative_control = nc, ratio = x$p / n$p, lower = quantile(r, 0.025),
               upper = quantile(r, 0.975), p_boot = p_boot(r), n_pairs = nrow(m))
  }))
}))
fwrite(nco, file.path(OUT, "negative_control_ratios.csv"))

# 2. MVPA vs light activity within the same replicate
ml <- rbindlist(lapply(DIS, function(o) {
  m <- get30(o, "primary_range_thr")
  l <- get30(o, "primary_range_thr", "light")
  b <- merge(m$b, l$b, by = "boot")
  r <- b$rr.x / b$rr.y
  data.table(disease = o, ratio = m$p / l$p, lower = quantile(r, 0.025), upper = quantile(r, 0.975),
             p_boot = p_boot(r))
}))
fwrite(ml, file.path(OUT, "mvpa_vs_light.csv"))

# 3. attenuation of log RR relative to the primary analysis
att <- rbindlist(lapply(DIS, function(o) {
  p <- get30(o, "primary_range_thr")$p
  rbindlist(lapply(c("s2_range_thr", "robust_healthuser", "lag3_range_thr"), function(a) {
    q <- get30(o, a)$p
    data.table(disease = o, analysis = a, rr = q, rr_primary = p, attenuation_pct = 100 * (1 - log(q) / log(p)))
  }))
}))
fwrite(att, file.path(OUT, "attenuation.csv"))

# 4. 99% CIs
ci99 <- rbindlist(lapply(DIS, function(o) {
  x <- get30(o, "primary_range_thr")
  data.table(disease = o, rr = x$p, lower_99 = quantile(x$b$rr, 0.005), upper_99 = quantile(x$b$rr, 0.995))
}))
fwrite(ci99, file.path(OUT, "ci99.csv"))

# 5. proportion able to take the reallocations
feas <- rbindlist(lapply(DIS, function(o) {
  p <- point(res_dir(o, "primary_range_thr"))
  data.table(disease = o, full_min = min(p[cohort == "full", prop_sub]), full_max = max(p[cohort == "full", prop_sub]),
             full_sleep_to_mvpa_30 = p[cohort == "full" & from == "mvpa" & duration == -30, prop_sub],
             all_groups_min = min(p$prop_sub))
}))
fwrite(feas, file.path(OUT, "feasibility.csv"))

for (x in list(nco, ml, att, ci99, feas)) print(x[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
