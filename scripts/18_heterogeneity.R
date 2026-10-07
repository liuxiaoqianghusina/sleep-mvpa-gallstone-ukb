# Between-disease heterogeneity of the isotemporal reallocation effects.
#
#   1. Cochran's Q / I^2 across the five diseases for every substitution and cohort,
#      on log RR with bootstrap standard errors.
#   2. Consistency of effect modification by sleep duration: Q across diseases of the
#      long-vs-normal and short-vs-normal ratios of RRs.
#   3. Pairwise ratios of RRs between diseases (Holm-adjusted bootstrap p).
#   4. Meta-regression of log RR on each disease's adiposity gradient (log HR per SD BMI,
#      age- and sex-adjusted Cox in the same analytic cohort), testing whether the
#      substitution effect scales with how "metabolic" the disease is.
#
# Bootstrap samples were drawn separately per disease, so estimates are treated as
# independent. The cohorts overlap and the true covariances are positive, so ignoring
# them overstates Var(difference): Q and the pairwise tests are conservative.
#
# Usage: Rscript scripts/18_heterogeneity.R

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))

diseases <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs")
res_dir <- function(o) file.path(if (o == "gallstone") file.path(proj, "results")
                                 else file.path(proj, "results", "outcomes", o), "primary_range_thr")
out_dir <- file.path(proj, "results", "heterogeneity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

key <- c("cohort", "from", "to", "duration")
boots <- rbindlist(lapply(diseases, function(o) {
  b <- rbindlist(lapply(list.files(file.path(res_dir(o), "boots"), full.names = TRUE), readRDS))
  b[, disease := o]
}))
point <- rbindlist(lapply(diseases, function(o) readRDS(file.path(res_dir(o), "point_estimate.rds"))[, disease := o]))
point[, to := NULL]  # point files carry `to` under the same substitution grid
point <- merge(point, unique(boots[, .(from, to, duration)]), by = c("from", "duration"))

# ---- 1. heterogeneity of RRs ----
se <- boots[, .(se = sd(log(rr))), by = c(key, "disease")]
est <- merge(point[, c(key, "disease", "rr", "prop_sub"), with = FALSE], se, by = c(key, "disease"))
cochran <- function(y, s) {
  w <- 1 / s^2
  mu <- sum(w * y) / sum(w)
  q <- sum(w * (y - mu)^2)
  df <- length(y) - 1
  list(q = q, df = df, p = pchisq(q, df, lower.tail = FALSE), i2 = max(0, (q - df) / q))
}
het <- est[, c(cochran(log(rr), se), n_reportable = sum(prop_sub >= 0.75)), by = key]
fwrite(het, file.path(out_dir, "q_rr_across_diseases.csv"))

# ---- 2. consistency of effect modification by sleep duration ----
ratio_long <- function(dt, a, b) {
  w <- dcast(dt, disease + from + to + duration + boot ~ cohort, value.var = "rr")
  w[, .(disease, from, to, duration, boot, ratio = get(a) / get(b))]
}
em <- rbindlist(lapply(list(c("long", "normal"), c("short", "normal")), function(p) {
  rb <- ratio_long(boots[cohort %in% p], p[1], p[2])
  rp <- ratio_long(copy(point)[, boot := 0L][cohort %in% p], p[1], p[2])
  s <- rb[, .(se = sd(log(ratio)), lower = quantile(ratio, 0.025), upper = quantile(ratio, 0.975)),
          by = .(disease, from, to, duration)]
  merge(rp[, !"boot"], s, by = c("disease", "from", "to", "duration"))[, contrast := paste(p, collapse = "_vs_")]
}))
fwrite(em, file.path(out_dir, "sleep_group_ratio_by_disease.csv"))
em_het <- em[, c(cochran(log(ratio), se),
                 n_same_direction = max(sum(ratio < 1), sum(ratio > 1)),
                 n_ci_excl_1 = sum(lower > 1 | upper < 1)),
             by = .(contrast, from, to, duration)]
fwrite(em_het, file.path(out_dir, "q_sleep_group_modification.csv"))

# ---- 3. pairwise ratios of RRs (independent bootstrap pairing) ----
pairs <- CJ(a = diseases, b = diseases)[match(a, diseases) < match(b, diseases)]
pw <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$a[i]; b <- pairs$b[i]
  ba <- boots[disease == a & cohort == "full"]; bb <- boots[disease == b & cohort == "full"]
  m <- merge(ba[, .(from, to, duration, boot, ra = rr)], bb[, .(from, to, duration, boot, rb = rr)],
             by = c("from", "to", "duration", "boot"))
  m[, ratio := ra / rb]
  pt <- merge(point[disease == a & cohort == "full", .(from, to, duration, ra = rr)],
              point[disease == b & cohort == "full", .(from, to, duration, rb = rr)],
              by = c("from", "to", "duration"))
  s <- m[, .(lower = quantile(ratio, 0.025), upper = quantile(ratio, 0.975),
             p_boot = pmin(1, 2 * min(mean(ratio <= 1), mean(ratio >= 1)))), by = .(from, to, duration)]
  merge(pt[, .(from, to, duration, ratio = ra / rb)], s, by = c("from", "to", "duration"))[, `:=`(disease_a = a, disease_b = b)]
}))
pw[, p_holm := p.adjust(p_boot, "holm"), by = .(from, to, duration)]
fwrite(pw, file.path(out_dir, "pairwise_ratio_full_cohort.csv"))

# ---- 4. adiposity gradient per disease and meta-regression ----
bmi_hr <- rbindlist(lapply(diseases, function(o) {
  Sys.setenv(OUTCOME = o)
  d <- load_cohort(proj)$df[!is.na(bmi)]
  d[, bmi_z := (bmi - mean(bmi)) / sd(bmi)]
  m <- coxph(Surv(time_to_gs, gs) ~ bmi_z + age_accel + sex, data = d)
  data.table(disease = o, n = m$n, events = m$nevent,
             log_hr_bmi = unname(coef(m)["bmi_z"]), se_log_hr_bmi = sqrt(vcov(m)["bmi_z", "bmi_z"]))
}))
fwrite(bmi_hr[, hr_bmi_per_sd := exp(log_hr_bmi)], file.path(out_dir, "bmi_hr_by_disease.csv"))

metareg <- function(y, s, x) {
  w <- 1 / s^2
  fit <- lm(y ~ x, weights = w)
  b <- coef(fit)[2]
  se_b <- sqrt(solve(t(cbind(1, x)) %*% diag(w) %*% cbind(1, x))[2, 2])  # fixed-effect SE
  q_resid <- sum(w * resid(fit)^2)
  # with only five diseases, also allow for residual heterogeneity: inflate the SE by the
  # dispersion (Q_resid / df, floored at 1) and use a t reference on k - 2 df
  phi <- max(1, q_resid / (length(y) - 2))
  list(slope = unname(b), se = se_b, p = 2 * pnorm(-abs(b / se_b)), q_resid = q_resid,
       p_dispersion = 2 * pt(-abs(b / (se_b * sqrt(phi))), length(y) - 2))
}
mr <- merge(est[cohort == "full"], bmi_hr, by = "disease")[
  , metareg(log(rr), se, log_hr_bmi), by = .(from, to, duration)]
fwrite(mr, file.path(out_dir, "metareg_on_bmi_gradient.csv"))

# ---- print the 30-min summary ----
lab <- function(dt) dt[, direction := ifelse(duration > 0, paste(from, "-> sleep"), paste("sleep ->", from))]
cat("Cochran's Q across 5 diseases, full cohort, 30-min substitutions:\n")
print(lab(het[cohort == "full" & abs(duration) == 30])[
  , .(direction, Q = round(q, 1), p = signif(p, 2), I2 = round(i2, 2))], row.names = FALSE)
cat("\nPer-disease RR (full cohort, 30 min):\n")
print(dcast(lab(est[cohort == "full" & abs(duration) == 30])[, rr := round(rr, 3)], direction ~ disease,
            value.var = "rr"), row.names = FALSE)
cat("\nSleep-group effect modification across diseases (30 min):\n")
print(lab(em_het[abs(duration) == 30])[
  , .(contrast, direction, Q = round(q, 1), p = signif(p, 2), I2 = round(i2, 2),
      same_dir = n_same_direction, ci_excl_1 = n_ci_excl_1)], row.names = FALSE)
cat("\nBMI gradient per disease (HR per SD):\n")
print(bmi_hr[, .(disease, events, hr = round(hr_bmi_per_sd, 2))], row.names = FALSE)
cat("\nMeta-regression of log RR on log HR(BMI), 30 min:\n")
print(lab(mr[abs(duration) == 30])[, .(direction, slope = round(slope, 3), p = signif(p, 2), p_dispersion = signif(p_dispersion, 2))], row.names = FALSE)
cat("\nPairwise, sleep -> MVPA / sleep -> light, 30 min (Holm p):\n")
print(pw[duration == -30 & from %in% c("mvpa", "light"),
         .(from, pair = paste(disease_a, "/", disease_b),
           est = sprintf("%.2f (%.2f, %.2f)", ratio, lower, upper), p_holm = round(p_holm, 3))], row.names = FALSE)
