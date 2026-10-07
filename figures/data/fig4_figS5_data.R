# Data for Fig. S5 and Fig. 4b (adiposity-related evidence). Writes CSVs to figures/data/ only.
#   Fig. S5: 30-min full-cohort log RR (sleep -> MVPA, MVPA -> sleep) vs log HR per SD BMI, WLS meta-regression
#   Fig. 4b: fat-depot HR per SD after imaging, five models per disease (Fig. 4a reads results/mri_* directly)
# Usage: Rscript figures/data/fig5_data.R

PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))
out <- file.path(PROJ, "figures/data")
diseases <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs")

# ---- Fig. S5: per-disease 30-min RR with bootstrap SE of log RR, and BMI gradient ----
est <- rbindlist(lapply(diseases, function(o) {
  x <- summarise_dir(res_dir(o, "primary_range_thr"))
  x[cohort == "full" & from == "mvpa" & minutes == 30,
    .(disease = o, direction, duration, rr, lower, upper, se_log, n_boot, prop_sub)]
}))
bmi <- fread(file.path(PROJ, "results/heterogeneity/bmi_hr_by_disease.csv"))
pts <- merge(est, bmi, by = "disease")
pts[, `:=`(log_rr = log(rr), w = 1 / se_log^2)]

# WLS fit exactly as metareg() in scripts/21_heterogeneity.R; intercept kept for drawing the line
fits <- pts[, {
  fit <- lm(log_rr ~ log_hr_bmi, weights = w)
  X <- cbind(1, log_hr_bmi)
  se_b <- sqrt(solve(t(X) %*% diag(w) %*% X)[2, 2])
  q <- sum(w * resid(fit)^2)
  phi <- max(1, q / (.N - 2))
  .(intercept = unname(coef(fit)[1]), slope_refit = unname(coef(fit)[2]), se_refit = se_b,
    q_resid_refit = q, p_disp_refit = 2 * pt(-abs(coef(fit)[2] / (se_b * sqrt(phi))), .N - 2))
}, by = .(direction, duration)]
mr <- fread(file.path(PROJ, "results/heterogeneity/metareg_on_bmi_gradient.csv"))[from == "mvpa"]
fits <- merge(fits, mr[, .(duration, slope, se, q_resid, p_dispersion)], by = "duration")
fits[, check_slope := abs(slope_refit - slope)]
print(fits)
stopifnot(all(fits$check_slope < 1e-8), all(abs(fits$p_disp_refit - fits$p_dispersion) < 1e-8))
fwrite(pts, file.path(out, "figS5_points.csv"))
fwrite(fits, file.path(out, "figS5_fits.csv"))

# ---- Fig. 4b: fat-depot Cox models ----
keep <- c("VAT (per SD)" = "VAT", "VAT, adjusted for BMI" = "VAT + BMI",
          "VAT, mutually adjusted for ASAT" = "VAT + ASAT",
          "Liver PDFF, log (per SD)" = "log PDFF", "Liver PDFF, adjusted for BMI" = "PDFF + BMI")
cox <- rbindlist(lapply(diseases, function(o) {
  f <- if (o == "gallstone") file.path(PROJ, "results/fat_depots_gallstone_cox.csv")
       else file.path(PROJ, "results/outcomes", o, "fat_depots_outcome_cox.csv")
  fread(f)[model %in% names(keep)][, `:=`(disease = o, label = keep[model])]
}))
cox[, `:=`(model_order = match(label, keep), disease_order = match(disease, diseases))]
setorder(cox, disease_order, model_order)
print(cox[, .(disease, label, n, events, hr = round(hr, 2), lower = round(lower, 2), upper = round(upper, 2))])
fwrite(cox, file.path(out, "fig4_fatcox.csv"))
