# Data for Fig. 2a and Fig. S4 (RR, RD and dose-response across the five diseases; primary analyses, full cohort).
# Reads results/primary_range_thr (gallstone) and results/outcomes/<o>/primary_range_thr (others),
# read-only, and writes figures/data/fig2_full_cohort.csv.
# Run: Rscript figures/data/fig2_data.R

PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))

diseases <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs")

out <- rbindlist(lapply(diseases, function(o) {
  x <- summarise_dir(res_dir(o, "primary_range_thr"))
  x <- x[cohort == "full"]
  x[, disease := o]
  x
}))

setcolorder(out, c("disease", "cohort", "direction", "from", "duration", "minutes",
                   "ref_per_1000", "rr", "lower", "upper", "rd_per_1000", "rd_lower", "rd_upper",
                   "se_log", "prop_sub", "reported", "n_boot"))
out <- out[, .(disease, cohort, direction, from, duration, minutes, ref_per_1000, rr, lower, upper,
               rd_per_1000, rd_lower, rd_upper, se_log, prop_sub, reported, n_boot)]

dir.create(file.path(PROJ, "figures/data"), showWarnings = FALSE, recursive = TRUE)
fwrite(out, file.path(PROJ, "figures/data/fig2_full_cohort.csv"))

# console check (30 min, full cohort)
chk <- out[minutes == 30, .(disease, direction, ref = round(ref_per_1000, 1),
                            RR = sprintf("%.2f (%.2f, %.2f)", rr, lower, upper),
                            RD = sprintf("%.1f (%.1f, %.1f)", rd_per_1000, rd_lower, rd_upper),
                            n_boot, prop_sub = round(prop_sub, 3))]
print(chk, nrows = 100)
cat("rows:", nrow(out), " min prop_sub:", round(min(out$prop_sub), 3),
    " any not reported:", any(!out$reported), "\n")
