# Absolute risk differences (cases per 1,000 over follow-up) and E-values for the
# primary isotemporal substitutions, plus E-values for the reference-composition contrasts.
#
# Usage: OUTCOME=<outcome> Rscript scripts/12_absolute_risk_evalues.R [result_dir=primary_range_thr]

suppressPackageStartupMessages(library(data.table))

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (Sys.getenv("OUTCOME", "gallstone") == "gallstone") file.path(proj, "results") else
  file.path(proj, "results", "outcomes", Sys.getenv("OUTCOME"))
res_dir <- file.path(res_root, if (length(args) >= 1) args[1] else "primary_range_thr")

# E-value (VanderWeele & Ding 2017); outcome is rare, so the risk ratio is used directly
evalue <- function(rr) {
  rr <- ifelse(rr < 1, 1 / rr, rr)
  rr + sqrt(rr * (rr - 1))
}
evalue_ci <- function(lower, upper) {
  # bound closest to the null; 1 if the interval crosses the null
  ifelse(lower <= 1 & upper >= 1, 1, ifelse(upper < 1, evalue(upper), evalue(lower)))
}

point <- readRDS(file.path(res_dir, "point_estimate.rds"))
boots <- rbindlist(lapply(list.files(file.path(res_dir, "boots"), full.names = TRUE), readRDS))
key <- c("cohort", "from", "duration")
boots[, rd := (sub_risk - ref_risk) * 1000]
ci <- boots[, .(rd_lower = quantile(rd, 0.025), rd_upper = quantile(rd, 0.975),
                rr_lower = quantile(rr, 0.025), rr_upper = quantile(rr, 0.975)), by = key]
out <- merge(point[, .(cohort, from, duration, ref_risk, sub_risk, rr, prop_sub)], ci, by = key)
out[, `:=`(ref_per_1000 = ref_risk * 1000, rd_per_1000 = (sub_risk - ref_risk) * 1000,
           e_value = evalue(rr), e_value_ci = evalue_ci(rr_lower, rr_upper))]
out[, direction := ifelse(duration > 0, paste(from, "-> sleep"), paste("sleep ->", from))]
out[, minutes := abs(duration)]
setorder(out, cohort, from, duration)
fwrite(out, file.path(res_dir, "absolute_risk_and_evalues.csv"))

cat("30-min substitutions: cases per 1,000 over ~9.8 years, and E-values\n")
print(out[minutes == 30 & cohort %in% c("short", "normal", "long", "full"),
          .(cohort, direction, ref = round(ref_per_1000, 1),
            rd = sprintf("%.1f (%.1f, %.1f)", rd_per_1000, rd_lower, rd_upper),
            rr = sprintf("%.2f (%.2f, %.2f)", rr, rr_lower, rr_upper),
            e_value = round(e_value, 2), e_value_ci = round(e_value_ci, 2))])

ideal <- file.path(res_root, "ideal_composition", "contrasts_vs_typical.csv")
if (file.exists(ideal)) {
  ic <- fread(ideal)
  ic[, `:=`(e_value = evalue(rr), e_value_ci = evalue_ci(lower, upper))]
  fwrite(ic, file.path(res_root, "ideal_composition", "contrasts_vs_typical_evalues.csv"))
  cat("\nReference compositions vs typical:\n")
  print(ic, digits = 3)
}
