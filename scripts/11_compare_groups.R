# Summarise bootstrap results: risk ratios with 95% CI per sleep group, and
# between-group ratios of risk ratios (effect modification by sleep duration).
#
# Usage: OUTCOME=<outcome> Rscript scripts/11_compare_groups.R [result_dir=primary_range_thr]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
res_name <- if (length(args) >= 1) args[1] else "primary_range_thr"
proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
res_root <- if (Sys.getenv("OUTCOME", "gallstone") == "gallstone") file.path(proj, "results") else
  file.path(proj, "results", "outcomes", Sys.getenv("OUTCOME"))
res_dir <- file.path(res_root, res_name)

point <- readRDS(file.path(res_dir, "point_estimate.rds"))
boots <- rbindlist(lapply(list.files(file.path(res_dir, "boots"), full.names = TRUE), readRDS))
n_boot <- uniqueN(boots$boot)

label_from <- c(mvpa = "MVPA", light = "Light activity", inactivity = "Inactivity")
fmt <- function(est, lo, hi) sprintf("%.2f (%.2f, %.2f)", est, lo, hi)

# ---- risk ratios by sleep group ----
ci <- boots[, .(lower = quantile(rr, 0.025), upper = quantile(rr, 0.975)),
            by = .(cohort, from, duration)]
rr_tab <- merge(point[, .(cohort, from, duration, ref_risk, rr, prop_sub)], ci,
                by = c("cohort", "from", "duration"))
rr_tab[, direction := ifelse(duration > 0,
                             paste(label_from[from], "-> sleep"),
                             paste("Sleep ->", label_from[from]))]
rr_tab[, minutes := abs(duration)]
rr_tab[, estimate := ifelse(prop_sub >= 0.75, fmt(rr, lower, upper), "NR (<75% substitutable)")]
rr_wide <- dcast(rr_tab, direction + minutes ~ cohort, value.var = "estimate")
setcolorder(rr_wide, c("direction", "minutes", "short", "normal", "long", "full"))
fwrite(rr_tab, file.path(res_dir, "rr_long_format.csv"))
fwrite(rr_wide, file.path(res_dir, "table_rr_by_sleep_group.csv"))

# ---- ratio of risk ratios between groups (reference = normal sleepers) ----
ratio_of <- function(dt) {
  w <- dcast(dt, from + duration ~ cohort, value.var = "rr")
  w[, `:=`(short_vs_normal = short / normal,
           long_vs_normal = long / normal,
           long_vs_short = long / short,
           male_vs_female = full_male / full_female,
           male_vs_female_short = short_male / short_female,
           male_vs_female_normal = normal_male / normal_female,
           male_vs_female_long = long_male / long_female)]
  melt(w[, .(from, duration, short_vs_normal, long_vs_normal, long_vs_short, male_vs_female,
             male_vs_female_short, male_vs_female_normal, male_vs_female_long)],
       id.vars = c("from", "duration"), variable.name = "contrast", value.name = "ratio")
}
rp <- ratio_of(point)
rb <- boots[, ratio_of(.SD), by = boot]
rb_ci <- rb[, .(lower = quantile(ratio, 0.025), upper = quantile(ratio, 0.975),
                # two-sided bootstrap p: how often the ratio falls on the other side of 1
                p_boot = pmin(1, 2 * min(mean(ratio <= 1), mean(ratio >= 1)))),
            by = .(from, duration, contrast)]
ratio_tab <- merge(rp, rb_ci, by = c("from", "duration", "contrast"))
ratio_tab[, direction := ifelse(duration > 0,
                                paste(label_from[from], "-> sleep"),
                                paste("Sleep ->", label_from[from]))]
ratio_tab[, `:=`(minutes = abs(duration), estimate = fmt(ratio, lower, upper))]
setorder(ratio_tab, direction, minutes, contrast)
fwrite(ratio_tab[, .(direction, minutes, contrast, ratio, lower, upper, p_boot, estimate)],
       file.path(res_dir, "table_ratio_between_groups.csv"))

cat(sprintf("Bootstrap replicates: %d\n\n30-min substitutions, RR (95%% CI):\n", n_boot))
print(rr_wide[minutes == 30], row.names = FALSE)
cat("\nBetween-group ratio of RRs, 30 min (CI excluding 1 = evidence of effect modification):\n")
print(ratio_tab[minutes == 30, .(direction, contrast, estimate, p_boot = round(p_boot, 3))],
      row.names = FALSE)
