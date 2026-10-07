# Shared helpers for the figure and table data scripts: summarise a g-formula result directory
# (point_estimate.rds + boots/*.rds) into one table with 95% percentile bootstrap CIs.
#
# Usage inside a data script that has set PROJ (the repository root):
#   source(file.path(PROJ, "figures/data/extract_helpers.R"))
#   x <- summarise_dir(res_dir("gerd", "primary_range_thr"))

suppressPackageStartupMessages(library(data.table))

# result directory for an outcome ("gallstone" lives directly under results/)
res_dir <- function(outcome, analysis) {
  root <- if (outcome == "gallstone") file.path(PROJ, "results") else file.path(PROJ, "results", "outcomes", outcome)
  file.path(root, analysis)
}

# duration < 0: sleep replaced by `from` (sleep -> X); duration > 0: X -> sleep
summarise_dir <- function(dir) {
  p <- readRDS(file.path(dir, "point_estimate.rds"))
  b <- rbindlist(lapply(list.files(file.path(dir, "boots"), full.names = TRUE), readRDS))
  key <- c("cohort", "from", "duration")
  b[, rd := (sub_risk - ref_risk) * 1000]
  ci <- b[, .(lower = quantile(rr, 0.025), upper = quantile(rr, 0.975),
              rd_lower = quantile(rd, 0.025), rd_upper = quantile(rd, 0.975),
              se_log = sd(log(rr)), n_boot = uniqueN(boot)), by = key]
  x <- merge(p[, .(cohort, from, duration, ref_risk, sub_risk, rr, prop_sub)], ci, by = key)
  x[, `:=`(ref_per_1000 = ref_risk * 1000, rd_per_1000 = (sub_risk - ref_risk) * 1000,
           direction = ifelse(duration < 0, paste0("sleep_to_", from), paste0(from, "_to_sleep")),
           minutes = abs(duration), reported = prop_sub >= 0.75)]
  setorder(x, cohort, from, duration)
  x[]
}
