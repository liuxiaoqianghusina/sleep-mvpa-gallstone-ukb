# Data for Additional file 1, Tables S3-S5, from the primary analyses of the five diseases:
#   summaries_primary.csv     every cohort x direction x duration: RR, RD per 1,000, 95% CI and the proportion of
#                             participants able to take the reallocation (Tables S3 and S4)
#   ratios_between_groups.csv sleep-group and sex ratios of RRs with paired bootstrap P values (Table S5, with
#                             figures/data/figS6_*.csv for the separate models by sleep-duration group)
# Table S6 uses results/mri_*/mean_differences.csv and figures/data/fig4_fatcox.csv directly.
# Outputs: tables/output/summaries_primary.csv, tables/output/ratios_between_groups.csv
# Usage:   Rscript tables/supp_tables_data.R

suppressPackageStartupMessages(library(data.table))
PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))
OUT <- file.path(PROJ, "tables", "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dis <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs")

s <- rbindlist(lapply(dis, function(o) {
  x <- summarise_dir(res_dir(o, "primary_range_thr"))
  x[, outcome := o]
}), fill = TRUE)
fwrite(s[, .(outcome, cohort, from, direction, duration, minutes, prop_sub, reported, ref_per_1000, rr, lower, upper,
             rd_per_1000, rd_lower, rd_upper, n_boot)], file.path(OUT, "summaries_primary.csv"))

r <- rbindlist(lapply(dis, function(o) {
  x <- fread(file.path(res_dir(o, "primary_range_thr"), "table_ratio_between_groups.csv"))
  x[, outcome := o]
}))
fwrite(r, file.path(OUT, "ratios_between_groups.csv"))

cat("cohorts:", paste(sort(unique(s$cohort)), collapse = " "), "\n")
cat("n_boot:", paste(sort(unique(s$n_boot)), collapse = " "), "\n")
print(s[cohort == "full" & from == "mvpa" & duration == -30,
        .(outcome, rr = round(rr, 2), lower = round(lower, 2), upper = round(upper, 2), rd = round(rd_per_1000, 1),
          prop = round(prop_sub, 4))])
