# Data for Fig. 2b and 2c (gallstone disease as the primary outcome):
#   fig2m_surgery.csv  cholecystectomy endpoint, sleep -> MVPA at 15-60 min (RR and RD per 1,000)
#   fig2m_compare.csv  30-min sleep -> MVPA and sleep -> light RR for the five diseases and the two
#                      negative control outcomes, plus the within-replicate ratio RR(MVPA) / RR(light)
# The gallstone dose-response curves come from fig2_full_cohort.csv (data/fig2_data.R).
#
# Run: Rscript figures/data/fig2_main_data.R

suppressPackageStartupMessages(library(data.table))
PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))
OUT <- file.path(PROJ, "figures/data")

s <- summarise_dir(res_dir("gallstone_surgery", "primary_range_thr"))[cohort == "full" & from == "mvpa"]
fwrite(s[, .(direction, duration, minutes, rr, lower, upper, rd_per_1000, rd_lower, rd_upper, ref_per_1000, n_boot)],
       file.path(OUT, "fig2m_surgery.csv"))

outs <- c("gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs", "nco_cyst", "nco_earwax")
cmp <- rbindlist(lapply(outs, function(o) {
  d <- res_dir(o, "primary_range_thr")
  x <- summarise_dir(d)[cohort == "full" & duration == -30 & from %in% c("mvpa", "light")]
  b <- rbindlist(lapply(list.files(file.path(d, "boots"), full.names = TRUE), readRDS))
  b <- dcast(b[cohort == "full" & duration == -30 & from %in% c("mvpa", "light")], boot ~ from, value.var = "rr")
  p <- readRDS(file.path(d, "point_estimate.rds"))[cohort == "full" & duration == -30]
  r <- b$mvpa / b$light
  x[, `:=`(outcome = o, ratio = p[from == "mvpa", rr] / p[from == "light", rr],
           ratio_lower = quantile(r, 0.025), ratio_upper = quantile(r, 0.975))]
  x[, .(outcome, from, rr, lower, upper, ratio, ratio_lower, ratio_upper, n_boot)]
}))
fwrite(cmp, file.path(OUT, "fig2m_compare.csv"))
print(s[, .(minutes, rr = round(rr, 3), rd = round(rd_per_1000, 2), rd_lower = round(rd_lower, 2), rd_upper = round(rd_upper, 2))])
print(cmp[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
