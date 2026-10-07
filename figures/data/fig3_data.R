# Data for Fig. 3 (robustness forest): 30-min sleep -> MVPA RR, full cohort, across the primary
# analysis, the three core sensitivity analyses, the disease-specific checks and the two negative
# control outcomes. Reads point_estimate.rds + boots/*.rds through extract_helpers.R and writes
# figures/data/fig3_robustness.csv. Read-only with respect to results/ and data/.
#
# Run: Rscript figures/data/fig3_data.R

suppressPackageStartupMessages(library(data.table))
PROJ <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
source(file.path(PROJ, "figures/data/extract_helpers.R"))

# One row per plotted estimate, in plotting order.
#   group  : disease block in the figure (negative controls form the last block)
#   tier   : primary / core (lag3, s2, healthy-user) / specific (disease-specific check) / nco
#   N, events : analytic N and events of each analysis (as reported in Fig. 3 and Table 2)
rows <- fread(text = "
group|outcome|analysis|tier|key|N|events
gallstone|gallstone|primary_range_thr|primary|primary|86596|1961
gallstone|gallstone|lag3_range_thr|core|lag3|85893|1258
gallstone|gallstone|s2_range_thr|core|s2|86596|1961
gallstone|gallstone|robust_healthuser|core|healthuser|86596|1961
gallstone|gallstone_surgery|primary_range_thr|specific|cholecystectomy|86596|1042
gallstone|gallstone|s3_range_thr|specific|s3|86596|1961
gallstone|gallstone|primary_range|specific|ml|86596|1961
gallstone|gallstone|robust_knots4|specific|knots4|86596|1961
gallstone|gallstone|robust_reduced|specific|reduced|86596|1961
gallstone|gallstone|robust_mrisub|specific|mrisub|19752|374
peptic_ulcer|peptic_ulcer|primary_range_thr|primary|primary|89540|983
peptic_ulcer|peptic_ulcer|lag3_range_thr|core|lag3|89165|608
peptic_ulcer|peptic_ulcer|s2_range_thr|core|s2|89540|983
peptic_ulcer|peptic_ulcer|robust_healthuser|core|healthuser|89540|983
peptic_ulcer|peptic_ulcer|robust_noacid|specific|noacid|81843|812
peptic_ulcer|peptic_ulcer|robust_nsaid|specific|nsaid|89540|983
gerd|gerd|primary_range_thr|primary|primary|82983|5577
gerd|gerd|lag3_range_thr|core|lag3|81022|3616
gerd|gerd|s2_range_thr|core|s2|82983|5577
gerd|gerd|robust_healthuser|core|healthuser|82983|5577
gerd|gerd|robust_noacid|specific|noacid|78935|4761
diverticular|diverticular|primary_range_thr|primary|primary|86008|7167
diverticular|diverticular|lag3_range_thr|core|lag3|83578|4737
diverticular|diverticular|s2_range_thr|core|s2|86008|7167
diverticular|diverticular|robust_healthuser|core|healthuser|86008|7167
ibs|ibs|primary_range_thr|primary|primary|83646|1031
ibs|ibs|lag3_range_thr|core|lag3|83257|642
ibs|ibs|s2_range_thr|core|s2|83646|1031
ibs|ibs|robust_healthuser|core|healthuser|83646|1031
nco|nco_cyst|primary_range_thr|nco|nco_cyst|89481|646
nco|nco_earwax|primary_range_thr|nco|nco_earwax|84632|577
", sep = "|")

out <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
  r <- rows[i]
  d <- res_dir(r$outcome, r$analysis)
  x <- summarise_dir(d)
  x <- x[cohort == "full" & from == "mvpa" & duration == -30]
  stopifnot(nrow(x) == 1)
  cbind(r, x[, .(rr, lower, upper, rd_per_1000, rd_lower, rd_upper, ref_per_1000,
                 prop_sub, n_boot, reported)],
        dir = sub(paste0(PROJ, "/"), "", d, fixed = TRUE))
}))
out[, order := .I]

stopifnot(all(out$reported))
fwrite(out, file.path(PROJ, "figures/data/fig3_robustness.csv"))
print(out[, .(group, key, N, events, rr = round(rr, 4), lower = round(lower, 4),
              upper = round(upper, 4), n_boot, prop_sub = round(prop_sub, 3))])
